# command line interface ------------------------------------------------------#
library(docopt)
library(tidyverse)
'Validate putatEVEs based on retroblast hits

Usage:
  validate.R [options] <in.bed> <in.fai> <out.tsv> <out.pdf>
  
Options
  -p --prefix=<string>          EVE ID prefix in output files.
  -b --min-bitscore-frac=<0:1>  Minimum bitscore relative to top hit per locus 
                                to include hit in validation [default: 0.5]
  -E --eve-score-high=<0:100>   Minimum eve-score for high-confidence validatEVEs
                                [default: 30]
  -e --eve-score-low=<0:100>    Minimum eve-score for low-confidence validatEVEs
                                [default: 10]
  -r --retro-score-low=<0:100>  Minimum retro-score for low-confidence validatEVEs
                                even if with high eve-score [default: 10]
  -m --maybe-score-frac=<0:1>   Relative weight of maybe-viral hints in eve-score
                                computation [default: 0.2]
' -> doc
opt <- docopt(doc)
# devel
#opt <- docopt(doc, "issues/17/xxy/results/QWLK01-retro.bed foo.tsv bar.tsv")
#opt <- docopt(doc, "~/Code/projects/detectEVE/issues/17/anopheles-lequime-2017_rvdb-notax/results/APHL01-retro.bed foo.tsv bar.tsv")
opt <- map_at(opt, ~str_detect(., "score"), as.numeric)

# analysis --------------------------------------------------------------------#
cols <- str_split("locus qstart qend sseqid bitscore strand pident length mismatch gapopen sstart send evalue qlen slen qcovhsp desc database taxid lineage", " ")[[1]]
r0 <- read_tsv(opt[["in.bed"]], col_names = cols)
s0 <- read_tsv(opt[["in.fai"]], col_names = c("contig_id", "contig_length"), col_types="cn---")

# dear king philip came over for good soup though
tidy_lineage <- function(lng, rank_order=c("d", "k", "K", "p", "c", "o", "f", "g", "s", "t")){
	# str_match_all on NA returns matrix with one row of NA,NA
	# no match to "" gives empty matrix, which is what we need
	lng <- lng |> replace_na("") 
	ss <- str_match_all(lng, "(\\b\\w)(?::|__)([^,;]+)")
	d <- tibble(.rows= length(ss))
	ii <- rep(seq_along(ss), map_int(ss, nrow))
	mm <- list_c(ss)
	kk <- factor(mm[,2])
	vv <- mm[,3]
	for (k in levels(kk)){
		d[[k]] <- NA
		d[[k]][ii[kk==k]] <- vv[kk==k]
	}
	relocate(d, any_of(rank_order))
}

hints <- c(
	false="ubiquitin|collagen|heat.?shock|NEDD8",
	retro="retro|reverse transcriptase|transposon|transposase",
	viral="virus|viral|virion|endogenous|envelope|coat|capsid",
	maybe="hypothetical|uncharacterized|glycoprotein|polyprotein"
)

tidy_hints <- function(x, hints, ignore=NULL){
	if(!is.null(ignore)) x[ignore] <- ""
	
	map2_dfc(hints, names(hints), function(pattern, type){
		w <- x |> tolower() |> str_extract(pattern)
 		tibble(
			"{type}_hints" := w,
		)
	})
}

clean_desc <- function(x){
	x |> str_remove("^\\S+\\s*") |> str_remove("^acc.*\\|") |> str_remove("(?i)taxid=\\S+")
}

stringify_table <- function(x){
	x <- sort(x, decreasing=T)
	str_c(names(x), " (", x,")", collapse=", ")
}

set_na <- function(x, i){x[i] <- NA; x}

max_count <- function(x){
    max <- names(which.max(table(x, useNA = "no")))
    if(is.null(max)){return(NA_character_)}
    max
}

r1 <- r0 |> 
	mutate(
		tidy_lineage(lineage),
		tidy_hints(desc, hints, ignore=database=="VDB" | k == "Viruses"),
		suggests=coalesce(
			str_c("viral/", set_na(database, database != "VDB")),
			str_c("retro/UDB ", set_na(K, str_detect(K, "pararnavir", negate=T))),
			str_c("viral/UDB ", set_na(k, k != "Viruses")),
			str_c("false-viral/", false_hints),
			str_c("retro/", retro_hints),
			str_c("viral/", viral_hints),
			str_c("maybe-viral/", maybe_hints, " protein"),
			str_c("non-viral/annotated protein of ", k)
		),
	) |> 
	separate(suggests, c("suggests", "because"), sep = "/", extra = "merge") |> 
	relocate(locus, suggests, bitscore, because, desc, database, k)


r2 <- r1 |> 
	group_by(locus) |> 
	filter(bitscore >= opt$min_bitscore_frac * max(bitscore)) |> 
	summarize(
		viral_score = sum(bitscore[suggests == "viral"])/sum(bitscore) * 100,
		maybe_score = sum(bitscore[suggests == "maybe-viral"])/sum(bitscore) * 100,
		retro_score = sum(bitscore[suggests == "retro"])/sum(bitscore) * 100,
		false_score = sum(bitscore[suggests == "false-viral"])/sum(bitscore) * 100,
		eve_score = (viral_score + maybe_score * opt$maybe_score_frac - false_score) |> round(digits = 0),
		top_evalue = evalue[1],
                top_pident = pident[1],
                # top_coverage = qcovhsp[1],
                top_desc = desc[1] |> clean_desc(),
		# top_viral_desc = desc[k == "Viruses"][1] |> clean_desc(),
		top_viral_desc = coalesce(
			desc[k == "Viruses"][1] |> clean_desc(),
			desc[database =="VDB"][1] |> clean_desc()),
                top_viral_lineage = lineage[k == "Viruses"][1],
                top_viral_coverage = round(100*length[k == "Viruses"][1]/slen[k == "Viruses"][1], 2),
		max_count_phylum = max_count(p),
		suggests = stringify_table(table(suggests)),
                because = stringify_table(table(because)),
                eve_length = qlen[1],
		# across(ends_with("hints")	, ~stringify_table(table(.))),
	) |> 
        mutate(
            contig_id = str_remove(locus, "_\\d+-\\d+:[=-]$"),
            confidence = case_when(
		eve_score > opt$eve_score_low & retro_score > opt$retro_score_low ~ "low",
		eve_score > opt$eve_score_high ~ "high",
		eve_score > opt$eve_score_low ~ "low",
		.default = NA)) |> 
	filter(!is.na(confidence)) |> 
        arrange(confidence, -eve_score) |>
        left_join(s0) |> # contig_id/length
        mutate(
            percent_contig = round(eve_length/contig_length * 100, digits=2),
            eve_id = sprintf("%s_EVE%03d", opt$prefix, row_number())
        ) |> 
	relocate(eve_id, confidence, eve_score, suggests, because, locus, eve_length, percent_contig) |> 
	select(-viral_score, -maybe_score, -retro_score, -false_score, -contig_length, -contig_id)

write_tsv(r2, opt[["out.tsv"]])

r3 <- r1 |> 
	inner_join(select(r2, eve_id, eve_score, locus, max_count_phylum, confidence)) |> 
	arrange(desc(eve_id)) |> 
	mutate(eve_id = factor(eve_id, levels=unique(eve_id)))


################################################################################

filtered_eves <- r3[!duplicated(r3$eve_id), ]
confidence_colors <- setNames(ifelse(filtered_eves$confidence == "high", "black", "gray45"), filtered_eves$eve_id)

p1 <- ggplot(r3) +
  geom_point(aes(bitscore, eve_id), size=.3, color="black", show.legend = TRUE) +
  geom_point(aes(bitscore, eve_id, shape=ifelse(k == "Viruses", suggests, NA),
                 color=ifelse(k == "Viruses" & !str_detect(k, "unclassified"), p, NA)), size=3, alpha=.7) +
  scale_x_sqrt() + scale_color_brewer("Viral phylum", palette="Dark2", na.translate = FALSE) +
  scale_shape_manual(values = c(16,NA), labels = c("viral", "all")) +
  guides(color = guide_legend(order = 1), shape = guide_legend(order = 2, title = "Taxonomy")) +
  theme(legend.key = element_blank(), axis.text.y = element_text(color = confidence_colors))
height <- n_distinct(r1$locus)/10+2
if(height >50) height <- 48
ggsave(opt[["out.pdf"]], p1, width=8, height=)

