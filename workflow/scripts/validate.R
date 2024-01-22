# command line interface ------------------------------------------------------#
library(docopt)
'Validate putatEVEs based on retroblast hits

Usage:
  validate.R <in.o6> <out.tsv> <out.pdf>
' -> doc
opt <- docopt(doc)

# devel
# opt[["in.o6"]] <- "results/QMGA01-retro.o6"
# analysis --------------------------------------------------------------------#
library(tidyverse)
cols <- str_split("qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle staxids sscinames sskingdoms skingdoms sphylums", " ")[[1]]
r0 <- read_tsv(opt[["in.o6"]], col_names = cols)

# rename to Sebastians scheme                 
r1 <- r0 |> 
  group_by(qseqid) |> 
  slice_min(evalue, n=10) |> 
  mutate(
    evidence=case_when(
      sskingdoms == "Viruses" ~ "1 exogenous virus hit",
      str_detect(stitle, "virus|viral|endogenous") ~ "2 eve-ish description",
      str_detect(stitle, "hypothetical|uncharacterized") ~ "3 hypothetical protein"
    ),
    confidence=case_when(
      any(str_sub(evidence, 1, 1) == "1") ~ "high",
      any(!is.na(evidence)) ~ "low")
    ) |> 
  # only keep putatEVEs with at least 1 eve-ish evidence hit
  filter(any(!is.na(confidence)))

# get top virus hits
r2 <- r1 |>
  group_by(qseqid) |> 
  filter(!is.na(evidence) & length > 83) |> 
  arrange(-length, evalue) |>
  slice_head(n=1) |> 
  arrange(confidence, evidence, skingdoms, sphylums, evalue)

# plot hit overview
#p1 <- 
  
r3 <- r1 |> 
  mutate(
    qseqid=factor(qseqid, levels=rev(unique(r2$qseqid)))) |>
  filter(!is.na(qseqid))

p1 <- ggplot(r3) +
  geom_point(aes(length, qseqid), size=.3, color="black") +
  geom_point(aes(length, qseqid, shape=evidence, color=ifelse(confidence == "high", sphylums, NA)), data=r2, size=3, alpha=.7) +
  scale_x_sqrt() + scale_color_brewer("Viral phylum", palette="Dark2", na.value="grey50")
ggsave(opt[["out.pdf"]], p1, width=8, height=n_distinct(r1$qseqid)/10+2)


write_tsv(r2, opt[["out.tsv"]])
