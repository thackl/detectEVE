#!/usr/bin/env perl

# NCBI SRA WGS ------------------------------
# visit at
# https://www.ncbi.nlm.nih.gov/Traces/wgs/
# download from
# ftp.ncbi.nlm.nih.gov/sra/wgs_aux/

my $rsync = "rsync --copy-links --times --verbose rsync://";
my $sra_wgs_url_pre = "ftp.ncbi.nlm.nih.gov/sra/wgs_aux/";

for my $acc (@ARGV){
    download_sra_wgs($acc)
}

sub download_sra_wgs{
    my ($acc) = @_;
    my $url = url_sra_wgs($acc);
    print $url, "\n";
        
    system("$rsync$url");
    system("cat $acc*.fsa_nt.gz | gzip -cd > $acc.fna");
    system("rm $acc*.fsa_nt.gz")
}

sub url_sra_wgs{
    my ($acc) = @_;
    sprintf("%s/%s/%s/%s/%s.*fsa_nt.gz .", $sra_wgs_url_pre,
            substr($acc, 0, 2), substr($acc, 2,2), $acc, $acc);
}
