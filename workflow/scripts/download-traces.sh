#!/usr/bin/env bash
# Usage: download-traces.sh <acc> > out.fna
acc=$1
htm=${2:-$acc.html}
url="https://www.ncbi.nlm.nih.gov/Traces/wgs/${acc}?display=download"

echo "Accession: $acc" >&2
curl -so $htm $url

# check for bad status accessions, e.g. status-replaced, status-suppressed
bad=$(grep -oP '"status-[^"]+"' $htm)
if [ ! -z "$bad" ]; then
    echo Bad accession: $bad >&2
    echo For more info, check $url >&2
    exit 1;
fi

fgz=$(grep -oP 'https://sra-download.ncbi.nlm.nih.gov/[^"]*.fsa_nt.gz' $htm)
if [ -z "$fgz" ]; then
    echo Bad accession: No download links found >&2
    echo For more info, check $url >&2
    exit 1;
fi

echo -e "Downloading:\n$fgz" >&2
curl $fgz | gzip -cd

