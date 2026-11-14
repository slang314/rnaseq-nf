#!/usr/bin/env bash
# Fetch a small subset of the "airway" dataset (Himes et al. 2014, GEO GSE52778,
# SRA SRP033351): 4 airway smooth muscle cell lines, untreated vs dexamethasone,
# Illumina 2x63 bp paired-end. Also fetches the GENCODE human transcriptome.
#
# Subsampling: we stream each FASTQ from ENA and keep the FIRST N read pairs
# (curl | gunzip | head). This avoids downloading ~22 GB, but it is not a random
# sample (reads come from the start of the run). Fine for pipeline testing; use
# seqtk sample on full files for anything quantitative you want to publish.
#
# Usage: scripts/fetch_airway.sh [outdir] [n_pairs]
#   default outdir ~/bioinfo-data/airway, n_pairs 1000000
set -uo pipefail

OUT=${1:-$HOME/bioinfo-data/airway}
N=${2:-1000000}
GENCODE=50
ENA=https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103

mkdir -p "$OUT/fastq" "$OUT/ref"

# run accession : sample name : condition
SAMPLES="
SRR1039508:N61311_untrt:untreated
SRR1039509:N61311_dex:dex
SRR1039512:N052611_untrt:untreated
SRR1039513:N052611_dex:dex
SRR1039516:N080611_untrt:untreated
SRR1039517:N080611_dex:dex
SRR1039520:N061011_untrt:untreated
SRR1039521:N061011_dex:dex
"

fetch_head() {   # url -> first N reads (4N lines) gzipped to $2
    local url=$1 dest=$2
    if [[ -s $dest ]]; then echo "  exists: $dest"; return; fi
    # curl error 23 = "failure writing output" = head closed the pipe: expected, so hidden
    curl -sSL --retry 3 "$url" 2> >(grep -v '(23)' >&2) | gunzip -c 2>/dev/null | head -n $((N * 4)) | pigz -p 2 > "$dest.part"
    local lines
    lines=$(pigz -dc "$dest.part" | wc -l)
    if (( lines != N * 4 )); then echo "  ERROR: $dest has $lines lines, expected $((N * 4))"; return 1; fi
    mv "$dest.part" "$dest"
    echo "  ok: $dest"
}

FAILED=0
echo "sample,fastq_1,fastq_2,condition" > "$OUT/samplesheet.csv"
for entry in $SAMPLES; do
    IFS=: read -r run name cond <<< "$entry"
    sub=${run: -1}   # ENA folder: 00 + last digit of the run accession
    for r in 1 2; do
        fetch_head "$ENA/00$sub/$run/${run}_$r.fastq.gz" "$OUT/fastq/${name}_R$r.fastq.gz" &
    done
    for job in $(jobs -p); do wait "$job" || FAILED=1; done
    echo "$name,fastq/${name}_R1.fastq.gz,fastq/${name}_R2.fastq.gz,$cond" >> "$OUT/samplesheet.csv"
done

TX="$OUT/ref/gencode.v${GENCODE}.transcripts.fa.gz"
if [[ ! -s $TX ]]; then
    curl -sSL --retry 3 -o "$TX" \
        "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${GENCODE}/gencode.v${GENCODE}.transcripts.fa.gz"
fi
echo "transcriptome: $TX ($(du -h "$TX" | cut -f1))"
echo "samplesheet:   $OUT/samplesheet.csv"
if (( FAILED )); then echo "Some downloads failed (truncated transfer?). Rerun: finished files are skipped."; exit 1; fi
