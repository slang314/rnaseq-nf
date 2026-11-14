# rnaseq-nf

A small bulk RNA-seq pipeline written in **Nextflow DSL2**, laid out the way
nf-core pipelines are (modules, a workflow, configs, profiles, tests and CI).
It is meant to be read, run and extended while learning how production
bioinformatics pipelines are built.

```
samplesheet.csv
      │
      ├─► FASTQC ───────────────────────────────┐
      │                                         │
      └─► FASTP (trim) ─► SALMON_QUANT ─┬───────┼─► MULTIQC ─► multiqc_report.html
                              ▲         │       │
          transcripts.fa ─► SALMON_INDEX│       │
          (tx2gene / MAKE_TX2GENE) ─────┘       │
                                        └─► MERGE_QUANT ─► gene/transcript count + TPM matrices
```

## What this practises

| Skill job postings ask for | Where it lives here |
|---|---|
| Workflow manager (Nextflow) | `main.nf`, `workflows/rnaseq.nf`, `modules/local/*.nf` |
| Containers (Docker) + conda | `container` / `conda` lines in every module, `-profile docker` / `-profile conda` |
| Reproducibility | Pinned tool versions, `software_versions.yml`, run reports in `results/pipeline_info/` |
| Testing | `nf-test` pipeline test that checks the *biology* (simulated fold changes), Python unit tests |
| CI / Git | `.github/workflows/ci.yml` runs every test on every push |
| Standard NGS tools | FastQC, fastp, Salmon, MultiQC |

## Setup on Windows (WSL2)

Everything runs inside WSL2 Ubuntu, so you avoid the Windows path and conda
problems.

```bash
# 1. In PowerShell (admin), once:
wsl --install -d Ubuntu-24.04
# 2. Install Docker Desktop -> Settings -> Resources -> WSL integration -> enable Ubuntu
# 3. Inside Ubuntu:
sudo apt update && sudo apt install -y git openjdk-17-jre-headless   # Nextflow needs Java 17+
curl -s https://get.nextflow.io | bash && sudo mv nextflow /usr/local/bin/
curl -fsSL https://get.nf-test.com | bash && sudo mv nf-test /usr/local/bin/
docker run --rm hello-world           # confirms Docker works from WSL
```

Keep the repo and data in the Linux filesystem (`~/projects/...`), not under
`/mnt/c/...`. Disk I/O across the Windows boundary is many times slower.

## Run

```bash
git clone <your-repo-url> rnaseq-nf && cd rnaseq-nf

# Synthetic test, takes a few minutes
nextflow run main.nf -profile test,docker

# The same test through nf-test, with the biology assertions
nf-test test

# Python unit tests only (no Nextflow needed)
python3 -m unittest discover -s tests -v
```

No Docker? Use `-profile test,conda` instead. Nextflow builds one conda
environment per tool the first time.

### Real data

Samplesheet (paths can be absolute, or relative to the samplesheet):

```csv
sample,fastq_1,fastq_2,condition
untreated_1,/data/untreated_1_R1.fastq.gz,/data/untreated_1_R2.fastq.gz,control
dex_1,/data/dex_1_R1.fastq.gz,/data/dex_1_R2.fastq.gz,treated
```

With a GENCODE human transcriptome (download `gencode.vXX.transcripts.fa.gz`
from gencodegenes.org), the transcript→gene map is built for you:

```bash
nextflow run main.nf -profile docker \
    --input samplesheet.csv \
    --transcriptome gencode.vXX.transcripts.fa.gz --gencode \
    --outdir results_real -resume
```

A good public dataset to start with is the "airway" study (GEO GSE52778: airway
smooth-muscle cells, with and without dexamethasone). Get the reads from SRA with
`fasterq-dump`, and use `seqtk sample` to take about 2M read pairs per sample for a fast first run.

`-resume` reuses every step that already finished, so you can fix one
step and rerun without starting over.

## Outputs

```
results/
├── fastqc/            per-sample raw-read QC
├── fastp/             trimming reports (json + html)
├── salmon/
│   ├── <sample>/      quant.sf, quant.genes.sf, aux_info/
│   ├── merged_gene_counts.tsv       -> input for DESeq2 / edgeR
│   ├── merged_gene_tpm.tsv
│   ├── merged_transcript_{counts,tpm}.tsv
│   └── mapping_summary.tsv
├── multiqc/multiqc_report.html
└── pipeline_info/     timeline, resource report, trace, DAG, software_versions.yml
```

## About the test data

`assets/test_data/` is fully synthetic (`scripts/simulate_test_data.py`):
8 genes, some with alternative isoforms that share exons, 4 samples × 5,000
read pairs, 0.5% sequencing errors, and 3% short fragments that run into the adapter.
`GENE_UP` is simulated 4× higher and `GENE_DOWN` 4× lower in "treated". The
nf-test checks that Salmon recovers those changes, so the test catches
results that are wrong as well as runs that crash.

## Container tags

The Docker image tags are pinned to biocontainers releases. If one fails to
pull (tags get rebuilt now and then), look up the current tag at
`https://quay.io/repository/biocontainers/<tool>?tab=tags`, or switch to `-profile conda`.

## Extensions (roughly in order of difficulty)

1. **Decoy-aware Salmon index**: add the genome as decoys (`--decoys`); this
   is what nf-core/rnaseq does, and it reduces false mappings.
2. **DESeq2 module**: an R process that reads `merged_gene_counts.tsv` plus the
   samplesheet `condition` column and writes a results table and a volcano plot.
3. **Samplesheet validation with the `nf-schema` plugin** and a `nextflow_schema.json`.
4. **STAR + featureCounts branch** behind a `--aligner star` parameter (the AstraZeneca
   posting lists both).
5. **Stranded-library check**: read `lib_format_counts.json` and warn when the
   detected library type doesn't match the expected one.
6. **Cloud profile**: an `awsbatch` profile, plus running from Seqera Platform.
7. **Port to single-cell**: swap Salmon for `simpleaf`/alevin-fry and add a Scanpy QC step.
