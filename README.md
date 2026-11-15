# rnaseq-nf

A small bulk RNA-seq pipeline written in **Nextflow DSL2**, laid out the way
nf-core pipelines are (modules, a workflow, configs, profiles, tests and CI).
It is meant to be read, run and extended while learning how production
bioinformatics pipelines are built. **[docs/GUIDE.md](docs/GUIDE.md)** explains every
step, flag, input and output.

```
samplesheet.csv
      │
      ├─► FASTQC ───────────────────────────────┐
      │                                         │
      └─► FASTP (trim) ─► SALMON_QUANT ─┬───────┼─► MULTIQC ─► multiqc_report.html
                              ▲         │       │
          transcripts.fa ─► SALMON_INDEX│       │
          (tx2gene / MAKE_TX2GENE) ─────┘       │
                                        ├─► MERGE_QUANT ─► gene/transcript count + TPM matrices
          samplesheet metadata ─────────┴─► DESEQ2 (tximport) ─► DE table, PCA, MA, volcano
```

## What this practises

| Skill job postings ask for | Where it lives here |
|---|---|
| Workflow manager (Nextflow) | `main.nf`, `workflows/rnaseq.nf`, `modules/local/*.nf` |
| Containers (Docker) + conda | `container` / `conda` lines in every module, `-profile docker` / `-profile conda` |
| Reproducibility | Pinned tool versions, `software_versions.yml`, run reports in `results/pipeline_info/` |
| Testing | `nf-test` pipeline test that checks the *biology* (simulated fold changes), Python unit tests |
| CI / Git | `.github/workflows/ci.yml` runs every test on every push |
| Standard NGS tools | FastQC, fastp, Salmon, MultiQC, tximport + DESeq2 |

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

Any extra columns (`cell_line`, `batch`, …) can be used in the DESeq2 design:
`--design '~ cell_line + condition' --contrast condition,dex,untreated`
(log2FC = dex / untreated). DESeq2 is skipped with a warning when there is only
one condition or no replicates; design mistakes stop the run before it starts.

Sample names must be unique and file-name safe (letters, digits, `.`, `_`, `-`);
the pipeline stops before running anything if they aren't. See
[docs/GUIDE.md](docs/GUIDE.md#4-sample-names-and-conditions).

#### Example: the airway dataset

`scripts/fetch_airway.sh` downloads a small piece of the public "airway" study
(GEO GSE52778: airway smooth-muscle cells from 4 donors, untreated vs dexamethasone,
2×63 bp). It streams the first 1M read pairs of each of the 8 samples from ENA
(about 1 GB instead of about 22 GB), fetches the GENCODE v50 transcriptome, and
writes a samplesheet to `~/bioinfo-data/airway`:

```bash
scripts/fetch_airway.sh
nextflow run main.nf -profile docker \
    --input ~/bioinfo-data/airway/samplesheet.csv \
    --transcriptome ~/bioinfo-data/airway/ref/gencode.v50.transcripts.fa.gz --gencode \
    --design '~ cell_line + condition' --contrast condition,dex,untreated \
    --outdir ~/bioinfo-data/airway/results -work-dir ~/bioinfo-data/airway/work -resume
```

The first reads of a run are not a random sample, which is fine for checking the
pipeline. For numbers you want to report, download the full files and subsample
with `seqtk sample`.

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
├── deseq2/            deseq2_results.tsv, normalized + VST counts, PCA / MA / volcano plots, dds.rds
├── multiqc/multiqc_report.html
└── pipeline_info/     timeline, resource report, trace, DAG, software_versions.yml
```

