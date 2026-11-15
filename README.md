# rnaseq-nf

[![CI](https://github.com/slang314/rnaseq-nf/actions/workflows/ci.yml/badge.svg)](https://github.com/slang314/rnaseq-nf/actions/workflows/ci.yml)

A bulk RNA-seq pipeline in **Nextflow DSL2** that takes paired-end FASTQ files to
gene-level counts, TPM matrices, and a DESeq2 differential expression analysis.
It runs every tool in a pinned Docker container or conda environment.

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

**Tools:** FastQC · fastp · Salmon · tximport + DESeq2 · MultiQC

**[docs/GUIDE.md](docs/GUIDE.md)** is the full technical guide: every step, flag, input and output.

## Features

- Follows nf-core conventions (module layout, resource labels, `ext.args`, test profiles), but it is an
  independent pipeline and not part of nf-core.
- Quantifies with Salmon against the transcriptome, then summarises to genes with a transcript→gene map.
  With `--gencode`, the map is built from GENCODE FASTA headers.
- DESeq2 with any design that uses samplesheet columns (for example `~ cell_line + condition`).
- Checks the samplesheet and design before any task starts: missing files, duplicate or unsafe
  sample names, and design columns that don't exist all stop the run straight away.
- Reproducible: exact container tags, a `software_versions.yml`, and run reports
  (timeline, trace, resource use) for every run.
- Tests check the results, not just that the run finished: the bundled synthetic data has known
  fold changes, and `nf-test` asserts that Salmon and DESeq2 recover them. CI runs on every push.

## Requirements

- Linux, macOS, or Windows via WSL2
- Java 17 or later and [Nextflow](https://www.nextflow.io/docs/latest/install.html) ≥ 24.04
  (tested with 26.10)
- [Docker](https://docs.docker.com/get-docker/) **or** conda/mamba
  ([Miniforge](https://github.com/conda-forge/miniforge))
- Optional: [nf-test](https://www.nf-test.com/) to run the pipeline tests

On WSL2, keep the repo and data in the Linux filesystem (`~/...`), not under `/mnt/c/`,
because I/O across the Windows boundary is much slower.

## Quick start

```bash
git clone https://github.com/slang314/rnaseq-nf.git
cd rnaseq-nf

# Run on the bundled synthetic data (a few minutes)
nextflow run main.nf -profile test,docker      # or -profile test,conda
```

Results are written to `results/`. Open `results/multiqc/multiqc_report.html` first.

## Usage

### 1. Samplesheet

```csv
sample,fastq_1,fastq_2,condition
untreated_1,/data/untreated_1_R1.fastq.gz,/data/untreated_1_R2.fastq.gz,control
dex_1,/data/dex_1_R1.fastq.gz,/data/dex_1_R2.fastq.gz,treated
```

- Paths can be absolute or relative to the samplesheet.
- Sample names must be unique and use only letters, digits, `.`, `_` and `-`.
- You can add extra columns (`cell_line`, `batch`, …) and use them in `--design`.

See [GUIDE.md §3–4](docs/GUIDE.md#3-inputs) for the full rules.

### 2. Run

With a GENCODE transcriptome (`gencode.vXX.transcripts.fa.gz` from
[gencodegenes.org](https://www.gencodegenes.org/)), `--gencode` builds the
transcript→gene map for you:

```bash
nextflow run main.nf -profile docker \
    --input samplesheet.csv \
    --transcriptome gencode.vXX.transcripts.fa.gz --gencode \
    --design '~ cell_line + condition' --contrast condition,treated,control \
    --outdir results -resume
```

For other transcriptomes, pass a two-column `transcript<TAB>gene` file with `--tx2gene`.
`-resume` reuses every step that has already finished, so after a failure or a parameter
change, only the affected steps rerun.

### Parameters

| Parameter | Default | Description |
|---|---|---|
| `--input` | required | Samplesheet CSV (`sample,fastq_1,fastq_2,condition[,…]`) |
| `--transcriptome` | required | Transcript FASTA (`.fa` or `.fa.gz`) |
| `--tx2gene` | none | `transcript<TAB>gene` map; optional with `--gencode` |
| `--gencode` | `false` | GENCODE headers: use `salmon index --gencode` and build tx2gene from the headers |
| `--salmon_k` | `31` | Salmon index k-mer size (use a smaller k for reads shorter than ~75 bp) |
| `--design` | `~ condition` | DESeq2 design formula; terms must be samplesheet columns |
| `--contrast` | none | `factor,numerator,denominator`; log2FC = numerator / denominator |
| `--deseq2_alpha` | `0.05` | Adjusted p-value cutoff for the summary and plots |
| `--skip_deseq2` | `false` | Stop after quantification |
| `--outdir` | `results` | Output folder |

DESeq2 is skipped with a warning when there is only one condition or no replicates.
To pass extra flags to a tool, set `ext.args` in `conf/modules.config` (see
[GUIDE.md §7](docs/GUIDE.md#7-parameters-profiles-and-resources)).

### Profiles

| Profile | What it does |
|---|---|
| `docker` | Runs each tool in its biocontainers / Seqera container |
| `conda` | Builds one conda environment per tool on first use |
| `test` | Uses the bundled synthetic data (combine it with `docker` or `conda`) |

Resources per process are set in `conf/base.config`. `resourceLimits` caps every task at
6 CPUs / 20 GB by default; raise it to fit your machine.

## Outputs

```
results/
├── fastqc/            per-sample raw-read QC
├── fastp/             trimming reports (json + html)
├── salmon/
│   ├── <sample>/      quant.sf, quant.genes.sf, aux_info/
│   ├── merged_gene_counts.tsv
│   ├── merged_gene_tpm.tsv
│   ├── merged_transcript_{counts,tpm}.tsv
│   └── mapping_summary.tsv
├── deseq2/            deseq2_results.tsv, normalized + VST counts, PCA / MA / volcano plots, dds.rds
├── multiqc/multiqc_report.html
└── pipeline_info/     timeline, resource report, trace, DAG, software_versions.yml
```

[GUIDE.md §6](docs/GUIDE.md#6-outputs) explains each file.

## Example: the airway dataset

`scripts/fetch_airway.sh` downloads part of the public "airway" study (GEO
[GSE52778](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52778): airway smooth-muscle cells
from 4 donors, untreated vs dexamethasone). It streams the first 1M read pairs of each of the
8 samples from ENA (about 1 GB), fetches the GENCODE v50 transcriptome, and writes a samplesheet:

```bash
scripts/fetch_airway.sh          # writes to ~/bioinfo-data/airway
nextflow run main.nf -profile docker \
    --input ~/bioinfo-data/airway/samplesheet.csv \
    --transcriptome ~/bioinfo-data/airway/ref/gencode.v50.transcripts.fa.gz --gencode \
    --design '~ cell_line + condition' --contrast condition,dex,untreated \
    --outdir ~/bioinfo-data/airway/results -work-dir ~/bioinfo-data/airway/work -resume
```

On this subset, 95–96% of reads map. The known glucocorticoid-response genes (DUSP1, KLF15,
TSC22D3, FKBP5, ZBTB16, PER1, CRISPLD2) all come out up-regulated and significant, and PC1
separates treated from untreated. The first reads of a run are not a random sample, so for
real analyses download the full files (and subsample with `seqtk sample` if you want a smaller set).

## Testing

```bash
nf-test test                                  # pipeline tests (uses -profile test,docker)
python3 -m unittest discover -s tests -v      # unit tests for the helper scripts in bin/
```

The test data in `assets/test_data/` is synthetic (`scripts/simulate_test_data.py`):
8 genes with alternative isoforms, 4 samples × 5,000 read pairs, 0.5% sequencing error, and
3% of fragments short enough to read into the adapter. `GENE_UP` is simulated 4× higher and
`GENE_DOWN` 4× lower in the treated samples. The tests check that:

- Salmon recovers both fold changes
- DESeq2 calls exactly those two genes (log2FC ≈ ±2)
- a samplesheet with only one condition skips DESeq2 cleanly

## Troubleshooting

- **A container fails to pull:** biocontainers tags are sometimes rebuilt. Find the current tag at
  `https://quay.io/repository/biocontainers/<tool>?tab=tags` and update the module, or use `-profile conda`.
- **A task fails:** the error message names the task's work directory. `cd` into it and read
  `.command.sh` (what ran), `.command.err` and `.command.log`.
  [GUIDE.md §9](docs/GUIDE.md#9-running-resuming-and-debugging) has more.
- **Tasks are killed for memory (exit 137):** each task is retried once with double memory. If it
  still fails, raise the label's memory or `resourceLimits` in `conf/base.config`.

## Limitations and roadmap

Currently: one contrast per run, no shrunken fold changes (`dds.rds` is saved, so you can run
`lfcShrink` in R), no genome alignment or BAM output, and no decoy-aware index. Full list:
[GUIDE.md §10](docs/GUIDE.md#10-what-it-does-not-do-yet).

Planned:
- Decoy-aware Salmon index (`--decoys`)
- `lfcShrink`, several contrasts per run, DE results in MultiQC
- Samplesheet and parameter validation with `nf-schema`
- A STAR + featureCounts branch (`--aligner star`)
- A warning when the library type Salmon detects doesn't match the expected one
- An AWS Batch profile

## Citations

If you use this pipeline, please cite the tools it runs:
[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/),
[fastp](https://doi.org/10.1093/bioinformatics/bty560),
[Salmon](https://doi.org/10.1038/nmeth.4197),
[tximport](https://doi.org/10.12688/f1000research.7563.2),
[DESeq2](https://doi.org/10.1186/s13059-014-0550-8),
[MultiQC](https://doi.org/10.1093/bioinformatics/btw354), and
[Nextflow](https://doi.org/10.1038/nbt.3820).

## License

MIT. See [LICENSE](LICENSE).
