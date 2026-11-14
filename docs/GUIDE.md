# rnaseq-nf technical guide

How the pipeline works from end to end: what goes in, what each tool does, what comes out,
and which flags control it. Read it next to the code; every section names the file it describes.

Contents
1. [What the pipeline is for](#1-what-the-pipeline-is-for)
2. [Nextflow in five minutes](#2-nextflow-in-five-minutes)
3. [Inputs](#3-inputs)
4. [Sample names and conditions](#4-sample-names-and-conditions)
5. [Step by step](#5-step-by-step)
6. [Outputs](#6-outputs)
7. [Parameters, profiles and resources](#7-parameters-profiles-and-resources)
8. [Tests](#8-tests)
9. [Running, resuming and debugging](#9-running-resuming-and-debugging)
10. [What it does not do (yet)](#10-what-it-does-not-do-yet)

---

## 1. What the pipeline is for

**Goal:** turn raw bulk RNA-seq FASTQ files into a gene × sample count matrix (for differential
expression) and a TPM matrix (for looking at expression levels), with a QC report that tells you
whether you should trust them.

**Approach:** *pseudo-alignment / selective alignment* with **Salmon** against the transcriptome,
not genome alignment. Salmon maps each read to the transcripts it is compatible with, then uses
an EM algorithm to split multi-mapping reads between isoforms. This is fast (minutes per sample on a
desktop) and is what nf-core/rnaseq's `--pseudo_aligner salmon` mode and most DE tutorials use.
The cost: no BAM file, so no genome browser views, no novel splice junctions, no variant calls.

```
samplesheet.csv ──► one [meta, reads] item per sample
                        │
                        ├──► FASTQC (raw reads) ─────────────────────────────┐
                        │                                                    │
                        └──► FASTP (trim/filter) ──► SALMON_QUANT ──┬────────┼──► MULTIQC
                                                         ▲          │        │
   transcripts.fa(.gz) ──► SALMON_INDEX (once) ──────────┘          │        │
   transcripts.fa(.gz) ──► MAKE_TX2GENE (if --gencode, no tx2gene) ─┘        │
                                                                    └──► MERGE_QUANT ──► matrices
```

Each box is one Nextflow **process** in `modules/local/`. Per-sample steps run in parallel.
The index is built once and shared by every sample.

---

## 2. Nextflow in five minutes

The pipeline uses only a few Nextflow ideas, and the rest of this guide depends on them.

| Concept | What it is | Where |
|---|---|---|
| **process** | One tool run: declared `input:`, `output:`, and a bash `script:`. Each run (a *task*) gets its own folder under `work/` | `modules/local/*.nf` |
| **channel** | A queue of items flowing between processes. Here usually `[meta, files]` tuples | `workflows/rnaseq.nf` |
| **meta map** | A small dictionary that travels with each sample's files: `[id: 'N61311_dex', condition: 'dex', single_end: false]`. Processes read `meta.id` to name outputs | built in `main.nf` |
| **workflow** | Wires processes together by passing one's outputs (`FASTP.out.reads`) as another's inputs | `workflows/rnaseq.nf` |
| **operators** | Channel transformations: `.map{}` reshape items, `.collect()` wait for all items and emit one list, `.mix()` merge channels | `workflows/rnaseq.nf` |
| **value channel** | A channel that can be read any number of times. The index and tx2gene are used by every sample, so each SALMON_QUANT task gets the same index | `SALMON_INDEX.out.index` |
| **label** | A tag on a process (`process_medium`) that `conf/base.config` maps to CPUs/memory | top of each module |
| **ext.args** | Extra tool flags injected from config, so you change flags without editing the module | `conf/modules.config` |
| **publishDir** | Copies selected outputs from `work/` to `results/`. Everything else stays in `work/` | `conf/modules.config` |
| **profile** | A named bundle of config. `-profile test,docker` = test data + run in containers | `nextflow.config` |
| **bin/** | Scripts here are put on `PATH` inside every task, so modules call `merge_salmon.py` directly | `bin/` |

How a task runs: Nextflow creates `work/ab/cdef…/`, symlinks the inputs in (*staging*), writes the
script to `.command.sh`, runs it inside the container (`docker run -u uid:gid -v …`), and checks that the
declared outputs exist. If they don't, the task fails, even if the tool exited 0.

---

## 3. Inputs

### 3.1 Samplesheet (`--input`)

CSV with a header. One row per sample:

```csv
sample,fastq_1,fastq_2,condition
N61311_untrt,fastq/N61311_untrt_R1.fastq.gz,fastq/N61311_untrt_R2.fastq.gz,untreated
N61311_dex,fastq/N61311_dex_R1.fastq.gz,fastq/N61311_dex_R2.fastq.gz,dex
```

| Column | Required | Meaning |
|---|---|---|
| `sample` | yes | Sample name. Becomes `meta.id`. See section 4 for the rules |
| `fastq_1` | yes | Read 1 (or the only read file for single-end) |
| `fastq_2` | no | Read 2. **Leave empty for single-end**; the sample then gets `meta.single_end = true` |
| `condition` | no | Experimental group. Stored as `meta.condition` (default `NA`). See section 4 |

Paths may be absolute, URLs (`s3://`, `https://`, which Nextflow downloads for you), or relative.
**Relative paths are resolved against the samplesheet's folder**, not your current directory
(`resolveReads()` in `main.nf`). That is why the test samplesheet uses bare filenames. Every file is
checked to exist (`checkIfExists: true`) before anything runs.

Parsing happens in `main.nf`. `sheet.splitCsv(header: true, strip: true)` reads the whole file up front
into a list of maps (`row.sample`, `row.fastq_1`…), so checks that need every row at once (duplicate
names) can stop the run before any task starts. `channel.fromList(rows)` then turns the list into a
channel, and `.map {}` turns each row into `[meta, [r1, r2]]`.

### 3.2 Transcriptome (`--transcriptome`)

A FASTA of transcript sequences (`.fa` or `.fa.gz`), e.g. GENCODE `gencode.v50.transcripts.fa.gz`
(all ~250k human transcripts, 176 MB). This is what Salmon maps to. Use the **transcripts** file, not
the genome.

### 3.3 Transcript-to-gene map (`--tx2gene` or `--gencode`)

Salmon quantifies transcripts. To also get gene-level numbers it needs to know which transcripts
belong to which gene: a two-column TSV, `transcript_id<TAB>gene_id`, no header.

- Give it yourself with `--tx2gene file.tsv` (the test data does this), **or**
- Use a GENCODE FASTA with `--gencode`. Then `MAKE_TX2GENE` builds it from the headers.

The transcript IDs in tx2gene must match the names Salmon puts in `quant.sf`, or the gene file
comes out wrong. `--gencode` handles this for GENCODE: see `SALMON_INDEX` and `MAKE_TX2GENE` below.

---

## 4. Sample names and conditions

### 4.1 What the sample name does

`sample` becomes `meta.id`, and `meta.id` is used as:

- the **file prefix** of every per-sample output (`N61311_dex.fastp.json`, `N61311_dex_1.trim.fastq.gz`)
- the **Salmon output folder** name (`results/salmon/N61311_dex/`)
- the **column header** in every merged matrix (`merge_salmon.py` takes the sample name from the folder name)
- the **task tag** in the log (`SALMON_QUANT (N61311_dex)`)
- the sample name in the MultiQC report

Rules that follow from this. `main.nf` checks the first two before any task starts:

1. **Unique.** Two rows with the same name would give two Salmon folders with the same name, so
   the merge would fail at the very end, and publishDir would silently overwrite one sample's
   reports with the other's. Duplicates stop the run: `Duplicate sample names in …: ctrl_rep1`.
2. **File-name safe.** Letters, digits, `.`, `_`, `-`, starting with a letter or digit
   (regex `[A-Za-z0-9][A-Za-z0-9._-]*`). A name like `ctrl rep1` stops the run.
3. **Not ending in `_R1`, `_R2` or `.fastp`.** MultiQC strips those (`extra_fn_clean_exts` in
   `assets/multiqc_config.yml`) to merge R1/R2 rows, so `sampleA_R1` would show up as `sampleA`.
4. **One sample = one row.** If a library was sequenced on several lanes, `cat` the lane FASTQs
   together first (gzip files can be concatenated directly). nf-core/rnaseq merges same-name rows
   automatically; this pipeline does not.

Matrix columns are sorted **alphabetically by sample name**, not in samplesheet order
(`samples = sorted(...)` in `merge_salmon.py`). Anything downstream should match columns by
name, never by position.

Names that encode the design (`N61311_dex`, `N61311_untrt`) make the matrices readable, but nothing in
the pipeline parses them. The nf-test is an exception: it groups samples by name prefix (`ctrl`, `trt`).

### 4.2 What the condition does: currently nothing

`condition` is read into `meta.condition` and carried along with every sample, **but no current
step uses it.** FastQC, fastp and Salmon process each sample independently, and per-sample quantification
does not depend on the experimental group (and shouldn't). The column exists for the planned
DESeq2 step (README extension 2), which will need the sample → group table.

So today, "indicating DE conditions" just means writing them in the samplesheet. Comparing the
groups is done after the pipeline, on `merged_gene_counts.tsv`.

### 4.3 More than one factor

Real designs often have more than one factor. Airway is a **paired design**: each of 4 cell lines
was split into untreated and dex. The right DE model is `~ cell_line + condition`, which compares
treated vs untreated *within* each cell line and removes the large differences between donors.
The samplesheet has no column for that today. The natural fix, when the DESeq2 module is added, is to
allow extra columns (e.g. `batch`, `cell_line`, `patient`), keep them in `meta`, and pass a design
formula as a parameter (`--design '~ cell_line + condition'`).

### 4.4 Only one condition

The pipeline runs normally: FastQC, trimming, quantification, matrices and MultiQC do not care
how many groups there are. Leave `condition` empty (it becomes `NA`) or give every row the same value.
What you get is **quantification and QC, not differential expression**. That is a perfectly good use:
expression atlases, checking a library prep, comparing samples by eye with TPM, feeding a
clustering/PCA, or quantifying samples that will be compared with others later.

What you **cannot** do with one condition is test for differential expression: there is nothing to
contrast. A future DESeq2 step should check for this and skip itself with a clear message (`design
needs at least 2 levels of condition`) instead of crashing. Two related limits to know:

- **Two conditions, one sample each:** DESeq2 can't estimate biological variance, so its p-values
  are not meaningful. You can still rank genes by fold change, as a description, not a test.
- **Comparing against samples from another run:** possible, but then batch and condition are
  confounded; you cannot tell biology from batch effects.

---

## 5. Step by step

The exact command for each step is in the module's `script:` block. `$task.cpus` comes from the
label, `$args` from `ext.args` in `conf/modules.config`.

### 5.1 FASTQC: raw read quality (`modules/local/fastqc.nf`)

| | |
|---|---|
| Purpose | Per-file quality overview of the **raw** reads: base quality by position, GC content, adapter content, duplication, overrepresented sequences |
| Input | `[meta, [R1, R2]]` from the samplesheet (untrimmed) |
| Command | `fastqc --threads N --quiet R1.fastq.gz R2.fastq.gz` |
| Output | `<file>_fastqc.html` (open in a browser) and `<file>_fastqc.zip` (raw data that MultiQC reads), one pair per FASTQ |
| Resources | `process_low` (2 CPUs, 4 GB) |

FastQC runs in parallel with fastp, not after it, so it describes the data as sequenced. fastp's own
report describes the trimmed data. For RNA-seq, some FastQC "fails" are normal and expected: high
duplication (highly expressed genes), biased base composition in the first ~10 bases (random hexamer
priming), and overrepresented sequences (mitochondrial/ribosomal transcripts).

### 5.2 FASTP: adapter trimming and read filtering (`modules/local/fastp.nf`)

| | |
|---|---|
| Purpose | Remove adapter sequence and drop low-quality or too-short reads, keeping R1/R2 in sync |
| Input | `[meta, [R1, R2]]` (raw) |
| Command (paired) | `fastp -i R1 -I R2 -o <id>_1.trim.fastq.gz -O <id>_2.trim.fastq.gz --detect_adapter_for_pe --thread N --json <id>.fastp.json --html <id>.fastp.html --qualified_quality_phred 20 --length_required 30` |
| Command (single) | same with `-i R1 -o <id>.trim.fastq.gz`, no `-I/-O`, no `--detect_adapter_for_pe` |
| Output | trimmed FASTQs (passed to Salmon, **not published**), `.fastp.json` (for MultiQC), `.fastp.html` |
| Resources | `process_medium` (4 CPUs, 8 GB) |

The flags:

| Flag | Effect |
|---|---|
| *(default)* adapter trimming | For paired reads fastp finds adapters by **overlapping R1 and R2**: if the insert is shorter than the read length, both reads run into adapter at the same point |
| `--detect_adapter_for_pe` | Also detect the adapter sequence itself (as in single-end mode). Catches adapters the overlap test misses |
| `--qualified_quality_phred 20` | Raises the threshold for a "qualified" base from Q15 to Q20. This is a **filter, not a trimmer**: a read pair is **dropped** if more than 40% of its bases (`--unqualified_percent_limit`, default 40) are below Q20. fastp does not trim low-quality read ends unless you add `--cut_front`/`--cut_tail`/`--cut_right` |
| `--length_required 30` | Drop reads shorter than 30 bp after trimming (default 15). Very short reads map ambiguously |
| *(default)* N filter | Drop reads with more than 5 `N` bases (`--n_base_limit`) |
| *(default)* poly-G trimming | Enabled automatically for NextSeq/NovaSeq data (2-colour chemistry reads "no signal" as G) |

In a pair, if either read fails, both are dropped, so R1 and R2 stay matched line for line. That
matters because Salmon reads the two files in lockstep. On airway (HiSeq 2000, 2×63 bp), 95-98% of
pairs pass.

### 5.3 MAKE_TX2GENE: build the transcript → gene map (`modules/local/make_tx2gene.nf`, `bin/make_tx2gene.py`)

Runs only when `--tx2gene` is not given (`if (tx2gene_file) … else MAKE_TX2GENE(...)` in the workflow).

| | |
|---|---|
| Input | the transcriptome FASTA |
| Command | `make_tx2gene.py transcripts.fa.gz > tx2gene.tsv` |
| Output | `tx2gene.tsv`: `ENST00000832824.1<TAB>ENSG00000290825.2` |

A GENCODE header has 8 `|`-separated fields:
`>ENST00000832824.1|ENSG00000290825.2|-|-|DDX11L16-260|DDX11L16|1379|lncRNA|`
(transcript ID, gene ID, Havana gene, Havana transcript, transcript name, **gene name**, length,
biotype). The script keeps fields 1 and 2. It exits with an error if a header has no `|`, i.e. if the FASTA isn't
GENCODE-style. Gene IDs keep their version suffix (`.2`). Gene **names** are not in the output, so
to look up `DUSP1` you need to map names to IDs from field 6.

### 5.4 SALMON_INDEX: index the transcriptome (`modules/local/salmon_index.nf`)

| | |
|---|---|
| Purpose | Build the k-mer index that Salmon maps against. Done once per run, shared by all samples |
| Input | the transcriptome FASTA |
| Command | `salmon index -t transcripts.fa.gz -i salmon_index -k 31 -p N [--gencode]` |
| Output | `salmon_index/` (published to `results/reference/`) |
| Resources | `process_high` (8 CPUs, 16 GB, capped to 6 CPUs by `resourceLimits`). Human GENCODE v50: 6 min, 3.1 GB peak RAM |

| Flag | Effect |
|---|---|
| `-k 31` (`--salmon_k`) | k-mer length. A read must share exact 31-mers with a transcript to be considered. Shorter k = more sensitive for short or error-prone reads, but slower and less specific. 31 is the default and suits reads ≥ ~75 bp. For 50 bp reads consider `--salmon_k 23` |
| `--gencode` | Cut each header at the first `|`, so transcripts are named `ENST…` instead of the full 8-field header. This is what makes the names in `quant.sf` match `tx2gene.tsv` |

The index has **no decoys**. Reads from intronic or intergenic regions (pre-mRNA, genomic DNA) can
be forced onto a similar transcript. A decoy-aware index adds the genome sequence so Salmon can
recognize and discard those reads (README extension 1).

### 5.5 SALMON_QUANT: quantify each sample (`modules/local/salmon_quant.nf`)

| | |
|---|---|
| Purpose | Assign reads to transcripts and estimate how many came from each |
| Input | `[meta, trimmed reads]` from FASTP, the index, `tx2gene.tsv` |
| Command | `salmon quant -i salmon_index -l A -1 <id>_1.trim.fastq.gz -2 <id>_2.trim.fastq.gz -g tx2gene.tsv -p N -o <id> --seqBias --gcBias` |
| Output | a folder named `<id>/` (published to `results/salmon/<id>/`) |
| Resources | `process_medium` (4 CPUs, 8 GB). Airway, 1M pairs: ~5 min, 8 GB peak |

How the index and tx2gene reach every sample: `SALMON_QUANT(FASTP.out.reads, SALMON_INDEX.out.index,
ch_tx2gene)`. The first input is a queue channel (one item per sample); the other two are single
values. Nextflow pairs every sample with the same index.

| Flag | Effect |
|---|---|
| `-l A` | Detect the **library type** automatically from the first reads: `I` (inward-facing pair) + `U` (unstranded), `SF`/`SR` (stranded, read 1 sense / antisense). The result is in `lib_format_counts.json`. Airway: `IU`, compatible ratio 1.0 |
| `-1/-2` or `-r` | Paired or single-end reads (chosen from `meta.single_end`) |
| `-g tx2gene.tsv` | Also write `quant.genes.sf`: transcripts summed per gene (counts and TPM summed; length = TPM-weighted mean of transcript lengths) |
| `--seqBias` | Model and correct for sequence-specific bias at fragment ends (random-hexamer priming bias) |
| `--gcBias` | Model and correct for fragment GC-content bias (PCR). Most useful for paired-end data |
| *(default)* selective alignment | Since 1.0, Salmon scores candidate mappings with a real alignment score (`--validateMappings` is on), which rejects poor matches that plain k-mer matching would accept |

Bias correction is the slow part: on airway it took ~3.5 of the ~5 min per sample, because it works over
every transcript in the index (~250k), however many reads there are. Mapping 1M pairs took ~1 min.

What is in the output folder:

| File | Content |
|---|---|
| `quant.sf` | Transcript level: `Name  Length  EffectiveLength  TPM  NumReads` |
| `quant.genes.sf` | The same columns, per gene |
| `aux_info/meta_info.json` | `num_processed`, `num_mapped`, `percent_mapped`, versions, detected library type |
| `lib_format_counts.json` | Library type evidence (how many fragments fit each orientation) |
| `libParams/flenDist.txt` | Fragment-length distribution (paired-end) |
| `logs/salmon_quant.log` | Full log with timestamps; first place to look if numbers seem off |
| `cmd_info.json` | The exact command Salmon recorded |

The two kinds of numbers:
- **NumReads**: estimated number of fragments from that transcript. **Fractional** (multi-mapping
  reads are split, e.g. `30.0020`). This is what DE tools model.
- **TPM**: transcripts per million, corrected for **effective length** (`EffectiveLength`: length
  minus the fragment-length distribution, so short transcripts aren't overcounted). Comparable
  **within** a sample (gene A vs gene B). Across samples it is reasonable for a quick look, but don't run
  statistics on it.

### 5.6 MERGE_QUANT: one matrix per measure (`modules/local/merge_quant.nf`, `bin/merge_salmon.py`)

| | |
|---|---|
| Input | **all** Salmon folders at once: `SALMON_QUANT.out.results.map { _meta, dir -> dir }.collect()` |
| Command | `merge_salmon.py quant/* --outdir .` |
| Output | `merged_gene_counts.tsv`, `merged_gene_tpm.tsv`, `merged_transcript_counts.tsv`, `merged_transcript_tpm.tsv`, `mapping_summary.tsv` |

`.collect()` is the "wait for everyone" step: it waits until every sample has finished, then emits one
list, so MERGE_QUANT runs once. `stageAs: 'quant/*'` puts all the folders in a `quant/` subfolder in
the task directory. The script takes each sample name from its folder name, sorts the samples
alphabetically, and writes a feature × sample table (missing features = 0).

### 5.7 MULTIQC: one report (`modules/local/multiqc.nf`, `assets/multiqc_config.yml`)

| | |
|---|---|
| Input | FastQC zips + fastp JSONs + Salmon folders, collected with `.mix().collect()` |
| Command | `multiqc --force --config multiqc_config.yml --filename multiqc_report.html .` |
| Output | `multiqc_report.html`, `multiqc_report_data/` (the parsed numbers as TSV/JSON) |

MultiQC scans the folder, recognizes each tool's output files, and builds one interactive report.
`stageAs: '?/*'` gives each input its own numbered subfolder, so files from different samples
with the same name (every Salmon folder has `meta_info.json`) don't collide. The config sets the
section order and strips `_R1`/`_R2`/`.fastp` so each sample appears as one row.

What to check in the report: fastp "% passed filter" and adapter content; Salmon "% mapped" and the
fragment-length distribution; FastQC per-base quality. A sample that is off from the rest on any of these
is worth a closer look before DE.

---

## 6. Outputs

```
results/
├── fastqc/          <sample>_R{1,2}_fastqc.{html,zip}     raw-read QC, per FASTQ
├── fastp/           <sample>.fastp.{html,json}            trimming report (trimmed reads are not kept)
├── make_tx2gene/    tx2gene.tsv                           only when built from GENCODE headers
├── reference/       salmon_index/                          reusable index (~1 GB for human)
├── salmon/
│   ├── <sample>/                    quant.sf, quant.genes.sf, aux_info/, logs/ …
│   ├── merged_gene_counts.tsv       gene × sample, estimated counts  → DE input
│   ├── merged_gene_tpm.tsv          gene × sample, TPM               → plots, eyeballing
│   ├── merged_transcript_counts.tsv transcript × sample
│   ├── merged_transcript_tpm.tsv
│   └── mapping_summary.tsv          sample, reads_processed, percent_mapped
├── multiqc/         multiqc_report.html, multiqc_report_data/
└── pipeline_info/   timeline_*.html  report_*.html  trace_*.txt  dag_*.html  software_versions.yml
```

`pipeline_info/` is about the **run**, not the biology: `timeline` shows when each task ran and for how
long (find the bottleneck), `report` shows CPU/memory per task (tune the labels), `trace` is the same as
a table, `dag` draws the workflow graph, and `software_versions.yml` lists every tool version (for
methods sections).

**Taking counts into DESeq2/edgeR:** the merged counts are fractional, and they do not carry Salmon's
per-sample effective lengths. The recommended route is `tximport` (R) on the per-sample `quant.sf` +
`tx2gene.tsv`, which passes length offsets to DESeq2 (`DESeqDataSetFromTximport`). Rounding
`merged_gene_counts.tsv` and using `DESeqDataSetFromMatrix` also works and is common, but it ignores
differences in isoform usage between samples.

---

## 7. Parameters, profiles and resources

### Parameters (`nextflow.config`, set with `--name value`)

| Parameter | Default | Meaning |
|---|---|---|
| `--input` | required | samplesheet CSV |
| `--transcriptome` | required | transcript FASTA (.fa/.fa.gz) |
| `--tx2gene` | none | transcript→gene TSV. Required unless `--gencode` |
| `--gencode` | false | GENCODE headers: `salmon index --gencode` and build tx2gene automatically |
| `--salmon_k` | 31 | k-mer size for the index |
| `--outdir` | `results` | where published files go |
| `--multiqc_config` | `assets/multiqc_config.yml` | MultiQC settings |

Single dash = Nextflow option (`-profile`, `-resume`, `-work-dir`), double dash = pipeline parameter.

### Changing tool flags without touching modules

`conf/modules.config` sets `ext.args` per process. To drop Salmon's bias correction for one run, use a
custom config file:

```groovy
// my.config
process { withName: 'SALMON_QUANT' { ext.args = '' } }
```
`nextflow run main.nf -profile docker -c my.config …`. Configs given with `-c` override the defaults.

### Profiles

| Profile | Effect |
|---|---|
| `docker` | Run each task in its `container` image as your user (`-u $(id -u):$(id -g)`), so outputs aren't owned by root |
| `conda` | Build one conda env per tool from the `conda` line instead (slower first run, no Docker needed) |
| `test` | Load `conf/test.config`: bundled synthetic data, resources capped at 2 CPUs / 6 GB |

Combine with commas: `-profile test,docker`.

### Resources (`conf/base.config`)

| Label | CPUs | Memory | Used by |
|---|---|---|---|
| `process_single` | 1 | 2 GB | MAKE_TX2GENE, MERGE_QUANT, MULTIQC |
| `process_low` | 2 | 4 GB | FASTQC |
| `process_medium` | 4 | 8 GB | FASTP, SALMON_QUANT |
| `process_high` | 8 | 16 GB | SALMON_INDEX |

- `resourceLimits = [cpus: 6, memory: 20.GB]` caps every request to what WSL2 is given.
- Memory is `× task.attempt`: if a task is killed for using too much memory (exit 137 etc.), it is retried once with double memory.
- Nextflow runs as many tasks at once as fit under the limits. With 4-CPU Salmon tasks and a 6-CPU cap,
  only **one quant runs at a time**. Setting `process_medium` to 3 CPUs would let two run side by side.

---

## 8. Tests

| Test | Command | What it proves |
|---|---|---|
| Python unit tests | `python3 -m unittest discover -s tests -v` | `make_tx2gene.py` and `merge_salmon.py` parse and write correctly (seconds, no Nextflow) |
| Pipeline test | `nextflow run main.nf -profile test,docker` | The whole pipeline runs on the synthetic data |
| nf-test | `nf-test test` | Runs the test profile **and checks the biology**: GENE_UP > 2.5× up, GENE_DOWN > 2.5× down, GENE_A unchanged (`tests/main.nf.test`) |
| CI | push to GitHub | `.github/workflows/ci.yml` runs the unit tests then nf-test with Docker |

The synthetic data has a known answer (`scripts/simulate_test_data.py`, `assets/test_data/truth.tsv`).
That is the point of it: a test that only checks "it finished" passes on a pipeline that gives wrong numbers.

Real-data check (airway, `scripts/fetch_airway.sh`, first 1M pairs × 8 samples): the known
glucocorticoid-response genes are up in dex in all 4 cell lines (log2 TPM ratio: KLF15 3.5-4.9,
ZBTB16 3.7-4.5, FKBP5 3.1-4.0, DUSP1 2.3-3.6, TSC22D3 2.5-4.2), GAPDH flat. Not yet an automated test.

---

## 9. Running, resuming and debugging

```bash
# synthetic test
nextflow run main.nf -profile test,docker

# real data (airway); keep data, work/ and results outside the repo
scripts/fetch_airway.sh                          # → ~/bioinfo-data/airway
nextflow run main.nf -profile docker \
    --input ~/bioinfo-data/airway/samplesheet.csv \
    --transcriptome ~/bioinfo-data/airway/ref/gencode.v50.transcripts.fa.gz --gencode \
    --outdir ~/bioinfo-data/airway/results -work-dir ~/bioinfo-data/airway/work -resume
```

**`-resume`**: every task has a hash of its inputs, script and container. On a rerun, tasks whose hash
is unchanged are taken from `work/` (shown as `cached`). Change a Salmon flag and only Salmon and
the steps after it rerun. Delete `work/` and everything reruns.

**A task failed.** The error message names the task's work folder. Go there:

```bash
cd work/ab/cdef1234…
cat .command.sh     # the exact script that ran (flags filled in)
cat .command.err    # the tool's stderr: usually the answer
cat .command.log    # stdout + stderr
cat .exitcode
ls -la              # were the inputs staged? which outputs exist?
bash .command.run   # rerun it by hand, in the same container
```

Common exit codes: `1` tool error (read `.command.err`), `126` script not executable (`chmod +x bin/*`),
`127` command not found (wrong container / not on PATH), `137` killed for memory (retried automatically
once), "Missing output file(s)": the tool ran but didn't write the file name the module declares.

`nextflow lint .` checks the syntax of every `.nf` and config file without running anything.
`nextflow log` lists past runs, and `nextflow log <run> -f name,status,workdir` maps tasks to folders.

---

## 10. What it does not do (yet)

- **No differential expression.** It stops at matrices; `condition` is carried but unused (section 4).
- **Basic samplesheet validation only**: required columns, files exist, unique and file-name-safe
  sample names. No schema for the other columns or parameters yet (nf-schema, README extension 3).
- **No lane merging**: one row per sample.
- **No decoy-aware index**, so some intronic/intergenic reads can be assigned to transcripts.
- **No genome alignment / BAM**, so no coverage tracks, junction QC, or rRNA/intron-rate QC (RSeQC, Qualimap).
- **No strandedness check**: Salmon detects the library type, but nothing warns you if it isn't what you expected.
- **Not tuned for full-size data**: 20-40M read pairs per sample will work, but quant time and
  `process_medium` memory should be rechecked on a real run.
