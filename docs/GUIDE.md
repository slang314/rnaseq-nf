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

**Goal:** turn raw bulk RNA-seq FASTQ files into a gene × sample count matrix and a TPM matrix,
test which genes differ between conditions (DESeq2), and produce a QC report that tells you
whether you should trust any of it.

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
                                                                    ├──► MERGE_QUANT ──► matrices
   samplesheet.csv (metadata columns) ──────────────────────────────┴──► DESEQ2 ──► DE table, PCA
                                                    (skipped with one condition or no replicates)
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

### 4.2 How conditions reach DESeq2

The samplesheet file itself is an input to the DESEQ2 step. Every column except `fastq_1`/`fastq_2`
is **sample metadata**, and the design formula can use any of them. You can add columns freely:

```csv
sample,fastq_1,fastq_2,condition,cell_line
N61311_untrt,…,…,untreated,N61311
N61311_dex,…,…,dex,N61311
```

Three parameters say what to test:

| Parameter | Default | Meaning |
|---|---|---|
| `--design` | `~ condition` | R model formula over samplesheet columns. **Put the variable you test last**, and the variables you adjust for (batch, donor, cell line) before it |
| `--contrast` | none | `factor,numerator,denominator`. log2FC = log2(numerator / denominator), so positive = higher in the numerator |
| `--skip_deseq2` | false | Don't run DESeq2 at all |

**Always give `--contrast` for real data.** Without it, a 2-level factor is compared with the
**alphabetically first level as the baseline** (DESeq2's default). For the test data that is `control`, which
is what you want. For airway it would be `dex`, and every fold change would come out with the wrong sign.
The DESeq2 log states which direction it used, and `deseq2_summary.txt` repeats it.

The other conditions in `meta.condition` are still not used by any per-sample step. FastQC, fastp and Salmon
treat each sample alone, and quantifying a sample should not depend on its group.

### 4.3 More than one factor

Airway is a **paired design**: each of 4 cell lines was split into untreated and dex. Donors differ
a lot from each other, much more than the dex effect for many genes. `~ cell_line + condition`
estimates the dex effect *within* each cell line, so donor differences are removed instead of being
counted as noise. Run it with:

```bash
--design '~ cell_line + condition' --contrast condition,dex,untreated
```

The same pattern covers batch (`~ batch + condition`), patients (`~ patient + treatment`), and so on.
Interaction designs (`~ genotype * treatment`) parse, but the `--contrast` form only tests a simple
two-level comparison of one factor. Interaction terms need R directly (load `dds.rds`).

### 4.4 What is checked before the run

`planDeseq2()` in `main.nf` reads the samplesheet before any task starts, so a mistake in the design
fails in seconds, not after an hour of quantification:

| Situation | What happens |
|---|---|
| `--design` names a column that isn't in the samplesheet | **error**, lists the available columns |
| A design column is empty for some samples | **error**, names them |
| `--contrast` is not `factor,numerator,denominator`, or names a level that doesn't exist | **error** |
| The factor has 3+ levels and no `--contrast` | **error**, asks which two to compare |
| The factor has **only one level** | **warning**, DESeq2 skipped, the rest runs |
| A compared level has **only one sample** | **warning**, DESeq2 skipped, the rest runs |

`bin/deseq2.R` repeats the essential checks, and also stops if the design has no residual degrees of freedom
(more model coefficients than the samples can support, e.g. `~ cell_line + condition` with one sample per cell line).

### 4.5 Only one condition

The pipeline runs normally: FastQC, trimming, quantification, matrices and MultiQC do not care
how many groups there are. DESeq2 is skipped with a warning. Leave `condition` empty or give every row the
same value. What you get is **quantification and QC, not differential expression**. That is a perfectly
good use: expression atlases, checking a library prep, comparing samples by eye with TPM, feeding a
clustering/PCA, or quantifying samples that will be compared with others later.

What you **cannot** do with one condition is test for differential expression: there is nothing to
contrast. Two related limits:

- **Two conditions, one sample each:** DESeq2 can't estimate biological variance, so the pipeline
  skips it. You can still rank genes by fold change from the TPM matrix, as a description, not a test.
- **Comparing against samples from another run:** possible, but then batch and condition are
  confounded: you cannot tell biology from batch effects, and no design formula can separate them.

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
| Command | `make_tx2gene.py transcripts.fa.gz --names gene_names.tsv > tx2gene.tsv` |
| Output | `tx2gene.tsv`: `ENST00000832824.1<TAB>ENSG00000290825.2` (one row per transcript, no header) |
| | `gene_names.tsv`: `gene_id<TAB>gene_name` (one row per gene, with header); used by DESEQ2 |

A GENCODE header has 8 `|`-separated fields:
`>ENST00000832824.1|ENSG00000290825.2|-|-|DDX11L16-260|DDX11L16|1379|lncRNA|`
(transcript ID, gene ID, Havana gene, Havana transcript, transcript name, **gene name**, length,
biotype). The script keeps fields 1 and 2 for tx2gene, and fields 2 and 6 for the name table. It exits with an
error if a header has no `|`, i.e. if the FASTA isn't GENCODE-style. Gene IDs keep their version suffix
(`.2`). If you give your own `--tx2gene`, there is no name table, and the DE results have IDs only.

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

### 5.7 DESEQ2: differential expression (`modules/local/deseq2.nf`, `bin/deseq2.R`)

| | |
|---|---|
| Purpose | Test, gene by gene, whether expression differs between two levels of a factor, adjusting for the other design variables |
| Input | all Salmon folders (`ch_quant_dirs`, the same collected list MERGE_QUANT gets), `tx2gene.tsv`, the samplesheet, `gene_names.tsv` or `[]` |
| Command | `deseq2.R --quant_dir quant --tx2gene tx2gene.tsv --samplesheet samplesheet.csv --design '~ cell_line + condition' --alpha 0.05 --contrast condition,dex,untreated --gene_names gene_names.tsv --outdir .` |
| Output | see the table below; published to `results/deseq2/` |
| Container | `community.wave.seqera.io/library/bioconductor-deseq2_bioconductor-tximport:9f6b62bb2b7e1cb9` (R 4.5.3, DESeq2 1.50.2, tximport 1.38.2) |
| Resources | `process_low` (2 CPUs, 4 GB). DESeq2 runs single-threaded here |

**Why this container is different.** Biocontainers images hold one tool each, and this step needs two
R packages. [Seqera Containers](https://seqera.io/containers/) builds an image from any list of conda
packages, and its tag is a hash of that list, so the same list always gives the same image. nf-core uses
it for every module with several tools. It was built with one API call (`POST https://wave.seqera.io/v1alpha2/container`,
packages `bioconductor-deseq2=1.50.2`, `bioconductor-tximport=1.38.2`, `freeze: true`).

`[]` for an optional input: when there's no gene-name table, the workflow passes an empty list as
`gene_names`. Nextflow stages nothing, `gene_names` is false in the script block, and the flag is left out.
This is the standard nf-core way to make a `path` input optional.

What `deseq2.R` does, in order:

1. **Sample table.** Reads the samplesheet, drops the FASTQ columns, turns the design variables into
   factors, and **relevels** the tested factor so the contrast's denominator is the reference level.
2. **tximport** (`type = "salmon"`, `dropInfReps = TRUE`). Reads each `quant.sf`, sums transcripts to genes
   with tx2gene, and records each gene's **average transcript length per sample**, weighted by isoform
   abundance. `dropInfReps` skips Salmon bootstrap replicates, which this pipeline doesn't create.
3. **`DESeqDataSetFromTximport`.** Rounds the counts and turns those lengths into per-gene, per-sample
   normalization offsets. A gene whose isoform mix shifts to longer transcripts gets more reads without
   being more expressed, and the offsets account for that.
4. **Low-count filter.** Keeps genes with ≥ 10 reads in at least as many samples as the smallest compared
   group. Genes with almost no reads can't be significant and only add to the multiple-testing burden.
5. **`estimateSizeFactors`.** The **median-of-ratios** normalization: each sample's scaling factor is the median,
   over genes, of its count divided by that gene's geometric mean across samples. It assumes most genes don't
   change, which makes it robust to a few very large changes. TPM is not, and that is why
   the synthetic unchanged genes look 0.82× in TPM but ≈ 1× in DESeq2.
6. **`estimateDispersions`.** For each gene, how much counts vary between replicates beyond Poisson noise
   (the negative binomial "dispersion"). DESeq2 fits a trend of dispersion against mean expression and
   **shrinks** each gene's estimate towards it, which matters when there are only 2-4 replicates. If the
   trend fit fails (it needs many genes; the synthetic data has 8), the script falls back to gene-wise estimates.
7. **`nbinomWaldTest` + `results(contrast = …, alpha = 0.05)`.** Fits the model for every gene and tests
   the log2 fold change ≠ 0. `padj` is Benjamini-Hochberg adjusted. `alpha` also tunes DESeq2's
   **independent filtering**: very low-mean genes get `padj = NA` when removing them gives more discoveries.
8. **VST** (`vst(blind = TRUE)`). Variance-stabilizing transform for PCA and clustering: on log2(counts) low
   counts are very noisy; VST evens that out. `blind = TRUE` ignores the design, the recommended setting for QC.

Outputs:

| File | Content |
|---|---|
| `deseq2_results.tsv` | One row per tested gene, **sorted by padj**: `gene_id`, `gene_name` (with GENCODE), `baseMean` (mean normalized count), `log2FoldChange`, `lfcSE`, `stat` (Wald statistic), `pvalue`, `padj` |
| `deseq2_summary.txt` | Design, contrast and direction, number of samples and genes, number up/down at `padj < alpha` |
| `normalized_counts.tsv` | Counts divided by the normalization factors: comparable across samples |
| `vst_counts.tsv` | VST values (log2-like scale) for PCA, heatmaps, clustering |
| `pca.pdf`, `pca.tsv` | PC1 vs PC2 of the 500 most variable genes, colored by the tested factor |
| `ma_plot.pdf` | log2FC against mean expression; significant genes highlighted |
| `volcano.pdf` | log2FC against −log10 p-value; the top 15 significant genes labeled |
| `dds.rds` | The fitted DESeq2 object: `dds <- readRDS("dds.rds")` in R to try other contrasts or plots |
| `R_sessionInfo.log` | Exact R and package versions |

Reading the results: `log2FoldChange` +2 means 4× higher in the numerator. `padj < 0.05` means the gene is
called significant at a 5% false discovery rate. A small padj with a tiny fold change can be real but
biologically uninteresting, and a large fold change with `padj = NA` was filtered as too lowly expressed.
The fold changes are **not shrunk**. For ranking genes by effect size, `lfcShrink()` (apeglm or ashr) on `dds.rds`
gives more stable estimates for low-count genes.

### 5.8 MULTIQC: one report (`modules/local/multiqc.nf`, `assets/multiqc_config.yml`)

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
├── deseq2/          deseq2_results.tsv, deseq2_summary.txt, normalized_counts.tsv, vst_counts.tsv,
│                    pca.{pdf,tsv}, ma_plot.pdf, volcano.pdf, dds.rds    (when DESeq2 runs)
├── multiqc/         multiqc_report.html, multiqc_report_data/
└── pipeline_info/   timeline_*.html  report_*.html  trace_*.txt  dag_*.html  software_versions.yml
```

`pipeline_info/` is about the **run**, not the biology: `timeline` shows when each task ran and for how
long (find the bottleneck), `report` shows CPU/memory per task (tune the labels), `trace` is the same as
a table, `dag` draws the workflow graph, and `software_versions.yml` lists every tool version (for
methods sections).

**Using the merged counts elsewhere** (edgeR, limma-voom): they are fractional, and they do not carry
Salmon's per-sample effective lengths. The DESEQ2 step avoids this by importing with tximport. For other tools,
tximport with `countsFromAbundance = "lengthScaledTPM"` is the recommended way to get counts.

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
| `--design` | `~ condition` | DESeq2 model formula over samplesheet columns (section 4.2) |
| `--contrast` | none | `factor,numerator,denominator` |
| `--deseq2_alpha` | 0.05 | padj cutoff for the summary, plots and independent filtering |
| `--skip_deseq2` | false | quantification and QC only |
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
| `process_low` | 2 | 4 GB | FASTQC, DESEQ2 |
| `process_medium` | 4 | 8 GB | FASTP, SALMON_QUANT |
| `process_high` | 8 | 16 GB | SALMON_INDEX |

- `resourceLimits = [cpus: 6, memory: 20.GB]` caps every request (defaults fit a 6-core, 24 GB workstation; raise them for bigger machines).
- Memory is `× task.attempt`: if a task is killed for using too much memory (exit 137 etc.), it is retried once with double memory.
- Nextflow runs as many tasks at once as fit under the limits. With 4-CPU Salmon tasks and a 6-CPU cap,
  only **one quant runs at a time**. Setting `process_medium` to 3 CPUs would let two run side by side.

---

## 8. Tests

| Test | Command | What it proves |
|---|---|---|
| Python unit tests | `python3 -m unittest discover -s tests -v` | `make_tx2gene.py` (including the gene-name table) and `merge_salmon.py` parse and write correctly (seconds, no Nextflow) |
| Pipeline test | `nextflow run main.nf -profile test,docker` | The whole pipeline runs on the synthetic data |
| nf-test | `nf-test test` | Three tests (`tests/main.nf.test`): **(1)** TPM: GENE_UP > 2.5× up, GENE_DOWN > 2.5× down, GENE_A unchanged. **(2)** DESeq2: GENE_UP log2FC in 1.5..2.5 and GENE_DOWN in −2.5..−1.5 (truth ±2), both padj < 0.01, the six unchanged genes \|log2FC\| < 0.5. **(3)** one-condition samplesheet: run succeeds, no DESEQ2 task, quantification done |
| CI | push to GitHub | `.github/workflows/ci.yml` runs the unit tests then nf-test with Docker |

On the synthetic data DESeq2 gives GENE_UP +2.03, GENE_DOWN −1.87, unchanged genes −0.05 to +0.27, and
exactly those 2 genes at padj < 0.05.

The synthetic data has a known answer (`scripts/simulate_test_data.py`, `assets/test_data/truth.tsv`).
That is the point of it: a test that only checks "it finished" passes on a pipeline that gives wrong numbers.

Real-data check (airway, `scripts/fetch_airway.sh`, first 1M pairs × 8 samples,
`--design '~ cell_line + condition' --contrast condition,dex,untreated`), run by hand, not an automated test:

- 9,477 genes tested, **512 at padj < 0.05** (335 up, 177 down). Full-depth airway gives far more; 1M
  reads per sample means much less power.
- The known glucocorticoid-response genes are all significant and up: ZBTB16 +6.4, KLF15 +5.8, FKBP5 +3.7,
  TSC22D3 +3.2, DUSP1 +2.9, PER1 +2.6, CRISPLD2 +2.6 (log2FC). GAPDH −0.22, padj 0.56.
- Top hits by padj: GPX3, ADAMTS1, CXCL12 (down), COL1A1 (down), NNMT, FKBP5.
- PCA: PC1 separates dex from untreated in every cell line; PC2 separates out cell line N080611, which is
  the donor effect the `cell_line` term removes.
- ACTB is also called up (+0.84, padj 2e-6). Unexplained: possibly real, possibly from using the first
  reads of each run instead of a random sample.

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
    --design '~ cell_line + condition' --contrast condition,dex,untreated \
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

- **One contrast per run**, no interaction tests, no shrunken fold changes (use `dds.rds` in R).
- **CI tests the Docker profile only.** The conda profile (including DESeq2) has been checked by hand.
- **Basic samplesheet validation only**: required columns, files exist, unique and file-name-safe
  sample names. No schema for the other columns or parameters yet (nf-schema, README extension 3).
- **No lane merging**: one row per sample.
- **No decoy-aware index**, so some intronic/intergenic reads can be assigned to transcripts.
- **No genome alignment / BAM**, so no coverage tracks, junction QC, or rRNA/intron-rate QC (RSeQC, Qualimap).
- **No strandedness check**: Salmon detects the library type, but nothing warns you if it isn't what you expected.
- **Not tuned for full-size data**: 20-40M read pairs per sample will work, but quant time and
  `process_medium` memory should be rechecked on a real run.
