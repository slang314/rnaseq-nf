#!/usr/bin/env nextflow
/*
 * rnaseq-nf: a small, nf-core-style bulk RNA-seq pipeline
 *   FastQC -> fastp -> Salmon (index + quant) -> merged matrices -> MultiQC
 *
 * Quick test (needs Docker or conda):
 *   nextflow run main.nf -profile test,docker
 */

include { RNASEQ } from './workflows/rnaseq'

// Resolve a samplesheet path: absolute paths and URLs are used as-is,
// relative paths are taken relative to the samplesheet's own folder.
def resolveReads(String p, sheetDir) {
    if (p.startsWith('/') || p.contains('://')) return file(p, checkIfExists: true)
    return file("${sheetDir}/${p}", checkIfExists: true)
}

// Decide whether DESeq2 can run, from the samplesheet alone (before any task starts).
// Errors for things that are clearly mistakes (a design column that doesn't exist);
// a warning + skip for designs that are valid but can't be tested (one condition,
// no replicates), so QC and quantification still finish.
def planDeseq2(List rows, header) {
    if (params.skip_deseq2) return false

    def vars = params.design.toString().replaceFirst(/^\s*~/, '').split(/[+*:]/)*.trim().findAll { v -> v && !(v in ['0', '1']) }
    if (!vars) error "--design '${params.design}' has no variables"
    def missing = vars.findAll { v -> !(v in header) }
    if (missing) error "--design uses column(s) not in the samplesheet: ${missing.join(', ')} (columns: ${header.join(', ')})"
    vars.each { v ->
        def empty = rows.findAll { row -> !row[v] }*.sample
        if (empty) error "Column '${v}' is in the design but empty for: ${empty.join(', ')}"
    }

    def factor = vars[-1]
    def compare
    if (params.contrast) {
        def c = params.contrast.toString().tokenize(',')*.trim()
        if (c.size() != 3) error "--contrast must be factor,numerator,denominator (got '${params.contrast}')"
        factor  = c[0]
        if (!(factor in vars)) error "--contrast factor '${factor}' is not in --design '${params.design}'"
        compare = c[1..2]
    }
    def counts = rows.countBy { row -> row[factor] }
    if (counts.size() < 2) {
        log.warn "Only one '${factor}' value (${counts.keySet().join('')}) in the samplesheet: skipping DESeq2. QC and quantification still run."
        return false
    }
    if (compare) {
        def unknown = compare.findAll { l -> !(l in counts.keySet()) }
        if (unknown) error "--contrast level(s) ${unknown.join(', ')} not found in column '${factor}' (values: ${counts.keySet().join(', ')})"
    } else if (counts.size() > 2) {
        error "'${factor}' has ${counts.size()} values (${counts.keySet().join(', ')}): say which two to compare with --contrast ${factor},<numerator>,<denominator>"
    } else {
        compare = counts.keySet().toList()
    }
    def unreplicated = compare.findAll { l -> counts[l] < 2 }
    if (unreplicated) {
        log.warn "No replicates for ${factor} = ${unreplicated.join(', ')}: DESeq2 can't estimate variance, skipping it. QC and quantification still run."
        return false
    }
    return true
}

workflow {
    main:
        if (!params.input)         error "Missing --input (samplesheet CSV: sample,fastq_1,fastq_2,condition)"
        if (!params.transcriptome) error "Missing --transcriptome (transcript FASTA)"
        if (!params.tx2gene && !params.gencode) {
            error "Give --tx2gene (transcript<TAB>gene), or use a GENCODE FASTA with --gencode so it can be built"
        }

        def sheet    = file(params.input, checkIfExists: true)
        def sheetDir = sheet.parent

        // Read the whole sheet up front (not inside a channel) so problems that need
        // every row at once, like duplicate names, stop the run before any task starts.
        def rows = sheet.splitCsv(header: true, strip: true)
        if (!rows) error "Samplesheet ${sheet} has no sample rows"
        rows.each { row ->
            if (!row.sample) error "Samplesheet row is missing 'sample': ${row}"
            // sample names become file and folder names: keep them file-name safe
            if (!(row.sample ==~ /[A-Za-z0-9][A-Za-z0-9._-]*/)) {
                error "Sample name '${row.sample}' may only contain letters, digits, '.', '_' and '-' (and must start with a letter or digit)"
            }
        }
        def dups = rows.countBy { row -> row.sample }.findAll { _name, n -> n > 1 }.keySet()
        if (dups) error "Duplicate sample names in ${sheet}: ${dups.join(', ')}. One row per sample (merge lane FASTQs first)"

        def runDeseq2 = planDeseq2(rows, rows[0].keySet())

        ch_samples = channel
            .fromList(rows)
            .map { row ->
                if (!row.fastq_1) error "Sample '${row.sample}' is missing fastq_1"
                def single = !row.fastq_2
                def meta   = [id: row.sample, condition: row.condition ?: 'NA', single_end: single]
                def reads  = single ? [resolveReads(row.fastq_1, sheetDir)]
                                    : [resolveReads(row.fastq_1, sheetDir), resolveReads(row.fastq_2, sheetDir)]
                [meta, reads]
            }

        RNASEQ(
            ch_samples,
            file(params.transcriptome, checkIfExists: true),
            params.tx2gene ? file(params.tx2gene, checkIfExists: true) : [],  // [] = build it
            sheet,
            runDeseq2
        )

    onComplete:
        log.info(workflow.success
            ? "\nDone. Results in: ${params.outdir}\n  MultiQC report: ${params.outdir}/multiqc/multiqc_report.html\n  DESeq2 (if run): ${params.outdir}/deseq2/deseq2_results.tsv\n"
            : "\nPipeline failed. Check .nextflow.log, or rerun with -resume after fixing.\n")
}
