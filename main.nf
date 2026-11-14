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
            params.tx2gene ? file(params.tx2gene, checkIfExists: true) : []   // [] = build it
        )

    onComplete:
        log.info(workflow.success
            ? "\nDone. Results in: ${params.outdir}\n  MultiQC report: ${params.outdir}/multiqc/multiqc_report.html\n"
            : "\nPipeline failed. Check .nextflow.log, or rerun with -resume after fixing.\n")
}
