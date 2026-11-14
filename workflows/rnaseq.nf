include { FASTQC        } from '../modules/local/fastqc'
include { FASTP         } from '../modules/local/fastp'
include { MAKE_TX2GENE  } from '../modules/local/make_tx2gene'
include { SALMON_INDEX  } from '../modules/local/salmon_index'
include { SALMON_QUANT  } from '../modules/local/salmon_quant'
include { MERGE_QUANT   } from '../modules/local/merge_quant'
include { DESEQ2        } from '../modules/local/deseq2'
include { MULTIQC       } from '../modules/local/multiqc'

workflow RNASEQ {
    take:
    ch_reads        // channel: [ meta, [ fastq_1, fastq_2? ] ]
    transcriptome   // path: transcript FASTA
    tx2gene_file    // path, or [] to build from GENCODE headers
    samplesheet     // path: the samplesheet (its metadata columns feed the DESeq2 design)
    run_deseq2      // boolean: decided in main.nf from the samplesheet and --skip_deseq2

    main:
    ch_versions = channel.empty()

    // 1. QC on raw reads
    FASTQC(ch_reads)
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    // 2. Adapter / quality trimming
    FASTP(ch_reads)
    ch_versions = ch_versions.mix(FASTP.out.versions.first())

    // 3. transcript -> gene map (given, or built from GENCODE headers)
    //    Gene names come only from GENCODE headers; [] means "no names" to DESEQ2
    if (tx2gene_file) {
        ch_tx2gene    = channel.value(tx2gene_file)
        ch_gene_names = channel.value([])
    } else {
        MAKE_TX2GENE(transcriptome)
        ch_tx2gene    = MAKE_TX2GENE.out.tsv
        ch_gene_names = MAKE_TX2GENE.out.names
    }

    // 4. Salmon: build the index once, quantify every sample against it
    SALMON_INDEX(transcriptome)
    SALMON_QUANT(FASTP.out.reads, SALMON_INDEX.out.index, ch_tx2gene)
    ch_versions = ch_versions.mix(SALMON_INDEX.out.versions, SALMON_QUANT.out.versions.first())

    // 5. One count matrix + one TPM matrix for the whole experiment
    //    .collect() waits for every sample and emits one list of folders
    ch_quant_dirs = SALMON_QUANT.out.results.map { _meta, dir -> dir }.collect()
    MERGE_QUANT(ch_quant_dirs)

    // 6. Differential expression (skipped with one condition, no replicates, or --skip_deseq2)
    if (run_deseq2) {
        DESEQ2(ch_quant_dirs, ch_tx2gene, samplesheet, ch_gene_names)
        ch_versions = ch_versions.mix(DESEQ2.out.versions)
    }

    // 7. One report for everything
    ch_multiqc = channel.empty()
        .mix(FASTQC.out.zip.map { _meta, f -> f })
        .mix(FASTP.out.json.map { _meta, f -> f })
        .mix(SALMON_QUANT.out.results.map { _meta, d -> d })
        .collect()
    MULTIQC(ch_multiqc, file(params.multiqc_config))

    // Software versions, one file, for the record
    ch_versions
        .unique()
        .collectFile(name: 'software_versions.yml', storeDir: "${params.outdir}/pipeline_info", newLine: false)

    emit:
    gene_tpm = MERGE_QUANT.out.gene_tpm
    report   = MULTIQC.out.report
}
