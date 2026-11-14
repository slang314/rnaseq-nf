process DESEQ2 {
    label 'process_low'

    conda 'bioconda::bioconductor-deseq2=1.50 bioconda::bioconductor-tximport=1.38'
    // Seqera Containers image with both packages (biocontainers images hold one tool each).
    // Rebuild: https://seqera.io/containers/ with bioconductor-deseq2=1.50.2, bioconductor-tximport=1.38.2
    container 'community.wave.seqera.io/library/bioconductor-deseq2_bioconductor-tximport:9f6b62bb2b7e1cb9'

    input:
    path quant_dirs, stageAs: 'quant/*'
    path tx2gene
    path samplesheet
    path gene_names     // optional: [] when there is no gene-name table

    output:
    path "deseq2_results.tsv"   , emit: results
    path "normalized_counts.tsv", emit: normalized
    path "vst_counts.tsv"       , emit: vst
    path "pca.tsv"              , emit: pca
    path "*.pdf"                , emit: plots
    path "deseq2_summary.txt"   , emit: summary
    path "dds.rds"              , emit: rds
    path "R_sessionInfo.log"    , emit: session
    path "versions.yml"         , emit: versions

    script:
    def contrast = params.contrast ? "--contrast '${params.contrast}'" : ''
    def names    = gene_names ? "--gene_names ${gene_names}" : ''
    """
    # R and its packages may write caches under ~; under docker -u there is no home
    export HOME=\$PWD

    deseq2.R \\
        --quant_dir quant \\
        --tx2gene $tx2gene \\
        --samplesheet $samplesheet \\
        --design '${params.design}' \\
        --alpha ${params.deseq2_alpha} \\
        $contrast \\
        $names \\
        --outdir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e 'cat(paste(R.version\$major, R.version\$minor, sep = "."))')
        bioconductor-deseq2: \$(Rscript -e 'cat(as.character(packageVersion("DESeq2")))')
        bioconductor-tximport: \$(Rscript -e 'cat(as.character(packageVersion("tximport")))')
    END_VERSIONS
    """
}
