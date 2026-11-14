process FASTP {
    tag "$meta.id"
    label 'process_medium'

    conda 'bioconda::fastp=0.23'
    container 'quay.io/biocontainers/fastp:0.23.4--h5f740d0_0'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.trim.fastq.gz"), emit: reads
    tuple val(meta), path("*.fastp.json")   , emit: json
    tuple val(meta), path("*.fastp.html")   , emit: html
    path "versions.yml"                     , emit: versions

    script:
    def p    = meta.id
    def args = task.ext.args ?: ''
    def io   = meta.single_end
        ? "-i ${reads[0]} -o ${p}.trim.fastq.gz"
        : "-i ${reads[0]} -I ${reads[1]} -o ${p}_1.trim.fastq.gz -O ${p}_2.trim.fastq.gz --detect_adapter_for_pe"
    """
    fastp \\
        $io \\
        --thread $task.cpus \\
        --json ${p}.fastp.json \\
        --html ${p}.fastp.html \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastp: \$(fastp --version 2>&1 | sed -e 's/fastp //')
    END_VERSIONS
    """
}
