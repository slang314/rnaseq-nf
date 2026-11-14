process SALMON_QUANT {
    tag "$meta.id"
    label 'process_medium'

    conda 'bioconda::salmon=1.10'
    container 'quay.io/biocontainers/salmon:1.10.3--h6dccd9a_2'

    input:
    tuple val(meta), path(reads)
    path index
    path tx2gene

    output:
    tuple val(meta), path("${meta.id}"), emit: results
    path "versions.yml"                , emit: versions

    script:
    def args   = task.ext.args ?: ''
    def inputs = meta.single_end ? "-r ${reads[0]}" : "-1 ${reads[0]} -2 ${reads[1]}"
    """
    salmon quant \\
        -i $index \\
        -l A \\
        $inputs \\
        -g $tx2gene \\
        -p $task.cpus \\
        -o ${meta.id} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        salmon: \$(salmon --version | sed 's/salmon //')
    END_VERSIONS
    """
}
