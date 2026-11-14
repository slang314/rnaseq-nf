process SALMON_INDEX {
    label 'process_high'

    conda 'bioconda::salmon=1.10'
    container 'quay.io/biocontainers/salmon:1.10.3--h6dccd9a_2'

    input:
    path fasta

    output:
    path "salmon_index" , emit: index
    path "versions.yml" , emit: versions

    script:
    // --gencode strips the |-separated extra fields from GENCODE headers
    def gencode = params.gencode ? '--gencode' : ''
    """
    salmon index \\
        -t $fasta \\
        -i salmon_index \\
        -k ${params.salmon_k} \\
        -p $task.cpus \\
        $gencode

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        salmon: \$(salmon --version | sed 's/salmon //')
    END_VERSIONS
    """
}
