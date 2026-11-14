process MAKE_TX2GENE {
    label 'process_single'

    conda 'conda-forge::python=3.12'
    container 'quay.io/biocontainers/python:3.12.12'

    input:
    path fasta

    output:
    path "tx2gene.tsv", emit: tsv

    script:
    """
    make_tx2gene.py $fasta > tx2gene.tsv
    """
}
