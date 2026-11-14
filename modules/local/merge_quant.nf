process MERGE_QUANT {
    label 'process_single'

    conda 'conda-forge::python=3.12'
    container 'quay.io/biocontainers/python:3.12.12'

    input:
    path quant_dirs, stageAs: 'quant/*'

    output:
    path "merged_gene_tpm.tsv"         , emit: gene_tpm
    path "merged_gene_counts.tsv"      , emit: gene_counts
    path "merged_transcript_*.tsv"     , emit: transcript
    path "mapping_summary.tsv"         , emit: mapping

    script:
    """
    merge_salmon.py quant/* --outdir .
    """
}
