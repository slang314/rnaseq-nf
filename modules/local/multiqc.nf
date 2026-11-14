process MULTIQC {
    label 'process_single'

    conda 'bioconda::multiqc=1.25'
    container 'quay.io/biocontainers/multiqc:1.25.1--pyhdfd78af_0'

    input:
    path multiqc_files, stageAs: '?/*'
    path config

    output:
    path "multiqc_report.html", emit: report
    path "*_data"            , emit: data

    script:
    """
    multiqc --force --config $config --filename multiqc_report.html .
    """
}
