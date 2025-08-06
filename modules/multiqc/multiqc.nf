process MULTIQC {
    tag "multiqc_${stage}"
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    // Docker container for MultiQC
    container 'quay.io/biocontainers/multiqc:1.30--pyhdfd78af_0'

    input:
    path (files)
    val (stage)

    output:
    path ("multiqc_report.html"), emit: html
    path ("multiqc_report_data/"), emit: data

    script:
    """
    multiqc . \\
    --title "Bacterial WGS ${stage} QC Report" \\
    --filename multiqc_report.html \\
    --dirs \\
    --dirs-depth 2
    """
}