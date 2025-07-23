process FASTQC {
    tag "$meta.id"
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    // Docker container for FASTQC
    container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.zip"), emit: zip

    script:
    """
    fastqc --outdir . --threads ${task.cpus} ${reads}
    """
}
