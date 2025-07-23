process FASTQC {
    tag "$meta.id"
    publishDir "${params.outdir}/fastqc", mode: 'copy'

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
