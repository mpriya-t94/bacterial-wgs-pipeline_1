process NANOPLOT {
    tag "$meta.id"
    publishDir "${params.outdir}/nanoplot", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.png"), emit: plots

    script:
    """
    NanoPlot --outdir . --threads ${task.cpus} --fastq ${reads}
    """
}