process NANOPLOT {
    tag "$meta.id"
    publishDir "${params.outdir}/${meta.id}/qc/${meta.stage}/nanoplot", mode: 'copy'
    
    container "quay.io/biocontainers/nanoplot:1.44.1--pyhdfd78af_0"
    
    input:
    tuple val(meta), path(reads)
    
    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.png") , emit: plots
    tuple val(meta), path("*.txt") , emit: stats
    
    script:
    """
    NanoPlot \\
        --fastq ${reads} \\
        --threads ${task.cpus} \\
        --prefix ${meta.id} \\
        --outdir .
    """
}
