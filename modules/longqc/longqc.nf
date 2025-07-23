process LONGQC {
    tag "$meta.id"
    publishDir "${params.outdir}/nanoplot", mode: 'copy'

    // Docker container for LongQC
    container 'quay.io/biocontainers/longqc:1.2.0c--hdfd78af_0'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.png"), emit: plots

    script:
    """
    python longqc sampleqc \\
    -x ont-rapid \\
    -o ${meta.id} \\
    -p ${task.cpus} \\
    ${reads}
    
    """
}