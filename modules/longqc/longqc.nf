process LONGQC {
    tag "$meta.id"
    publishDir "${params.outdir}/longqc", mode: 'copy'

    // Conda package for LongQC
    conda 'bioconda::longqc'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.png"), emit: plots

    script:
    """
    python \$(which longQC.py) sampleqc \\
    -x ont-rapid \\
    -o ${meta.id} \\
    -p ${task.cpus} \\
    ${reads}
    
    """
}