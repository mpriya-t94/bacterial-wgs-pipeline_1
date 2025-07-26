process SKESA {
    tag "$meta.id"
    publishDir "${params.outdir}/skesa", mode: 'copy'

    container 'quay.io/biocontainers/skesa:2.4.0--he1b5a44_0'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.fasta"), emit:contigs

    script:
    """

    skesa\\
        --reads ${reads.join(' --reads')} \\
        --use_paired_ends \\
        --cores ${params.skesa_cores ?: 12} \\
        --memory ${params.skesa_memory ?: 32} \\
        --min_contig ${params.skesa_min_contig ?: 500} \\
        --contigs_out "${meta.id}.fasta"

    """

}