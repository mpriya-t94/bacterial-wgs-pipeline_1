process SKESA {
    tag "$meta.id"
    publishDir "${params.outdir}/${meta.id}/assembly/skesa", mode: 'copy'

    container 'quay.io/biocontainers/skesa:2.5.1--h077b44d_3'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.fasta"), emit:contigs

    script:

    def prefix = "${meta.id}"
    def reads_list = reads instanceof List ? reads.join(',') : reads

    """
    skesa \\
        --reads ${reads_list} \\
        --use_paired_ends \\
        --cores ${task.cpus} \\
        --memory ${task.memory.toGiga()} \\
        --min_contig ${params.skesa_min_contig ?: 200} \\
        --contigs_out ${prefix}.fasta
    """

}