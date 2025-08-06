process QUAST {
    tag "$meta.id"
    publishDir "${params.outdir}/${meta.id}/assembly/quast", mode: 'copy'

    // Docker container for QUAST
    container 'quay.io/biocontainers/quast:5.2.0--py39pl5321h4e691d4_3'

    input:
    tuple val(meta), path(assemblies)
    path reference, stageAs: 'reference.fasta'
    path genes, stageAs: 'genes.gff'
    tuple val(meta2), path(reads, stageAs: 'null')

    output:
    tuple val(meta), path("quast_output/"), emit: results
    tuple val(meta), path("quast_output/report.html"), emit: html
    tuple val(meta), path("quast_output/report.tsv"), emit: tsv
    tuple val(meta), path("quast_output/report.txt"), emit: txt
    tuple val(meta), path("quast_output/report.pdf"), emit: pdf, optional: true
    tuple val(meta), path("quast_output/icarus.html"), emit: icarus, optional: true
    tuple val(meta), path("quast_output/contigs_reports/"), emit: contigs_reports, optional: true
    tuple val(meta), path("quast_output/reads_stats/"), emit: reads_stats, optional: true

    script:
    def prefix = "${meta.id}"
    
    // Handle optional reference genome
    def reference_arg = reference.name != 'NO_FILE' ? "-r ${reference}" : ""
    
    // Handle optional gene annotations
    def genes_arg = genes.name != 'NO_FILE' ? "-g ${genes}" : ""
    
    // Handle optional reads for mapping-based metrics
    def reads_arg = ""
    if (reads && reads.size() >= 2) {
        reads_arg = "-1 ${reads[0]} -2 ${reads[1]}"
    } else if (reads && reads.size() == 1) {
        reads_arg = "--single ${reads[0]}"
    }
    
    // Additional QUAST parameters
    def min_contig = params.quast_min_contig ? "--min-contig ${params.quast_min_contig}" : ""
    def labels = params.quast_labels ? "--labels ${params.quast_labels}" : ""
    def eukaryote = params.quast_eukaryote ? "--eukaryote" : ""
    def large = params.quast_large ? "--large" : ""
    def no_plots = params.quast_no_plots ? "--no-plots" : ""
    def no_html = params.quast_no_html ? "--no-html" : ""
    def no_icarus = params.quast_no_icarus ? "--no-icarus" : ""
    def no_sv = params.quast_no_sv ? "--no-sv" : ""
    def fast = params.quast_fast ? "--fast" : ""

    """
    quast.py \\
        --threads ${task.cpus} \\
        --output-dir quast_output \\
        ${reference_arg} \\
        ${genes_arg} \\
        ${reads_arg} \\
        ${min_contig} \\
        ${labels} \\
        ${eukaryote} \\
        ${large} \\
        ${no_plots} \\
        ${no_html} \\
        ${no_icarus} \\
        ${no_sv} \\
        ${fast} \\
        ${assemblies}
    """
}
