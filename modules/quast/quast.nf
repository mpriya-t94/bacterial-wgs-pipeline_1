// Parameters
params.assemblies = './assemblies/*.fasta'
params.outdir = './quast_results'
params.reference = null
params.genes = null
params.reads1 = null
params.reads2 = null
params.threads = 4

process QUAST {
    container 'quay.io/biocontainers/quast:5.2.0--py39pl5321h4e691d4_3'
    publishDir params.outdir, mode: 'copy'
    cpus params.threads
    
    input:
    path assemblies
    path reference
    path genes
    path reads1
    path reads2
    
    output:
    path "quast_output"
    path "quast_output/report.html"
    path "quast_output/report.tsv"
    path "quast_output/report.txt"
    path "quast_output/report.pdf", optional: true
    path "quast_output/icarus.html", optional: true
    path "quast_output/contigs_reports/", optional: true
    path "quast_output/reads_stats/", optional: true
    
    script:
    def reference_arg = reference.name != 'NO_FILE' ? "-r ${reference}" : ""
    def genes_arg = genes.name != 'NO_FILE' ? "-g ${genes}" : ""
    def reads_arg = ""
    if (reads1.name != 'NO_FILE' && reads2.name != 'NO_FILE') {
        reads_arg = "-1 ${reads1} -2 ${reads2}"
    } else if (reads1.name != 'NO_FILE') {
        reads_arg = "--single ${reads1}"
    }
    
    """
    quast.py \\
        --threads ${task.cpus} \\
        --output-dir quast_output \\
        ${reference_arg} \\
        ${genes_arg} \\
        ${reads_arg} \\
        ${assemblies}
    """
}