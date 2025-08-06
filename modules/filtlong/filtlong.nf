process FILTLONG {
    tag "$meta.id"
    publishDir "${params.outdir}/${meta.id}/trimming/filtlong", mode: 'copy'

    container 'quay.io/biocontainers/filtlong:0.2.1--h9a82719_0'

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(illumina_reads, stageAs: 'null')
    path assembly, stageAs: 'null'

    output:
    tuple val(meta), path("*_filtered.fastq.gz"), emit: filtered_reads
    tuple val(meta), path("*.log"), emit: log

    script:
    def prefix = "${meta.id}"
    
    // Reference parameters - handle optional inputs
    def illumina_ref = ""
    if (illumina_reads && illumina_reads.size() >= 2) {
        illumina_ref = "-1 ${illumina_reads[0]} -2 ${illumina_reads[1]}"
    }
    
    def assembly_ref = assembly ? "-a ${assembly}" : ""
    
    // Filtering parameters
    def min_length = params.min_length ? "--min_length ${params.min_length}" : ""
    def max_length = params.max_length ? "--max_length ${params.max_length}" : ""
    def keep_percent = params.keep_percent ? "--keep_percent ${params.keep_percent}" : ""
    def target_bases = params.target_bases ? "--target_bases ${params.target_bases}" : ""
    def min_mean_q = params.min_mean_q ? "--min_mean_q ${params.min_mean_q}" : ""
    def min_window_q = params.min_window_q ? "--min_window_q ${params.min_window_q}" : ""
    
    // Scoring weight parameters
    def length_weight = (params.length_weight && params.length_weight != 1) ? "--length_weight ${params.length_weight}" : ""
    def mean_q_weight = (params.mean_q_weight && params.mean_q_weight != 1) ? "--mean_q_weight ${params.mean_q_weight}" : ""
    def window_q_weight = (params.window_q_weight && params.window_q_weight != 1) ? "--window_q_weight ${params.window_q_weight}" : ""
    
    // Read manipulation parameters
    def trim = params.trim ? "--trim" : ""
    def split = params.split ? "--split ${params.split}" : ""
    def window_size = (params.window_size && params.window_size != 250) ? "--window_size ${params.window_size}" : ""
    def verbose = params.verbose ? "--verbose" : ""

    """
    filtlong \\
        ${illumina_ref} \\
        ${assembly_ref} \\
        ${min_length} \\
        ${max_length} \\
        ${keep_percent} \\
        ${target_bases} \\
        ${min_mean_q} \\
        ${min_window_q} \\
        ${length_weight} \\
        ${mean_q_weight} \\
        ${window_q_weight} \\
        ${trim} \\
        ${split} \\
        ${window_size} \\
        ${verbose} \\
        ${reads} \\
        2> ${prefix}_filtlong.log \\
        | gzip > ${prefix}_filtered.fastq.gz
    """
}
