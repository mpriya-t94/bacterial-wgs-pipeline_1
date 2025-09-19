process HYBRACTER {
    tag "$meta.id"
    publishDir "${params.outdir}/${meta.id}/assembly/hybracter", mode: 'copy'

    container 'quay.io/gbouras13/hybracter:0.11.0'

    input:
    tuple val(meta), path(long_reads)
    tuple val(meta2), path(short_reads)

output:
    tuple val(meta), path("*/FINAL_OUTPUT/complete/*_final.fasta"), emit: complete_assemblies, optional: true
    tuple val(meta), path("*/FINAL_OUTPUT/incomplete/*_final.fasta"), emit: incomplete_assemblies, optional: true
    tuple val(meta), path("*/FINAL_OUTPUT/complete/*_plasmid.fasta"), emit: plasmids, optional: true
    tuple val(meta), path("*/FINAL_OUTPUT/hybracter_summary.tsv"), emit: summary
    tuple val(meta), path("*/FINAL_OUTPUT/processing_summary.tsv"), emit: processing_summary

    script:
    def prefix = "${meta.id}"
    
    // Assembly mode - hybrid if short reads available, long-only if not
    def assembly_mode = (short_reads && short_reads.size() >= 2) ? "hybrid-single" : "long-single"
    
    // // Parameters
    // def threads = task.cpus ?: 4
    // def memory = task.memory ? "${task.memory.toGiga()}G" : "16G"
    // def chromosome_size = params.chromosome_size ?: "auto"
    // def chr_size_param = chromosome_size == "auto" ? "" : "-c ${chromosome_size}"
    
    // // Optional parameters
    // def min_length = params.min_length ? "--min_length ${params.min_length}" : ""
    // def min_quality = params.min_quality ? "--min_quality ${params.min_quality}" : ""
    // def contaminants = params.contaminants ? "--contaminants ${params.contaminants}" : ""
    // def logic = params.logic ? "--logic ${params.logic}" : ""
    // def depth_filter = params.depth_filter ? "--depth_filter ${params.depth_filter}" : ""
    // def no_pypolca = params.no_pypolca ? "--no_pypolca" : ""
    // // def no_medaka = params.no_medaka ? "--no_medaka" : ""
    // // def medaka_model = params.medaka_model ? "--medaka_model ${params.medaka_model}" : ""

    """
    mkdir -p ${prefix}
    
    if [[ "${assembly_mode}" == "hybrid-single" ]]; then
        hybracter hybrid-single \\
            -l ${long_reads} \\
            -1 ${short_reads[0]} \\
            -2 ${short_reads[1]} \\
            -s ${prefix} \\
            --auto
    else
        hybracter long-single \\
            -l ${long_reads} \\
            -s ${prefix} \\
            --auto
    fi
    """
}
