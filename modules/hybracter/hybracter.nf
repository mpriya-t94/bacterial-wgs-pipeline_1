process HYBRACTER {
    tag "$meta.id"
    label 'process_high'

    // Use the latest Hybracter container (v0.11.0)
    container 'biocontainers/hybracter:0.7.0--pyhdfd78af_0'


    input:
    tuple val(meta), path(long_reads), path(short_reads), val(mode), val(chromosome_size)

    output:
    tuple val(meta), path("*/FINAL_OUTPUT/complete/*.fasta"), emit: complete_assemblies
    tuple val(meta), path("*/FINAL_OUTPUT/incomplete/*.fasta"), emit: incomplete_assemblies, optional: true
    tuple val(meta), path("*/FINAL_OUTPUT/plasmids/*.fasta"), emit: plasmids, optional: true
    tuple val(meta), path("*/hybracter_summary.tsv"), emit: summary
    tuple val(meta), path("*/processing_summary.tsv"), emit: processing_summary
    tuple val(meta), path("*/FINAL_OUTPUT"), emit: final_output_dir
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def assembly_mode = mode ?: 'hybrid'
    def chr_size = chromosome_size ?: 'auto'
    def chr_size_param = chr_size == 'auto' ? '' : "--chromosome ${chr_size}"
    def threads = task.cpus ?: 4
    def memory = task.memory ? "${task.memory.toGiga()}G" : "16G"

    """
    # Create output directory
    mkdir -p ${prefix}

    # Run Hybracter based on mode
    case "${assembly_mode}" in
        "hybrid")
            hybracter hybrid \\
                -i ${long_reads} \\
                --out ${prefix} \\
                --threads ${threads} \\
                --ram ${memory} \\
                ${chr_size_param} \\
                ${args}
            ;;
        "hybrid-single")
            hybracter hybrid-single \\
                -l ${long_reads} \\
                -1 ${short_reads[0]} \\
                -2 ${short_reads.size() > 1 ? short_reads[1] : ''} \\
                -s ${meta.id} \\
                -c ${chr_size} \\
                --out ${prefix} \\
                --threads ${threads} \\
                --ram ${memory} \\
                ${args}
            ;;
        "long")
            hybracter long \\
                -i ${long_reads} \\
                --out ${prefix} \\
                --threads ${threads} \\
                --ram ${memory} \\
                ${chr_size_param} \\
                ${args}
            ;;
        "long-single")
            hybracter long-single \\
                -l ${long_reads} \\
                -s ${meta.id} \\
                -c ${chr_size} \\
                --out ${prefix} \\
                --threads ${threads} \\
                --ram ${memory} \\
                ${args}
            ;;
        *)
            echo "Error: Unsupported assembly mode: ${assembly_mode}"
            echo "Supported modes: hybrid, hybrid-single, long, long-single"
            exit 1
            ;;
    esac

    # Generate versions file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hybracter: \$(hybracter --version 2>&1 | grep -E '^hybracter' | sed 's/hybracter //')
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}/FINAL_OUTPUT/{complete,incomplete,plasmids}
    touch ${prefix}/FINAL_OUTPUT/complete/${prefix}_complete.fasta
    touch ${prefix}/FINAL_OUTPUT/plasmids/${prefix}_plasmid.fasta
    touch ${prefix}/hybracter_summary.tsv
    touch ${prefix}/processing_summary.tsv
    touch versions.yml
    """
}