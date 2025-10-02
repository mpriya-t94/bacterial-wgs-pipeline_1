process MASH_SKETCH {
    tag "$meta.id"
    label 'medium'
    publishDir "${params.outdir}/comparative/${analysis_type}/sketches", mode: 'copy'

    container 'quay.io/biocontainers/mash:2.3--he348c14_1'

    input:
    tuple val(meta), path(assembly), val(analysis_type)

    output:
    tuple val(meta), path("${meta.id}.msh"), val(analysis_type), emit: mash_sketch
    path "${meta.id}.msh", emit: sketch_file

    script:
    """
    mash sketch \\
        -o ${meta.id} \\
        -k ${params.mash_kmer_size} \\
        -s ${params.mash_sketch_size} \\
        ${assembly}
    """
}

process MASH_PASTE_SKETCHES {
    tag "${analysis_type}"
    label 'medium'
    publishDir "${params.outdir}/comparative/${analysis_type}", mode: 'copy'

    container 'quay.io/biocontainers/mash:2.3--he348c14_1'

    input:
    tuple val(analysis_type), path(sketch_files)

    output:
    tuple val(analysis_type), path("${analysis_type}_combined.msh"), emit: combined_sketch

    script:
    """
    mash paste \\
        ${analysis_type}_combined \\
        ${sketch_files.join(' ')}
    """
}

process MASH_DIST_MATRIX {
    tag "${analysis_type}"
    label 'high'
    publishDir "${params.outdir}/comparative/${analysis_type}", mode: 'copy'

    container 'quay.io/biocontainers/mash:2.3--he348c14_1'

    input:
    tuple val(analysis_type), path(combined_sketch)

    output:
    tuple val(analysis_type), path("${analysis_type}_distance_matrix.tsv"), emit: distance_matrix
    tuple val(analysis_type), path("${analysis_type}_similarity_matrix.tsv"), emit: similarity_matrix

    script:
    """
    mash dist \\
        -t \\
        ${combined_sketch} \\
        ${combined_sketch} \\
        > ${analysis_type}_distance_matrix.tsv
        
    awk 'BEGIN{OFS="\\t"} 
         NR==1{print \$0} 
         NR>1{for(i=2;i<=NF;i++) \$i=1-\$i; print \$0}' \\
        ${analysis_type}_distance_matrix.tsv > ${analysis_type}_similarity_matrix.tsv
    """
}

process MASHTREE_WITH_MATRIX {
    tag "${analysis_type}"
    label 'high'
    publishDir "${params.outdir}/comparative/${analysis_type}", mode: 'copy'

    container 'quay.io/biocontainers/mashtree:1.4.6--pl5321hdfd78af_0'

    input:
    tuple val(analysis_type), path(distance_matrix), path(assemblies)

    output:
    tuple val(analysis_type), path("${analysis_type}_tree.dnd"), emit: tree

    script:
    """
    mashtree \\
        --numcpus ${task.cpus} \\
        --truncLength ${params.truncate_length} \\
        --sort-order ${params.sort_order} \\
        --genomesize ${params.mash_genome_size} \\
        --kmerlength ${params.mash_kmer_size} \\
        --sketch-size ${params.mash_sketch_size} \\
        ${assemblies.join(' ')} \\
        > ${analysis_type}_tree.dnd
    """
}
