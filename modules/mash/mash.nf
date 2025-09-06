process MASH_SKETCH {
    tag "$meta.id"
    label 'medium'
    publishDir "${params.outdir}/sketches/mash/individual", mode: 'copy', overwrite: true

    cpus params.cpus_mash
    memory params.memory_mash

    container 'quay.io/biocontainers/mash:2.3--he348c14_1'

    input:
    tuple val(meta), path(sequence_file)

    output:
    tuple val(meta), path("${meta.id}.msh"), emit: mash_sketch
    path "${meta.id}.msh", emit: sketch_files
    path "mash_sketch_${meta.id}.log", emit: mash_sketch_log
    path "${meta.id}.mash.txt", emit: mash_screen

    script:
    def min_copies_opt = params.mash_min_copies ? "-m ${params.mash_min_copies}" : ""
    """
    mash sketch \\
        -o ${meta.id} \\
        -k ${params.mash_kmer_size} \\
        -s ${params.mash_sketch_size} \\
        -g ${params.mash_genome_size} \\
        ${min_copies_opt} \\
        ${sequence_file} \\
        2>&1 | tee mash_sketch_${meta.id}.log
    """
}

process MASH_PASTE_SKETCHES {
    tag "mash_paste_sketches"
    label 'medium'
    publishDir "${params.outdir}/sketches", mode: 'copy', overwrite: true

    cpus params.cpus_mash
    memory params.memory_mash

    container 'quay.io/biocontainers/mash:2.3--he348c14_1'

    input:
    path sketch_files

    output:
    path "all_samples.msh", emit: all_samples_msh
    path "mash_paste_sketches.log", emit: mash_paste_log

    script:
    """
    mash paste \\
        -o all_samples \\
        ${sketch_files.join(' ')} \\
        2>&1 | tee mash_paste_sketches.log
    """
}

process MASH_DIST_MATRIX {
    tag "mash_dist_matrix"
    label 'high'

    publishDir "${params.outdir}/distances", mode: 'copy', overwrite: true

    cpus params.cpus_mash
    memory params.memory_mash

    container 'quay.io/biocontainers/mash:2.3--he348c14_1'

    input:
    path all_samples_msh

    output:
    path "distance_matrix.tsv", emit: distance_matrix_tsv
    path "similarity_matrix.tsv", emit: similarity_score_tsv
    path "mash_dist_matrix.log", emit: mash_dist_log

    script:
    def max_dist_opt = params.mash_max_distance < 1.0 ? "-d ${params.mash_max_distance}" : ""
    def pvalue_opt = params.pvalue_threshold ? "-v ${params.pvalue_threshold}" : ""
    """
    # Compute pairwise distance matrix
    mash dist \\
        -t \\
        ${max_dist_opt} \\
        ${pvalue_opt} \\
        ${all_samples_msh} \\
        ${all_samples_msh} \\
        > distance_matrix.tsv \\
        2>&1 | tee mash_dist_matrix.log
        
    # Compute similarity score matrix (1 - distance)
    awk 'BEGIN{OFS="\\t"} 
         NR==1{print \$0} 
         NR>1{for(i=2;i<=NF;i++) \$i=1-\$i; print \$0}' \\
        distance_matrix.tsv > similarity_matrix.tsv
    """
}

process MASH_SCREEN {
    tag "$meta.id"
    label 'medium'

    publishDir "${params.outdir}/screen", mode: 'copy', overwrite: true

    cpus params.cpus_mash
    memory params.memory_mash

    container 'quay.io/biocontainers/mash:2.3--he348c14_1'

    when:
    params.reference_sketch != null

    input:
    tuple val(meta), path(sequence_file)
    path reference_sketch

    output:
    tuple val(meta), path("${meta.id}_screen_results.tsv"), emit: mash_screen
    path "mash_screen_${meta.id}.log", emit: mash_screen_log

    script:
    def pvalue_opt = params.pvalue_threshold ? "-p ${params.pvalue_threshold}" : ""
    def winner_opt = params.winner_takes_all ? "-w" : ""
    def identity_opt = params.min_identity ? "-i ${params.min_identity}" : ""
    """
    mash screen \\
        ${winner_opt} \\
        ${identity_opt} \\
        ${pvalue_opt} \\
        ${reference_sketch} \\
        ${sequence_file} \\
        > ${meta.id}_screen_results.tsv \\
        2>&1 | tee mash_screen_${meta.id}.log
    """
}

process MASHTREE_BASIC {
    tag "mashtree_basic"
    label 'high'

    publishDir "${params.outdir}/trees", mode: 'copy', overwrite: true

    cpus params.cpus_mash
    memory params.memory_mash

    container 'quay.io/biocontainers/mashtree:1.4.6--pl5321hdfd78af_0'

    input:
    path sequence_files

    output:
    path "mashtree_basic.dnd", emit: mashtree_treefile
    path "mashtree_basic.log", emit: mashtree_log

    script:
    def mindepth_opt = params.min_depth == 0 ? "--mindepth 0" : "--mindepth ${params.min_depth}"
    """
    mashtree \\
        --numcpus ${task.cpus} \\
        ${mindepth_opt} \\
        --truncLength ${params.truncate_length} \\
        --sort-order ${params.sort_order} \\
        --genomesize ${params.mash_genome_size} \\
        --kmerlength ${params.mash_kmer_size} \\
        --sketch-size ${params.mash_sketch_size} \\
        ${sequence_files.join(' ')} \\
        > mashtree_basic.dnd \\
        2>&1 | tee mashtree_basic.log
    """
}

process MASHTREE_WITH_MATRIX {
    tag "mashtree_with_matrix"
    label 'high'

    publishDir "${params.outdir}/trees", mode: 'copy', overwrite: true

    cpus params.cpus_mash
    memory params.memory_mash

    container 'quay.io/biocontainers/mashtree:1.4.6--pl5321hdfd78af_0'

    input:
    path sequence_files
    path distance_matrix_tsv

    output:
    path "mashtree_with_matrix.dnd", emit: mashtree_with_matrix_treefile
    path "mashtree_with_matrix.tsv", emit: mashtree_with_matrix_tsv
    path "mashtree_with_matrix.log", emit: mashtree_with_matrix_log

    script:
    def mindepth_opt = params.min_depth == 0 ? "--mindepth 0" : "--mindepth ${params.min_depth}"
    """
    mashtree \\
        --numcpus ${task.cpus} \\
        ${mindepth_opt} \\
        --truncLength ${params.truncate_length} \\
        --sort-order ${params.sort_order} \\
        --genomesize ${params.mash_genome_size} \\
        --kmerlength ${params.mash_kmer_size} \\
        --sketch-size ${params.mash_sketch_size} \\
        --distancematrix ${distance_matrix_tsv} \\
        ${sequence_files.join(' ')} \\
        > mashtree_with_matrix.dnd \\
        2>&1 | tee mashtree_with_matrix.log
    
    # Copy the input distance matrix as output for downstream processes
    cp ${distance_matrix_tsv} mashtree_with_matrix.tsv
    """
}

process MASHTREE_BOOTSTRAP {
    tag "mashtree_bootstrap"
    label 'high'

    publishDir "${params.outdir}/trees", mode: 'copy', overwrite: true

    cpus params.cpus_mash
    memory params.memory_mash

    container 'quay.io/biocontainers/mashtree:1.4.6--pl5321hdfd78af_0'

    when:
    params.confidence_method == 'bootstrap'

    input:
    path sequence_files

    output:
    path "mashtree_bootstrap.dnd", emit: mashtree_bootstrap_treefile
    path "mashtree_bootstrap.log", emit: mashtree_bootstrap_log

    script:
    def mindepth_opt = params.min_depth == 0 ? "--mindepth 0" : "--mindepth ${params.min_depth}"
    """
    mashtree_bootstrap.pl \\
        --reps ${params.bootstrap_reps} \\
        --numcpus ${task.cpus} \\
        ${sequence_files.join(' ')} \\
        -- \\
        ${mindepth_opt} \\
        --truncLength ${params.truncate_length} \\
        --sort-order ${params.sort_order} \\
        --genomesize ${params.mash_genome_size} \\
        --kmerlength ${params.mash_kmer_size} \\
        --sketch-size ${params.mash_sketch_size} \\
        > mashtree_bootstrap.dnd \\
        2>&1 | tee mashtree_bootstrap.log
    """
}

process MASHTREE_JACKKNIFE {
    tag "mashtree_jackknife"
    label 'high'

    publishDir "${params.outdir}/trees", mode: 'copy', overwrite: true

    cpus params.cpus_mash
    memory params.memory_mash

    container 'quay.io/biocontainers/mashtree:1.4.6--pl5321hdfd78af_0'

    when:
    params.confidence_method == 'jackknife'

    input:
    path sequence_files

    output:
    path "mashtree_jackknife.dnd", emit: mashtree_jackknife_treefile
    path "mashtree_jackknife.log", emit: mashtree_jackknife_log

    script:
    def mindepth_opt = params.min_depth == 0 ? "--mindepth 0" : "--mindepth ${params.min_depth}"
    """
    mashtree_jackknife.pl \\
        --reps ${params.jackknife_reps} \\
        --numcpus ${task.cpus} \\
        ${sequence_files.join(' ')} \\
        -- \\
        ${mindepth_opt} \\
        --truncLength ${params.truncate_length} \\
        --sort-order ${params.sort_order} \\
        --genomesize ${params.mash_genome_size} \\
        --kmerlength ${params.mash_kmer_size} \\
        --sketch-size ${params.mash_sketch_size} \\
        > mashtree_jackknife.dnd \\
        2>&1 | tee mashtree_jackknife.log
    """
}

process COMBINE_SCREEN_RESULTS {
    tag "combine_screen_results"
    label 'low'

    publishDir "${params.outdir}/screen", mode: 'copy', overwrite: true

    cpus 1
    memory '2 GB'

    container 'quay.io/biocontainers/csvtk:0.25.0--h9ee0642_0'

    when:
    params.reference_sketch != null

    input:
    path mash_screen_files

    output:
    path "combined_mash_screen_results.tsv", emit: combined_mash_screen_results
    path "combine_screen_results.log", emit: combine_screen_log
    path "screen_summary.csv", emit: screen_summary_csv
    
    script:
    """
    # Get the first file to extract header
    first_file=\$(echo ${mash_screen_files.join(' ')} | cut -d' ' -f1)
    
    # Extract header from first file
    head -1 \$first_file > combined_mash_screen_results.tsv

    # Combine all files (skip headers except from first file)
    for file in ${mash_screen_files.join(' ')}; do
        if [ "\$file" = "\$first_file" ]; then
            # Include entire first file
            cat \$file >> combined_mash_screen_results.tsv
        else
            # Skip header for subsequent files
            tail -n +2 \$file >> combined_mash_screen_results.tsv
        fi
    done

    # Create summary statistics
    total_lines=\$(wc -l < combined_mash_screen_results.tsv)
    num_files=\$(echo ${mash_screen_files.join(' ')} | wc -w)
    
    echo "Sample,Results_Count" > screen_summary.csv
    for file in ${mash_screen_files.join(' ')}; do
        sample_name=\$(basename \$file _screen_results.tsv)
        result_count=\$(tail -n +2 \$file | wc -l)
        echo "\$sample_name,\$result_count" >> screen_summary.csv
    done

    echo "Combined \$total_lines lines from \$num_files files." | tee combine_screen_results.log
    echo "Created summary with per-sample result counts." | tee -a combine_screen_results.log
    """
}

process CREATE_SUMMARY_REPORT {
    tag "create_summary_report"
    label 'low'

    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    cpus 1
    memory '2 GB'

    container 'python:3.9-slim'

    input:
    path distance_matrix, stageAs: "mash_distance_matrix.tsv"
    path similarity_scores, stageAs: "mash_similarity_scores.tsv"
    path mashtree_matrix, stageAs: "mashtree_distance_matrix.tsv", optional: true
    path basic_tree, stageAs: "basic_tree.dnd", optional: true
    path tree_with_matrix, stageAs: "tree_with_matrix.dnd", optional: true
    path bootstrap_tree, stageAs: "bootstrap_tree.dnd", optional: true
    path jackknife_tree, stageAs: "jackknife_tree.dnd", optional: true
    path combined_screen, stageAs: "combined_screen_results.tsv", optional: true
    
    output:
    path "pipeline_summary_report.html", emit: html_report
    path "pipeline_summary.txt", emit: text_summary
    path "analysis_stats.json", emit: stats_json
    
    script:
    """
    cat > create_report.py << 'EOF'
import os
import json
from pathlib import Path
import datetime

def analyze_files():
    results = {
        "pipeline": "${params.mode ?: 'mash_analysis'}",
        "timestamp": datetime.datetime.now().isoformat(),
        "parameters": {
            "sketch_size": ${params.mash_sketch_size ?: 1000},
            "kmer_size": ${params.mash_kmer_size ?: 21},
            "genome_size": ${params.mash_genome_size ?: '5000000'},
            "min_depth": ${params.min_depth ?: 5},
            "confidence_method": "${params.confidence_method ?: 'none'}",
        },
        "outputs": {}
    }

    # Check which files exist and analyze them
    files_to_check = [
        ("mash_distance_matrix.tsv", "mash_distance_matrix"),
        ("mash_similarity_scores.tsv", "mash_similarity_scores"),
        ("mashtree_distance_matrix.tsv", "mashtree_distance_matrix"),
        ("basic_tree.dnd", "basic_tree"),
        ("tree_with_matrix.dnd", "tree_with_matrix"),
        ("bootstrap_tree.dnd", "bootstrap_tree"),
        ("jackknife_tree.dnd", "jackknife_tree"),
        ("combined_screen_results.tsv", "combined_screen_results"),
    ]

    for filename, description in files_to_check:
        if os.path.exists(filename):
            file_size = os.path.getsize(filename)
            results["outputs"][description] = {
                "description": description,
                "size_bytes": file_size,
                "exists": True
            }
            
            # Count lines for text files
            if filename.endswith(('.tsv', '.dnd', '.csv')):
                try:
                    with open(filename, 'r') as f:
                        line_count = sum(1 for line in f)
                    results["outputs"][description]["line_count"] = line_count
                except Exception as e:
                    results["outputs"][description]["error"] = str(e)
        else:
            results["outputs"][description] = {
                "description": description,
                "exists": False
            }
    
    return results

def create_html_report(data):
    html_content = f'''
    <!DOCTYPE html>
    <html>
    <head>
        <title>MASH Pipeline Summary Report</title>
        <style>
            body {{ font-family: Arial, sans-serif; margin: 40px; }}
            h1, h2 {{ color: #333; }}
            table {{ border-collapse: collapse; width: 100%; }}
            th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
            th {{ background-color: #f2f2f2; }}
            .exists {{ color: green; }}
            .missing {{ color: red; }}
        </style>
    </head>
    <body>
        <h1>MASH Pipeline Summary Report</h1>
        <p><strong>Generated:</strong> {data["timestamp"]}</p>
        <p><strong>Pipeline Mode:</strong> {data["pipeline"]}</p>
        
        <h2>Parameters</h2>
        <table>
            <tr><th>Parameter</th><th>Value</th></tr>
    '''
    
    for param, value in data["parameters"].items():
        html_content += f'<tr><td>{param}</td><td>{value}</td></tr>'
    
    html_content += '''
        </table>
        
        <h2>Output Files</h2>
        <table>
            <tr><th>File</th><th>Status</th><th>Size (bytes)</th><th>Lines</th></tr>
    '''
    
    for output_name, output_data in data["outputs"].items():
        status = "exists" if output_data["exists"] else "missing"
        status_class = "exists" if output_data["exists"] else "missing"
        size = output_data.get("size_bytes", "N/A")
        lines = output_data.get("line_count", "N/A")
        
        html_content += f'''
        <tr>
            <td>{output_name}</td>
            <td class="{status_class}">{status}</td>
            <td>{size}</td>
            <td>{lines}</td>
        </tr>
        '''
    
    html_content += '''
        </table>
    </body>
    </html>
    '''
    
    return html_content

def create_text_summary(data):
    summary = f"""MASH Pipeline Summary Report
Generated: {data["timestamp"]}
Pipeline Mode: {data["pipeline"]}

Parameters:
"""
    for param, value in data["parameters"].items():
        summary += f"  {param}: {value}\\n"
    
    summary += "\\nOutput Files:\\n"
    for output_name, output_data in data["outputs"].items():
        status = "EXISTS" if output_data["exists"] else "MISSING"
        size = output_data.get("size_bytes", "N/A")
        lines = output_data.get("line_count", "N/A")
        summary += f"  {output_name}: {status} (Size: {size} bytes, Lines: {lines})\\n"
    
    return summary

# Main execution
data = analyze_files()

# Write JSON stats
with open("analysis_stats.json", "w") as f:
    json.dump(data, f, indent=2)

# Write HTML report
html_report = create_html_report(data)
with open("pipeline_summary_report.html", "w") as f:
    f.write(html_report)

# Write text summary
text_summary = create_text_summary(data)
with open("pipeline_summary.txt", "w") as f:
    f.write(text_summary)

print("Report generation completed successfully!")
EOF

    python create_report.py
    """
}