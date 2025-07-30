process FILTLONG {
    tag "${reads.baseName}"
    publishDir "${params.outdir}/filtered_reads", mode: 'copy'
    
    container 'quay.io/biocontainers/filtlong:0.2.1--h9a82719_0'
    
    input:
    path reads
    path illumina_1, stageAs: 'illumina_1.fastq.gz'
    path illumina_2, stageAs: 'illumina_2.fastq.gz'  
    path assembly, stageAs: 'assembly.fasta'
    
    output:
    path "${reads.baseName}_filtered.fastq.gz", emit: filtered_reads
    path "${reads.baseName}_filtlong.log", emit: log
    
    script:
    """
    # Set up reference options
    ILLUMINA_REF=""
    if [[ "${illumina_1.name}" != "input.1" && "${illumina_2.name}" != "input.2" ]]; then
        ILLUMINA_REF="-1 illumina_1.fastq.gz -2 illumina_2.fastq.gz"
    fi
    
    ASSEMBLY_REF=""
    if [[ "${assembly.name}" != "input.3" ]]; then
        ASSEMBLY_REF="-a assembly.fasta"
    fi
    
    # Set up filtering parameters
    MIN_LENGTH=""
    if [[ "${params.min_length}" != "null" ]]; then
        MIN_LENGTH="--min_length ${params.min_length}"
    fi
    
    MAX_LENGTH=""
    if [[ "${params.max_length}" != "null" ]]; then
        MAX_LENGTH="--max_length ${params.max_length}"
    fi
    
    KEEP_PERCENT=""
    if [[ "${params.keep_percent}" != "null" ]]; then
        KEEP_PERCENT="--keep_percent ${params.keep_percent}"
    fi
    
    TARGET_BASES=""
    if [[ "${params.target_bases}" != "null" ]]; then
        TARGET_BASES="--target_bases ${params.target_bases}"
    fi
    
    MIN_MEAN_Q=""
    if [[ "${params.min_mean_q}" != "null" ]]; then
        MIN_MEAN_Q="--min_mean_q ${params.min_mean_q}"
    fi
    
    MIN_WINDOW_Q=""
    if [[ "${params.min_window_q}" != "null" ]]; then
        MIN_WINDOW_Q="--min_window_q ${params.min_window_q}"
    fi
    
    # Set up scoring weights
    LENGTH_WEIGHT=""
    if [[ "${params.length_weight}" != "1" ]]; then
        LENGTH_WEIGHT="--length_weight ${params.length_weight}"
    fi
    
    MEAN_Q_WEIGHT=""
    if [[ "${params.mean_q_weight}" != "1" ]]; then
        MEAN_Q_WEIGHT="--mean_q_weight ${params.mean_q_weight}"
    fi
    
    WINDOW_Q_WEIGHT=""
    if [[ "${params.window_q_weight}" != "1" ]]; then
        WINDOW_Q_WEIGHT="--window_q_weight ${params.window_q_weight}"
    fi
    
    # Set up read manipulation options
    TRIM=""
    if [[ "${params.trim}" == "true" ]]; then
        TRIM="--trim"
    fi
    
    SPLIT=""
    if [[ "${params.split}" != "null" ]]; then
        SPLIT="--split ${params.split}"
    fi
    
    WINDOW_SIZE=""
    if [[ "${params.window_size}" != "250" ]]; then
        WINDOW_SIZE="--window_size ${params.window_size}"
    fi
    
    VERBOSE=""
    if [[ "${params.verbose}" == "true" ]]; then
        VERBOSE="--verbose"
    fi
    
    echo "Starting Filtlong filtering for ${reads.baseName}" > ${reads.baseName}_filtlong.log
    echo "Input file: ${reads}" >> ${reads.baseName}_filtlong.log
    echo "Parameters used:" >> ${reads.baseName}_filtlong.log
    
    # Count input reads and bases
    echo "Input statistics:" >> ${reads.baseName}_filtlong.log
    if [[ ${reads} == *.gz ]]; then
        INPUT_READS=\$(zcat ${reads} | awk 'NR%4==1' | wc -l)
        INPUT_BASES=\$(zcat ${reads} | awk 'NR%4==2 {sum+=length(\$0)} END {print sum}')
    else
        INPUT_READS=\$(awk 'NR%4==1' ${reads} | wc -l)
        INPUT_BASES=\$(awk 'NR%4==2 {sum+=length(\$0)} END {print sum}' ${reads})
    fi
    echo "  Input reads: \$INPUT_READS" >> ${reads.baseName}_filtlong.log
    echo "  Input bases: \$INPUT_BASES" >> ${reads.baseName}_filtlong.log
    echo "" >> ${reads.baseName}_filtlong.log
    
    # Run Filtlong
    filtlong \\
        \$ILLUMINA_REF \\
        \$ASSEMBLY_REF \\
        \$MIN_LENGTH \\
        \$MAX_LENGTH \\
        \$KEEP_PERCENT \\
        \$TARGET_BASES \\
        \$MIN_MEAN_Q \\
        \$MIN_WINDOW_Q \\
        \$LENGTH_WEIGHT \\
        \$MEAN_Q_WEIGHT \\
        \$WINDOW_Q_WEIGHT \\
        \$TRIM \\
        \$SPLIT \\
        \$WINDOW_SIZE \\
        \$VERBOSE \\
        ${reads} \\
        2>> ${reads.baseName}_filtlong.log \\
        | gzip > ${reads.baseName}_filtered.fastq.gz
    
    # Count output reads and bases
    echo "" >> ${reads.baseName}_filtlong.log
    echo "Output statistics:" >> ${reads.baseName}_filtlong.log
    OUTPUT_READS=\$(zcat ${reads.baseName}_filtered.fastq.gz | awk 'NR%4==1' | wc -l)
    OUTPUT_BASES=\$(zcat ${reads.baseName}_filtered.fastq.gz | awk 'NR%4==2 {sum+=length(\$0)} END {print sum}')
    echo "  Output reads: \$OUTPUT_READS" >> ${reads.baseName}_filtlong.log  
    echo "  Output bases: \$OUTPUT_BASES" >> ${reads.baseName}_filtlong.log
    echo "  Reads retained: \$(echo "scale=2; \$OUTPUT_READS * 100 / \$INPUT_READS" | bc)%" >> ${reads.baseName}_filtlong.log
    echo "  Bases retained: \$(echo "scale=2; \$OUTPUT_BASES * 100 / \$INPUT_BASES" | bc)%" >> ${reads.baseName}_filtlong.log
    
    echo "Filtlong filtering completed successfully" >> ${reads.baseName}_filtlong.log
    """
}

/*
 * Process: Generate filtering summary report
 */
process FILTLONG_REPORT {
    publishDir "${params.outdir}/reports", mode: 'copy'
    
    container 'quay.io/biocontainers/python:3.9'
    
    input:
    path filtered_reads
    path log_file
    
    output:
    path "filtlong_summary.html"
    path "filtlong_stats.txt"
    
    script:
    """
    # Create summary statistics
    cat > filtlong_stats.txt << EOF
Filtlong Filtering Summary
=========================
Run Date: \$(date)
Input File: ${filtered_reads}

Parameters Used:
EOF
    
    # Extract parameters from command line
    echo "  Minimum length: ${params.min_length}" >> filtlong_stats.txt
    echo "  Keep percent: ${params.keep_percent}%" >> filtlong_stats.txt
    if [ "${params.max_length}" != "null" ]; then
        echo "  Maximum length: ${params.max_length}" >> filtlong_stats.txt
    fi
    if [ "${params.target_bases}" != "null" ]; then
        echo "  Target bases: ${params.target_bases}" >> filtlong_stats.txt
    fi
    if [ "${params.illumina_1}" != "null" ]; then
        echo "  Illumina reference: Yes" >> filtlong_stats.txt
    fi
    if [ "${params.assembly}" != "null" ]; then
        echo "  Assembly reference: Yes" >> filtlong_stats.txt
    fi
    if [ "${params.trim}" == "true" ]; then
        echo "  Trimming: Enabled" >> filtlong_stats.txt
    fi
    if [ "${params.split}" != "null" ]; then
        echo "  Splitting: ${params.split} bases" >> filtlong_stats.txt
    fi
    
    echo "" >> filtlong_stats.txt
    
    # Add statistics from log file
    grep -A 10 "Input statistics:" ${log_file} >> filtlong_stats.txt || echo "No input statistics found" >> filtlong_stats.txt
    grep -A 10 "Output statistics:" ${log_file} >> filtlong_stats.txt || echo "No output statistics found" >> filtlong_stats.txt
    
    # Create HTML report
    python3 << 'EOF'
import os
from datetime import datetime

# Read stats file
with open('filtlong_stats.txt', 'r') as f:
    stats_content = f.read()

# Create HTML report
html_content = f'''
<!DOCTYPE html>
<html>
<head>
    <title>Filtlong Filtering Report</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 40px; }}
        .header {{ background-color: #f0f0f0; padding: 20px; border-radius: 5px; }}
        .stats {{ background-color: #f8f8f8; padding: 15px; margin: 20px 0; border-radius: 5px; }}
        .parameter {{ margin: 5px 0; }}
        pre {{ background-color: #f5f5f5; padding: 10px; border-radius: 3px; }}
    </style>
</head>
<body>
    <div class="header">
        <h1>Filtlong Quality Filtering Report</h1>
        <p>Generated on: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}</p>
    </div>
    
    <div class="stats">
        <h2>Filtering Summary</h2>
        <pre>{stats_content}</pre>
    </div>
    
    <div class="stats">
        <h2>About Filtlong</h2>
        <p>Filtlong is a tool for filtering long reads by quality. It uses both read length 
        (longer is better) and read identity (higher is better) when choosing which reads 
        pass the filter.</p>
        
        <h3>Key Features:</h3>
        <ul>
            <li>Quality-based filtering using Phred scores or external reference</li>
            <li>Length-based filtering with customizable thresholds</li>
            <li>Percentage-based filtering to keep best reads</li>
            <li>Target base count filtering for consistent output sizes</li>
            <li>Read trimming and splitting for improved quality</li>
        </ul>
    </div>
</body>
</html>
'''

with open('filtlong_summary.html', 'w') as f:
    f.write(html_content)
EOF
    """
}