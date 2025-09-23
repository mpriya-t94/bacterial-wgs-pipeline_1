process ABRICATE {
    container 'quay.io/biocontainers/abricate:1.0.1--pl526hdfd78af_0'
    publishDir params.outdir, mode: 'copy'
    
    input:
    path genomes
    
    output:
    path "abricate_output", emit: abricate_results
    path "abricate_output/summary.txt", emit: abricate_summary
    
    script:
    """
    # Create output directory
    mkdir -p abricate_output
    
    # Create input directory
    mkdir -p genome_dir
    
    # Copy all genome files to input directory
    for file in ${genomes}; do
        if [ -f "\$file" ]; then
            cp "\$file" genome_dir/
        fi
    done
    
    # Validate database selection
    DB=${params.abricate_db ?: 'ncbi'}
    echo "Using ABRICATE database: \$DB"
    
    # Run ABRICATE on genome files
    # Handle multiple possible file extensions
    for genome in genome_dir/*; do
        if [ -f "\$genome" ]; then
            basename_genome=\$(basename "\$genome")
            # Remove file extension for cleaner output names
            sample_name=\${basename_genome%.*}
            
            echo "Processing \$basename_genome with database \$DB"
            
            abricate --threads ${task.cpus ?: params.threads ?: 1} \\
                     --db \$DB \\
                     --minid ${params.abricate_minid ?: 80} \\
                     --mincov ${params.abricate_mincov ?: 80} \\
                     "\$genome" > "abricate_output/\${sample_name}_abricate.tsv"
        fi
    done
    
    # Generate summary report
    if ls abricate_output/*.tsv 1> /dev/null 2>&1; then
        echo "Generating ABRICATE summary report..."
        abricate --summary abricate_output/*_abricate.tsv > abricate_output/summary.txt
    else
        echo "No ABRICATE results found to summarize" > abricate_output/summary.txt
    fi
    
    # Optional: Create a combined results file with headers
    if ls abricate_output/*_abricate.tsv 1> /dev/null 2>&1; then
        echo "Creating combined results file..."
        # Get header from first file
        head -1 abricate_output/*_abricate.tsv | head -1 > abricate_output/combined_results.tsv
        
        # Append all results without headers
        for result in abricate_output/*_abricate.tsv; do
            tail -n +2 "\$result" >> abricate_output/combined_results.tsv
        done
    fi
    """
}