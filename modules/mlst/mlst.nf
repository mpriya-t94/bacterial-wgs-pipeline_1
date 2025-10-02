process MLST {
    container 'quay.io/biocontainers/mlst:2.23.0--hdfd78af_0'
    punlishDir params.outdir, mode: 'copy'

    input:
    path genomes

    output:
    path "mlst_output", emit: mlst_results
    path "mlst_output/mlst_results.csv", emit: mlst_json, optional: true
    path "mlst_output/mlst_results.json", emmit: mlst_json, optional: true
    path "mlst_output/novel_alleles.fasta", emit: novel_alleles, optional: true

    script:
    def scheme_option = params.mlst_scheme ? " --scheme ${params.mlst_scheme}" : ""
    def legacy_option = params.mlst_legacy ? " --legacy" : ""
    def minscore_option = params.mlst_minscore ? " --minscore ${params.mlst_minscore}" : ""
    def minid_option = params.mlst_minid ? " --minid ${params.mlst_minid}" : ""
    def mincov_option = params.mlst_mincov ? " --mincov ${params.mlst_mincov}" : ""
    
    """
    # Create output directory
    mkdir -p mlst_output

    # Create input directory
    mkdir -p genome_dir

    # Copy genome files to input directory
    for files in ${genomes} ; do
        if [ -f "\$files" ]; the
            cp "\$file" genome_dir/
        fi
    done

    # Run MLST on all genomes
    echo "Processing genomes with MLST..."

    mlst ${scheme_option} \\
        ${legacy_option} \\
        ${minscore_option} \\
        ${minid_option} \\
        ${mincov_option} \\
        --nopath \\
        --threads ${task.cpus ?: 1} \\
        genome_dir/* > mlst_output/mlst_results.csv

    # Generate JSON output if requested
    if [ "${params.mlst_output_json}" = "true" ]; then
        echo "Generating JSON output..."
        mlst ${scheme_option} \\
            ${minscore_option} \\
            ${minid_option} \\
            ${mincov_option} \\
            --json mlst_outptut/mlst_results.json \\
            --quiet \\
            --nopath \\
            --threads ${task.cpus ?: 1} \\
            genome_dir/*
    fi

    # Extract novel alleles if requested
    if [ "${params.mlst_extract_novel_alleles}" = "true" ]; then
        echo "Extracting novel alleles..."
        mlst ${scheme_option} \\
            ${minscore_option} \\
            ${minid_option} \\
            ${mincov_option} \\
            --novel_alleles mlst_output/novel_alleles.fasta \\
            --quiet \\
            --nopath \\
            --threads ${task.cpus ?: 1} \\
            genome_dir/*
    fi

        # Check if novel alleles file was created and has content
        if [ ! -s mlst_output/novel_alleles.fasta ]; then
            rm -f mlst_output/novel_alleles.fasta
        fi
    fi
    
    # Generate summary statistics
    echo "Generating summary statistics..."
    {
        echo "=== MLST Analysis Summary ==="
        echo "Date: \$(date)"
        echo "Total genomes processed: \$(wc -l < mlst_output/mlst_results.tsv | tail -1)"
        echo ""
        echo "Scheme distribution:"
        tail -n +2 mlst_output/mlst_results.tsv | cut -f2 | sort | uniq -c | sort -nr
        echo ""
        echo "ST distribution (top 10):"
        tail -n +2 mlst_output/mlst_results.tsv | cut -f3 | grep -v '^-\$' | sort | uniq -c | sort -nr | head -10
        echo ""
        if [ -f mlst_output/novel_alleles.fasta ]; then
            echo "Novel alleles found: \$(grep -c '^>' mlst_output/novel_alleles.fasta)"
        else
            echo "Novel alleles found: 0"
        fi
    } > mlst_output/summary_stats.txt
    
    echo "MLST analysis completed successfully"
    """
}