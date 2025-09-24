process AMRFINDERPLUS {
    container 'quay.io/biocontainers/ncbi-amrfinderplus:3.12.8--h283d18e_0'
    publishDir params.outdir, mode: 'copy'

    input:
    path genomes

    output:
    path "amrfinderplus_output", emit: amrfinderplus_results
    path "amrfinderplus_output/amrfinderplus_report.csv", emit: amrfinderplus_csv
    path "amrfinderplus_output/combined_results.csv", emit: amrfinderplus_combined_csv
    path "amrfinderplus_output/pointmutations.csv", emit: amrfinderplus_pointmutations_csv

    script:
    def organims_option = params.amrfinder_organism ? "--organism ${params.amrfinder_organism}" : ""
    def plus_option = params.amrfinder_plus ? "--plus" : ""
    def ident_min = params.amrfinder_ident_min ?: 0.9
    def coverage_min = params.amrfinder_coverage_min ?: 0.5
    def translate_table = params.amrfinder_translate_table ?: 11
    def database_version = params.amrfinder_db_version ? "--database_version ${params.amrfinder_db_version}" : ""

    """
    # create output directory
    mkdir -p amrfinderplus_output

    # create input directory
    mkdir -p genomes_dir

    # copy all genomes files to input directory
    for files in ${genomes}; do
        if [ -f "\$file" ]; then
            cp "\$file" genomes_dir/
        fi
        
        # Update AMRFinderPlus database (if needed)
        echo "Updating AMRFinderPlus database..."
        amrfinder_update ${database_version} || echo "Database update failed or not needed"
    
        # Display database information
        echo "AMRFinderPlus database information:"
        amrfinder --database_version 2>/dev/null || echo "Database version information not available"
    
        # List supported organisms if organism-specific analysis is requested
        if [ -n "${params.amrfinder_organism}" ]; then
            echo "Supported organisms for --organism option:"
            amrfinder --list_organisms 2>/dev/null || echo "Organism list not available"
            echo ""
        fi

        # Process each genome individually
        for genome in genomes_dir/*; do
            if [ -f "\$genome" ]; then
                basename_genome=\$(basename "\$genome")
                sample_name=\${basename_genome%.*}  # Extract sample name without extension

                echo "Processing \$basename_genome..."

                # Run AMRFinderPlus
                amrfinder --nucleotide "\$genome" \\
                    ${organims_option} \\
                    ${plus_option} \\
                    --ident_min ${ident_min} \\
                    --coverage_min ${coverage_min} \\
                    --translation_table ${translation_table} \\
                    --threads ${task.cpus ?: 1} \\
                    --name \$sample_name \\
                    --output "amrfinderplus_output/\${sample_name}_amrfinder.csv" \\
                    2> "amrfinderplus_output/\${sample_name}_amrfinder.log"

                # Check if output was generated
            if [ ! -f "amrfinderplus_output/\${sample_name}_amrfinder.tsv" ]; then
                echo "Warning: No output generated for \$sample_name"
                # Create empty file with headers
                echo -e "Protein identifier\\tContig id\\tStart\\tStop\\tStrand\\tGene symbol\\tSequence name\\tScope\\tElement type\\tElement subtype\\tClass\\tSubclass\\tMethod\\tTarget length\\tReference sequence length\\tCoverage\\tIdentity\\tAlignment length\\tAccession of closest sequence\\tName of closest sequence\\tHMM id\\tHMM description" > "amrfinderplus_output/\${sample_name}_amrfinder.tsv"
            fi
        fi
    done

    echo "AMRFinderPlus analysis completed."
    """

}