nextflow.enable.dsl=2

// Parameters
params.genomes = './genomes/*.fasta'
params.outdir = './checkm2_results'
params.threads = 4

// Main workflow
workflow {
    genome_ch = Channel.fromPath(params.genomes)
    CHECKM2(genome_ch.collect())
}

// CheckM2 process
process CHECKM2 {
    container 'quay.io/biocontainers/checkm2:1.0.2--pyh7cba7a3_0'
    publishDir params.outdir, mode: 'copy'
    
    input:
    path genomes
    
    output:
    path "checkm2_output"
    path "checkm2_output/quality_report.tsv"
    
    script:
    """
    mkdir -p genome_dir
    
    # Copy all genome files to input directory
    for file in ${genomes}; do
        cp \$file genome_dir/
    done
    
    # Download database if not present
    export CHECKM2DB=\${PWD}/checkm2_db
    checkm2 database --download --path \${CHECKM2DB} || echo "Database already exists or download failed"
    
    # Run CheckM2
    checkm2 predict --threads ${params.threads} --input genome_dir --output-directory checkm2_output
    """
}