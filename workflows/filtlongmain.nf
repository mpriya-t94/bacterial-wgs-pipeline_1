nextflow.enable.dsl = 2

// Include modules
include { FILTLONG; FILTLONG_REPORT } from '../modules/filtlong/filtlong.nf'

params.input = null                    // Input long reads (FASTQ/FASTA)
params.outdir = "results"             // Output directory
params.illumina_1 = null              // Illumina R1 reads (optional reference)
params.illumina_2 = null              // Illumina R2 reads (optional reference)
params.assembly = null                // Reference assembly (optional)

// Filtlong filtering parameters
params.min_length = 1000              // Minimum read length
params.max_length = null              // Maximum read length  
params.keep_percent = 90              // Keep this percentage of best reads
params.target_bases = null            // Target number of bases to keep
params.min_mean_q = null              // Minimum mean quality threshold
params.min_window_q = null            // Minimum window quality threshold

// Filtlong scoring weights
params.length_weight = 1              // Weight for length score
params.mean_q_weight = 1              // Weight for mean quality score
params.window_q_weight = 1            // Weight for window quality score

// Read manipulation options
params.trim = false                   // Trim non-matching bases from ends
params.split = null                   // Split reads at N consecutive non-matching bases
params.window_size = 250              // Window size for quality assessment

// Other options
params.verbose = false                // Verbose output
params.help = false

workflow {
    // Show help message if requested
    if (params.help) {
        log.info"""
        Usage:
            nextflow run main.nf --input <long_reads.fastq.gz> [options]

        Required arguments:
            --input                 Path to input long reads (FASTQ/FASTA format)

        Output options:
            --outdir                Output directory (default: results)

        External reference options (choose one):
            --illumina_1            Illumina R1 reads for reference
            --illumina_2            Illumina R2 reads for reference  
            --assembly              Reference assembly in FASTA format

        Filtering thresholds:
            --min_length            Minimum read length (default: 1000)
            --max_length            Maximum read length
            --keep_percent          Keep this percentage of best reads (default: 90)
            --target_bases          Target number of bases to keep
            --min_mean_q            Minimum mean quality threshold
            --min_window_q          Minimum window quality threshold

        Scoring weights:
            --length_weight         Weight for length score (default: 1)
            --mean_q_weight         Weight for mean quality score (default: 1)
            --window_q_weight       Weight for window quality score (default: 1)

        Read manipulation:
            --trim                  Trim non-matching bases from read ends
            --split                 Split reads at N consecutive non-matching bases
            --window_size           Window size for quality assessment (default: 250)

        Other:
            --verbose               Enable verbose output
            --help                  Show this help message

        Examples:
            # Basic filtering by length and percentage
            nextflow run main.nf --input long_reads.fastq.gz

            # Using Illumina reads as reference with trimming/splitting
            nextflow run main.nf --input long_reads.fastq.gz \\
                --illumina_1 illumina_R1.fastq.gz \\
                --illumina_2 illumina_R2.fastq.gz \\
                --trim --split 500

            # Using assembly reference with custom parameters  
            nextflow run main.nf --input long_reads.fastq.gz \\
                --assembly reference.fasta \\
                --min_length 2000 --keep_percent 85 \\
                --length_weight 2
        """.stripIndent()
        exit 0
    }

    // Validate required parameters
    if (!params.input) {
        log.error "Error: --input parameter is required"
        log.info "Use --help to see usage information"
        exit 1
    }
    
    // Create input channel for long reads
    input_reads = Channel.fromPath(params.input, checkIfExists: true)
    
    // Create reference channels if provided
    illumina_1_ch = params.illumina_1 ? Channel.fromPath(params.illumina_1, checkIfExists: true) : Channel.empty()
    illumina_2_ch = params.illumina_2 ? Channel.fromPath(params.illumina_2, checkIfExists: true) : Channel.empty()
    assembly_ch = params.assembly ? Channel.fromPath(params.assembly, checkIfExists: true) : Channel.empty()
    
    // Run Filtlong filtering
    FILTLONG(
        input_reads,
        illumina_1_ch,
        illumina_2_ch, 
        assembly_ch
    )
    
    // Generate summary report
    FILTLONG_REPORT(FILTLONG.out.filtered_reads, FILTLONG.out.log)
    
    // Workflow completion and error handling
    workflow.onComplete {
        log.info """
        Pipeline completed!
        Results saved to: ${params.outdir}
        
        Output files:
        - Filtered reads: ${params.outdir}/filtered_reads/
        - Summary report: ${params.outdir}/reports/filtlong_summary.html
        - Statistics: ${params.outdir}/reports/filtlong_stats.txt
        """.stripIndent()
    }

    workflow.onError {
        log.error "Pipeline failed with error: ${workflow.errorMessage}"
    }
}