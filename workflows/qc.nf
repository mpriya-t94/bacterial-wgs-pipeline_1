// Import modules
include { FASTQC as FASTQC_PRE } from '../modules/fastqc/fastqc.nf'
include { FASTQC as FASTQC_POST } from '../modules/fastqc/fastqc.nf'
include { LONGQC as LONGQC_PRE } from '../modules/longqc/longqc.nf'
include { LONGQC as LONGQC_POST } from '../modules/longqc/longqc.nf'
include { FASTP } from '../modules/fastp/fastp.nf'
include { MULTIQC as MULTIQC_PRE } from '../modules/multiqc/multiqc.nf'
include { MULTIQC as MULTIQC_POST } from '../modules/multiqc/multiqc.nf'

//
// PRE-QC SUBWORKFLOW
//
workflow PRE_QC {
    take:
    short_reads  // Channel: [[meta], [reads]]
    long_reads   // Channel: [[meta], [reads]]

    main:
    
    // Run FastQC on short reads in parallel
    FASTQC_PRE(short_reads)
    
    // Run LongQC on long reads in parallel (only if long reads exist)
    ch_long_filtered = long_reads
        .filter { _meta, reads -> reads != null && reads.toString() != 'null' }
    
    LONGQC_PRE(ch_long_filtered)
    ch_longqc_results = LONGQC_PRE.out
    
    // Collect FastQC results for MultiQC (short reads only)
    ch_pre_multiqc = FASTQC_PRE.out.zip
        .map { _meta, files -> files }
        .flatten()
        .collect()
    
    // Run Pre-QC MultiQC report
    MULTIQC_PRE(ch_pre_multiqc, "pre")

    emit:
    // FastQC outputs
    fastqc_html = FASTQC_PRE.out.html
    fastqc_zip = FASTQC_PRE.out.zip
    
    // LongQC outputs
    longqc_html = ch_longqc_results.html
    longqc_plots = ch_longqc_results.plots
    
    // MultiQC outputs
    multiqc_html = MULTIQC_PRE.out.html
    multiqc_data = MULTIQC_PRE.out.data
    
    // Pass through channels for downstream processing
    short_reads_passthrough = short_reads
    long_reads_passthrough = long_reads
}

//
// TRIMMING SUBWORKFLOW
//
workflow TRIMMING {
    take:
    short_reads     // Channel: [[meta], [reads]]
    fastqc_zip      // Channel: [[meta], zip_files] from PRE_QC

    main:
    
    // Join short reads with their corresponding FastQC results
    ch_trimming_input = short_reads
        .join(fastqc_zip, by: 0)
        .map { meta, reads, zip -> 
            [meta, reads, zip]
        }
    
    // Run adaptive trimming with FASTP
    FASTP(ch_trimming_input)

    emit:
    // Trimmed reads for POST-QC
    trimmed_reads = FASTP.out.trimmed_reads
    
    // FASTP QC outputs
    fastp_html = FASTP.out.html
    fastp_json = FASTP.out.json
}

//
// POST-QC SUBWORKFLOW
//
workflow POST_QC {
    take:
    trimmed_short_reads  // Channel: [[meta], [trimmed_reads]]
    long_reads          // Channel: [[meta], [reads]] - original long reads

    main:
    
    // Run FastQC on trimmed short reads
    FASTQC_POST(trimmed_short_reads)
    
    // Run LongQC on original long reads (same as pre-QC for comparison)
    ch_long_filtered = long_reads
        .filter { _meta, reads -> reads != null && reads.toString() != 'null' }
    
    LONGQC_POST(ch_long_filtered)
    ch_longqc_post_results = LONGQC_POST.out
    
    // Collect Post-QC FastQC results for MultiQC
    ch_post_multiqc = FASTQC_POST.out.zip
        .map { _meta, files -> files }
        .flatten()
        .collect()
    
    // Run Post-QC MultiQC report
    MULTIQC_POST(ch_post_multiqc, "post")

    emit:
    // FastQC outputs
    fastqc_html = FASTQC_POST.out.html
    fastqc_zip = FASTQC_POST.out.zip
    
    // LongQC outputs
    longqc_html = ch_longqc_post_results.html
    longqc_plots = ch_longqc_post_results.plots
    
    // MultiQC outputs
    multiqc_html = MULTIQC_POST.out.html
    multiqc_data = MULTIQC_POST.out.data
}

//
// MAIN QC WORKFLOW
//
workflow QC {
    take:
    samplesheet

    main:
    
    // Read the samplesheet and create channels
    ch_samples = Channel
        .fromPath(samplesheet)
        .splitCsv(header: true)

    // Split into short and long reads channels with better organization
    ch_samples
        .branch { row -> 
            short_reads: row.fastq_1 && row.fastq_2
                return [[id: row.sample_id, read_type: 'short'],
                        [file(row.fastq_1), file(row.fastq_2)]] 
            
            long_reads: row.long_reads
                return [[id: row.sample_id, read_type: 'long'],
                        [file(row.long_reads)]]
        }
        .set { read_channels }

    // Execute PRE-QC subworkflow
    PRE_QC(
        read_channels.short_reads,
        read_channels.long_reads
    )
    
    // Execute TRIMMING subworkflow (short reads only)
    TRIMMING(
        PRE_QC.out.short_reads_passthrough,
        PRE_QC.out.fastqc_zip
    )
    
    // Execute POST-QC subworkflow
    POST_QC(
        TRIMMING.out.trimmed_reads,
        PRE_QC.out.long_reads_passthrough
    )

    emit:
    // PRE-QC outputs
    fastqc_pre_html = PRE_QC.out.fastqc_html
    longqc_pre_html = PRE_QC.out.longqc_plots
    fastqc_pre_zip = PRE_QC.out.fastqc_zip
    longqc_pre_zip = PRE_QC.out.longqc_html
    multiqc_pre_html = PRE_QC.out.multiqc_html
    multiqc_pre_data = PRE_QC.out.multiqc_data

    // Trimming outputs
    fastp_trimmed_reads = TRIMMING.out.trimmed_reads
    fastp_html = TRIMMING.out.fastp_html
    fastp_json = TRIMMING.out.fastp_json

    // POST-QC outputs
    fastqc_post_html = POST_QC.out.fastqc_html
    longqc_post_html = POST_QC.out.longqc_plots
    fastqc_post_zip = POST_QC.out.fastqc_zip
    longqc_post_zip = POST_QC.out.longqc_html
    multiqc_post_html = POST_QC.out.multiqc_html
    multiqc_post_data = POST_QC.out.multiqc_data
}