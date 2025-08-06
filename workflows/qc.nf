// Import modules
include { FASTQC as FASTQC_PRE } from '../modules/fastqc/fastqc.nf'
include { FASTQC as FASTQC_POST } from '../modules/fastqc/fastqc.nf'
include { NANOPLOT as NANOPLOT_PRE } from '../modules/nanoplot/nanoplot.nf'
include { NANOPLOT as NANOPLOT_POST } from '../modules/nanoplot/nanoplot.nf'
include { FASTP } from '../modules/fastp/fastp.nf'
include { FILTLONG } from '../modules/filtlong/filtlong.nf'
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
    
    // Run nanoplot on long reads in parallel (only if long reads exist)
    ch_long_filtered = long_reads
        .filter { _meta, reads -> reads != null && reads.toString() != 'null' }
    
    NANOPLOT_PRE(ch_long_filtered)
    ch_nanoplot_results = NANOPLOT_PRE.out
    
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
    
    // nanoplot outputs
    nanoplot_html = ch_nanoplot_results.html
    nanoplot_plots = ch_nanoplot_results.plots
    
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
    long_reads      // Channel: [[meta], [reads]]
    fastqc_zip      // Channel: [[meta], zip_files] from PRE_QC
    illumina_refs   // Channel: [[meta], [illumina_R1, illumina_R2]]
    assembly        // Channel: [[meta], assembly_file]

    main:
    
    // Join short reads with their corresponding FastQC results
    ch_trimming_input = short_reads
        .join(fastqc_zip, by: 0)
        .map { meta, reads, zip -> 
            [meta, reads, zip]
        }
    
    // Run adaptive trimming with FASTP
    FASTP(ch_trimming_input)

    // Long read processing through FILTLONG
    ch_long_for_filtlong = long_reads
        .filter { _meta, reads -> reads != null && reads.toString() != 'null' }

    // Prepare FILTLONG input
    // Change from combine to join:
    ch_filtlong_input = ch_long_for_filtlong
        .join(illumina_refs, by: 0, remainder: true)
        .join(assembly, by: 0, remainder: true)
        .map { meta, reads, illum, asm -> 
            [meta, reads, illum ?: [meta, []], asm ?: null]
    }

    // Run FILTLONG on long reads
    FILTLONG(
        ch_filtlong_input.map { meta, reads, illum, asm -> [meta, reads] },
        ch_filtlong_input.map { meta, reads, illum, asm -> illum },
        ch_filtlong_input.map { meta, reads, illum, asm -> asm }
    )

    emit:
    // Trimmed reads for POST-QC
    trimmed_reads = FASTP.out.trimmed_reads
    
    // FASTP QC outputs
    fastp_html = FASTP.out.html
    fastp_json = FASTP.out.json

    // FILTLONG outputs
    filtered_long_reads = FILTLONG.out.filtered_reads
    filtlong_log        = FILTLONG.out.log
}


//
// POST-QC SUBWORKFLOW
//
workflow POST_QC {
    take:
    trimmed_short_reads  // Channel: [[meta], [trimmed_reads]]
    filtered_long_reads  // Channel: [[meta], [reads]] - original long reads

    main:
    
    // Run FastQC on trimmed short reads
    FASTQC_POST(trimmed_short_reads)
    
    // Run nanoplot on FILTERED long reads
    ch_long_filtered = filtered_long_reads
        .filter { _meta, reads -> reads != null && reads.toString() != 'null' }
    
    NANOPLOT_POST(ch_long_filtered)
    ch_nanoplot_post_results = NANOPLOT_POST.out
    
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
    
    // nanoplot outputs
    nanoplot_html = ch_nanoplot_post_results.html
    nanoplot_plots = ch_nanoplot_post_results.plots
    
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
    ch_samples_split = ch_samples
    .multiMap { row -> 
        short_reads: 
            (row.fastq_1 && row.fastq_2) ? 
            [[id: row.sample_id, read_type: 'short'], [file(row.fastq_1), file(row.fastq_2)]] : null
        
        long_reads: 
            row.long_reads ? 
            [[id: row.sample_id, read_type: 'long'], [file(row.long_reads)]] : null
        
        illumina_refs: 
            (row.illumina_ref_1 && row.illumina_ref_2) ? 
            [[id: row.sample_id, read_type: 'illumina'], [file(row.illumina_ref_1), file(row.illumina_ref_2)]] : null
        
        assembly: 
            row.assembly ? 
            [[id: row.sample_id, read_type: 'assembly'], [file(row.assembly)]] : null
    }

    // Filter out nulls - this is the correct syntax
    ch_short_reads = ch_samples_split.short_reads.filter { it != null }
    ch_long_reads = ch_samples_split.long_reads.filter { it != null }
    ch_illumina_refs = ch_samples_split.illumina_refs.filter { it != null }
    ch_assembly = ch_samples_split.assembly.filter { it != null }

    // Execute PRE-QC subworkflow
    PRE_QC(
        ch_short_reads,
        ch_long_reads
    )
    
    // Execute TRIMMING subworkflow (short reads only)
    TRIMMING(
        PRE_QC.out.short_reads_passthrough,
        PRE_QC.out.long_reads_passthrough,
        PRE_QC.out.fastqc_zip,
        ch_illumina_refs,
        ch_assembly
    )
    
    // Execute POST-QC subworkflow
    POST_QC(
        TRIMMING.out.trimmed_reads,
        TRIMMING.out.filtered_long_reads
    )

    emit:
    // PRE-QC outputs
    fastqc_pre_html = PRE_QC.out.fastqc_html
    nanoplot_pre_html = PRE_QC.out.nanoplot_html
    fastqc_pre_zip = PRE_QC.out.fastqc_zip
    nanoplot_pre_plots = PRE_QC.out.nanoplot_plots
    multiqc_pre_html = PRE_QC.out.multiqc_html
    multiqc_pre_data = PRE_QC.out.multiqc_data

    // Trimming outputs
    fastp_trimmed_reads = TRIMMING.out.trimmed_reads
    fastp_html = TRIMMING.out.fastp_html
    fastp_json = TRIMMING.out.fastp_json
    filtlong_filtered_reads = TRIMMING.out.filtered_long_reads
    filtlong_log = TRIMMING.out.filtlong_log

    // POST-QC outputs
    fastqc_post_html = POST_QC.out.fastqc_html
    nanoplot_post_html = POST_QC.out.nanoplot_html
    fastqc_post_zip = POST_QC.out.fastqc_zip
    nanoplot_post_plots = POST_QC.out.nanoplot_plots
    multiqc_post_html = POST_QC.out.multiqc_html
    multiqc_post_data = POST_QC.out.multiqc_data
}