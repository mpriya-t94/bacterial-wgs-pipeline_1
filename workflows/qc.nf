// Import FASTQC and NANOPLOT modules
include { FASTQC } from '../modules/fastqc/fastqc.nf'
include { LONGQC } from '../modules/longqc/longqc.nf'

workflow QC {
    take:
    samplesheet

    main:

    // Read the samplesheet and create channels
    ch_samples = Channel
    .fromPath(samplesheet)
    .splitCsv(header: true)

    // Split into short and long reads channels
    ch_samples
    .branch { row -> short_reads: row.fastq_1 && row.fastq_2
    return [[id: row.sample_id, read_type: 'short'],
    [file(row.fastq_1), file(row.fastq_2)]] 
    
    long_reads: row.long_reads
    return [[id: row.sample_id, read_type: 'long'],
    [file(row.long_reads)]]

    }
    .set { reads_channels }

    // Run FASTQC for short reads
    FASTQC(reads_channels.short_reads)

    // Run NANOPLOT for long reads
    LONGQC(reads_channels.long_reads)

    emit:
    // Emit QC results for downstream use
    fastqc_html = FASTQC.out.html
    fastqc_zip = FASTQC.out.zip
    nanoplot_html = LONGQC.out.html
    nanoplot_plots = LONGQC.out.plots
}