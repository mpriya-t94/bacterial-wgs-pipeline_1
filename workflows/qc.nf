// Import FASTQC and NANOPLOT modules
include { FASTQC as FASTQC_PRE } from '../modules/fastqc/fastqc.nf'
include { FASTQC as FASTQC_POST } from '../modules/fastqc/fastqc.nf'
include {LONGQC as LONGQC_PRE } from '../modules/longqc/longqc.nf'
include { LONGQC as LONGQC_POST } from '../modules/longqc/longqc.nf'
include { FASTP } from '../modules/fastp/fastp.nf'
include {MULTIQC as MULTIQC_PRE} from '../modules/multiqc/multiqc.nf'
include { MULTIQC as MULTIQC_POST } from '../modules/multiqc/multiqc.nf'

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
    .set { read_channels }

    // Run PRE-QC for short and long reads
    FASTQC_PRE(read_channels.short_reads)
    LONGQC_PRE(read_channels.long_reads)

    // Run PRE-QC MultiQC for short reads only
    ch_pre_multiqc = FASTQC_PRE.out.zip
    .map { _meta, files -> files }
    .flatten()
    .collect()

    MULTIQC_PRE(ch_pre_multiqc, "pre")

    // Run FASTP adaptive trimming for short reads
    ch_trimming = read_channels.short_reads
    .join(FASTQC_PRE.out.zip, by: 0)

    FASTP(ch_trimming)

    // Run POST-QC for short and long reads
    FASTQC_POST(FASTP.out.trimmed_reads)
    LONGQC_POST(read_channels.long_reads)

    // Run POST-QC MultiQC for short reads only
    ch_post_multiqc = FASTQC_POST.out.zip
    .map { _meta, files -> files }
    .flatten()
    .collect()

    MULTIQC_POST(ch_post_multiqc, "post")

    emit:
    // PRE-QC outputs
    fastqc_pre_html = FASTQC_PRE.out.html
    longqc_pre_html = LONGQC_PRE.out.plots
    fastqc_pre_zip = FASTQC_PRE.out.zip
    longqc_pre_zip = LONGQC_PRE.out.html
    multiqc_pre_html = MULTIQC_PRE.out.html
    multiqc_pre_data = MULTIQC_PRE.out.data

    // Trimming outputs
    fastp_trimmed_reads = FASTP.out.trimmed_reads
    fastp_html = FASTP.out.html
    fastp_json = FASTP.out.json

    // POST-QC outputs
    fastqc_post_html = FASTQC_POST.out.html
    longqc_post_html = LONGQC_POST.out.plots
    fastqc_post_zip = FASTQC_POST.out.zip
    longqc_post_zip = LONGQC_POST.out.html
    multiqc_post_html = MULTIQC_POST.out.html
    multiqc_post_data = MULTIQC_POST.out.data
}