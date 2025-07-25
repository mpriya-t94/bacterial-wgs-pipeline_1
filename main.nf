// Import QC workflow
include { QC } from './workflows/qc.nf'

workflow {
    // Parameter check
    if (!params.samplesheet) {
        error "Please provide a samplesheet using the '--samplesheet' parameter."
    }

    // Call the QC workflow with the provided samplesheet
    QC(params.samplesheet)

    // Print messages
    QC.out.fastqc_pre_html.view { meta, _files -> 
    "FastQC completed for sample: ${meta.id}" }

    QC.out.longqc_pre_html.view { meta, _files -> 
    "LongQC completed for sample: ${meta.id}" }

    QC.out.multiqc_pre_html.view { file ->
    "Pre-QC MultiQC report generated: ${file}" }

    QC.out.fastp_trimmed_reads.view { meta, _files ->
    "FASTP adaptive trimming completed for sample: ${meta.id}" }

    QC.out.fastqc_post_html.view { meta, _files ->
    "Post-QC FastQC completed for sample: ${meta.id}" }

    QC.out.longqc_post_html.view { meta, _files ->
    "Post-QC LongQC completed for sample: ${meta.id}" }

    QC.out.multiqc_post_html.view { file ->
    "Post-QC MultiQC report generated: ${file}" }

}