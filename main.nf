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
    QC.out.fastqc_html.view { meta, _files -> 
    "FastQC completed for sample: ${meta.id}" }

    QC.out.nanoplot_html.view { meta, _files -> 
    "NanoPlot completed for sample: ${meta.id}" }

}