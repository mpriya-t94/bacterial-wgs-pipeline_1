// Import workflows and modules
include { QC } from './workflows/qc.nf'
include { ASSEMBLY } from './workflows/assembly.nf'
include { QUAST } from './modules/quast/quast.nf'

workflow {
    // Parameter checks
    if (!params.samplesheet) {
        error "Please provide a samplesheet using the '--samplesheet' parameter."
    }

    // Call the QC workflow with the provided samplesheet
    QC(params.samplesheet)

    // Print QC completion messages
    QC.out.fastqc_pre_html.view { meta, _files -> 
        "FastQC completed for sample: ${meta.id}" 
    }

    QC.out.nanoplot_pre_html.view { meta, _files -> 
        "Nanoplot completed for sample: ${meta.id}" 
    }

    QC.out.multiqc_pre_html.view { file ->
        "Pre-QC MultiQC report generated: ${file}" 
    }

    QC.out.fastp_trimmed_reads.view { meta, _files ->
        "FASTP adaptive trimming completed for sample: ${meta.id}" 
    }

    QC.out.fastqc_post_html.view { meta, _files ->
        "Post-QC FastQC completed for sample: ${meta.id}" 
    }

    QC.out.nanoplot_post_html.view { meta, _files ->
        "Post-QC Nanoplot completed for sample: ${meta.id}" 
    }

    QC.out.multiqc_post_html.view { file ->
        "Post-QC MultiQC report generated: ${file}" 
    }

    // Call the ASSEMBLY workflow
    ASSEMBLY(
        QC.out.fastp_trimmed_reads,
        QC.out.filtlong_filtered_reads
    )

    // Print assembly completion messages
    ASSEMBLY.out.skesa_contigs.view { meta, _files ->
        "SKESA assembly completed for sample: ${meta.id}"
    }

    ASSEMBLY.out.hybracter_complete.view { meta, _files ->
        "HYBRACTER complete assembly generated for sample: ${meta.id}"
    }

    ASSEMBLY.out.hybracter_incomplete.view { meta, _files ->
        "HYBRACTER incomplete assembly generated for sample: ${meta.id}"
    }

    // Collect all assemblies for QUAST
    ch_all_assemblies = ASSEMBLY.out.skesa_contigs
        .mix(ASSEMBLY.out.hybracter_complete)
        .mix(ASSEMBLY.out.hybracter_incomplete)

    // Prepare optional reference files
    ch_reference = params.reference_genome ? 
        Channel.value(file(params.reference_genome)) : 
        Channel.value(file('NO_FILE'))
    
    ch_genes = params.gene_annotations ? 
        Channel.value(file(params.gene_annotations)) : 
        Channel.value(file('NO_FILE'))

    // Run QUAST assessment directly
    QUAST(
        ch_all_assemblies,
        ch_reference,
        ch_genes,
        QC.out.fastp_trimmed_reads
    )

    // Print QUAST completion messages
    QUAST.out.html.view { meta, _files ->
        "QUAST assessment completed for sample: ${meta.id}"
    }

    QUAST.out.tsv.view { meta, _files ->
        "QUAST TSV report generated for sample: ${meta.id}"
    }
}