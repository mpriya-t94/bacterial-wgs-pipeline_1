nextflow.enable.dsl = 2

// Import modules (uncomment and adjust paths as needed)
// include { HYBRACTER } from '../modules/hybracter/hybracter.nf'
// include { SKESA } from '../modules/skesa/skesa.nf'

// If modules are in the same directory, include them directly
include { HYBRACTER } from '../modules/hybracter/hybracter.nf'
include { SKESA } from '../modules/skesa/skesa.nf'

//
// HYBRID/LONG-READ ASSEMBLY SUBWORKFLOW
//
workflow HYBRID_LONG_ASSEMBLY {
    take:
    long_reads      // Channel: [[meta], [long_reads]]
    short_reads     // Channel: [[meta], [short_reads]] - optional for hybrid
    assembly_mode   // val: 'hybrid', 'hybrid-single', 'long', 'long-single'
    chromosome_size // val: chromosome size or 'auto'

    main:
    
    // Prepare input channel for HYBRACTER
    // Combine long and short reads based on sample ID
    ch_hybracter_input = long_reads
        .map { meta, reads -> 
            def new_meta = meta.clone()
            new_meta.assembly_type = 'hybracter'
            return [new_meta, reads, [], assembly_mode, chromosome_size]
        }
    
    // For hybrid modes, join with short reads
    if (assembly_mode in ['hybrid', 'hybrid-single']) {
        ch_hybracter_input = long_reads
            .join(short_reads, by: 0, remainder: true)
            .map { meta, long_r, short_r -> 
                def new_meta = meta.clone()
                new_meta.assembly_type = 'hybracter'
                return [new_meta, long_r ?: [], short_r ?: [], assembly_mode, chromosome_size]
            }
    }
    
    // Run HYBRACTER
    HYBRACTER(ch_hybracter_input)

    emit:
    complete_assemblies = HYBRACTER.out.complete_assemblies
    incomplete_assemblies = HYBRACTER.out.incomplete_assemblies
    plasmids = HYBRACTER.out.plasmids
    summary = HYBRACTER.out.summary
    processing_summary = HYBRACTER.out.processing_summary
    final_output_dir = HYBRACTER.out.final_output_dir
    versions = HYBRACTER.out.versions
}

//
// SHORT-READ ASSEMBLY SUBWORKFLOW
//
workflow SHORT_READ_ASSEMBLY {
    take:
    short_reads  // Channel: [[meta], [reads]]

    main:
    
    // Prepare input for SKESA
    ch_skesa_input = short_reads
        .map { meta, reads -> 
            def new_meta = meta.clone()
            new_meta.assembly_type = 'skesa'
            return [new_meta, reads]
        }
    
    // Run SKESA
    SKESA(ch_skesa_input)

    emit:
    contigs = SKESA.out.contigs
}

//
// MAIN ASSEMBLY WORKFLOW
//
workflow ASSEMBLY {
    take:
    samplesheet         // Path to samplesheet
    assembly_modes      // Map: [sample_id: assembly_mode] - optional
    chromosome_sizes    // Map: [sample_id: chromosome_size] - optional

    main:
    
    // Read the samplesheet and create channels
    ch_samples = Channel
        .fromPath(samplesheet)
        .splitCsv(header: true)

    // Branch samples based on available read types and requested assembly modes
    ch_samples
        .branch { row -> 
            // Determine assembly mode for this sample
            def sample_mode = assembly_modes?.get(row.sample_id) ?: 'auto'
            def chr_size = chromosome_sizes?.get(row.sample_id) ?: 'auto'
            
            // Classify based on available reads and mode
            hybrid_assembly: (row.long_reads && row.fastq_1 && row.fastq_2) && 
                           (sample_mode in ['hybrid', 'hybrid-single', 'auto'])
                return [
                    [id: row.sample_id, read_type: 'hybrid'], 
                    [file(row.long_reads)],
                    [file(row.fastq_1), file(row.fastq_2)],
                    sample_mode == 'auto' ? 'hybrid-single' : sample_mode,
                    chr_size
                ]
            
            long_only_assembly: row.long_reads && 
                              (!row.fastq_1 || !row.fastq_2 || sample_mode in ['long', 'long-single'])
                return [
                    [id: row.sample_id, read_type: 'long'], 
                    [file(row.long_reads)],
                    [],
                    sample_mode == 'auto' ? 'long-single' : sample_mode,
                    chr_size
                ]
            
            short_only_assembly: (row.fastq_1 && row.fastq_2) && 
                               (!row.long_reads || sample_mode == 'short')
                return [
                    [id: row.sample_id, read_type: 'short'],
                    [file(row.fastq_1), file(row.fastq_2)]
                ]
        }
        .set { assembly_channels }

    // Run hybrid/long-read assemblies
    ch_hybrid_long_input = assembly_channels.hybrid_assembly
        .mix(assembly_channels.long_only_assembly)
    
    if (!ch_hybrid_long_input.empty) {
        // Split the channel for HYBRID_LONG_ASSEMBLY inputs
        ch_hybrid_long_input
            .map { meta, long_reads, short_reads, mode, chr_size ->
                [
                    [meta, long_reads], 
                    [meta, short_reads], 
                    mode, 
                    chr_size
                ]
            }
            .transpose()
            .set { ch_split_inputs }
        
        HYBRID_LONG_ASSEMBLY(
            ch_split_inputs[0], // long_reads
            ch_split_inputs[1], // short_reads  
            ch_split_inputs[2], // assembly_mode
            ch_split_inputs[3]  // chromosome_size
        )
    }

    // Run short-read assemblies
    if (!assembly_channels.short_only_assembly.empty) {
        SHORT_READ_ASSEMBLY(assembly_channels.short_only_assembly)
    }

    // Collect versions (only from HYBRACTER)
    ch_versions = Channel.empty()
    if (!ch_hybrid_long_input.empty) {
        ch_versions = ch_versions.mix(HYBRID_LONG_ASSEMBLY.out.versions)
    }

    emit:
    // HYBRACTER/Long-read outputs (if any)
    hybracter_complete = !ch_hybrid_long_input.empty ? HYBRID_LONG_ASSEMBLY.out.complete_assemblies : Channel.empty()
    hybracter_incomplete = !ch_hybrid_long_input.empty ? HYBRID_LONG_ASSEMBLY.out.incomplete_assemblies : Channel.empty()
    hybracter_plasmids = !ch_hybrid_long_input.empty ? HYBRID_LONG_ASSEMBLY.out.plasmids : Channel.empty()
    hybracter_summary = !ch_hybrid_long_input.empty ? HYBRID_LONG_ASSEMBLY.out.summary : Channel.empty()
    hybracter_processing_summary = !ch_hybrid_long_input.empty ? HYBRID_LONG_ASSEMBLY.out.processing_summary : Channel.empty()
    hybracter_final_output = !ch_hybrid_long_input.empty ? HYBRID_LONG_ASSEMBLY.out.final_output_dir : Channel.empty()
    
    // SKESA outputs (if any)
    skesa_contigs = !assembly_channels.short_only_assembly.empty ? SHORT_READ_ASSEMBLY.out.contigs : Channel.empty()
    
    // Combined outputs
    versions = ch_versions
}

//
// INTEGRATION WITH QC WORKFLOW
//
workflow QC_ASSEMBLY {
    take:
    samplesheet
    assembly_modes      // optional
    chromosome_sizes    // optional

    main:
    
    // Import QC workflow (adjust path as needed)
    // include { QC } from './qc.nf'
    
    // Run QC first
    // QC(samplesheet)
    
    // For now, run assembly directly on original reads
    // In a full pipeline, you'd use QC.out.fastp_trimmed_reads for short reads
    ASSEMBLY(samplesheet, assembly_modes, chromosome_sizes)

    emit:
    // Assembly outputs
    hybracter_complete = ASSEMBLY.out.hybracter_complete
    hybracter_incomplete = ASSEMBLY.out.hybracter_incomplete
    hybracter_plasmids = ASSEMBLY.out.hybracter_plasmids
    hybracter_summary = ASSEMBLY.out.hybracter_summary
    hybracter_processing_summary = ASSEMBLY.out.hybracter_processing_summary
    hybracter_final_output = ASSEMBLY.out.hybracter_final_output
    skesa_contigs = ASSEMBLY.out.skesa_contigs
    versions = ASSEMBLY.out.versions
}

//
// DEFAULT WORKFLOW FOR TESTING
//
workflow {
    // Example parameters
    def assembly_modes = [
        'sample1': 'hybrid-single',
        'sample2': 'long-single', 
        'sample3': 'short'
    ]
    
    def chromosome_sizes = [
        'sample1': 'auto',
        'sample2': '5000000',
        'sample3': 'auto'
    ]

    // Example samplesheet path
    def samplesheet = params.input ?: 'samplesheet.csv'
    
    // Run the assembly workflow
    ASSEMBLY(samplesheet, assembly_modes, chromosome_sizes)

    // View outputs
    ASSEMBLY.out.hybracter_complete.view { "HYBRACTER complete: $it" }
    ASSEMBLY.out.skesa_contigs.view { "SKESA contigs: $it" }
    ASSEMBLY.out.versions.collectFile(name: 'assembly_versions.yml')
}