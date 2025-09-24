// Import assembly modules
include { SKESA } from '../modules/skesa/skesa.nf'
include { HYBRACTER } from '../modules/hybracter/hybracter.nf'

workflow ASSEMBLY {
    take:
    trimmed_short_reads  // Channel: [[meta], [reads]] from FASTP
    filtered_long_reads  // Channel: [[meta], [reads]] from FILTLONG
    
    main:
    
    // 1. Run SKESA on all short reads
    SKESA(trimmed_short_reads)
    
    // 2. For HYBRACTER - simpler approach using join operations
    
    // Try to join short and long reads by sample ID
    ch_joined = trimmed_short_reads
        .join(filtered_long_reads, by: 0, remainder: true)
        .map { meta, short_reads, long_reads ->
            [meta, short_reads ?: [], long_reads ?: []]
        }
        .filter { meta, short_reads, long_reads ->
            // Only process samples that have long reads
            long_reads.size() > 0
        }
    
    // Also handle samples that have only long reads (no short reads)
    ch_long_only = filtered_long_reads
        .join(trimmed_short_reads, by: 0, remainder: true)
        .filter { meta, long_reads, short_reads ->
            short_reads == null
        }
        .map { meta, long_reads, _short_reads ->
            [meta, [], long_reads]  // Empty short reads, keep long reads
        }
    
    // Combine both scenarios
    ch_all_for_hybracter = ch_joined.mix(ch_long_only)
    
    // Prepare HYBRACTER inputs
    ch_hybracter_long = ch_all_for_hybracter
        .map { meta, short_reads, long_reads -> [meta, long_reads] }
    
    ch_hybracter_short = ch_all_for_hybracter
        .map { meta, short_reads, long_reads -> [meta, short_reads] }
    
    // Run HYBRACTER
    HYBRACTER(
        ch_hybracter_long,
        ch_hybracter_short
    )
    
    emit:
    // SKESA outputs (short-read assemblies)
    skesa_contigs = SKESA.out.contigs
    
    // HYBRACTER outputs (hybrid and long-read assemblies)
    hybracter_complete = HYBRACTER.out.complete_assemblies
    hybracter_incomplete = HYBRACTER.out.incomplete_assemblies
    // hybracter_plasmids = HYBRACTER.out.plasmids
    hybracter_summary = HYBRACTER.out.summary
    // hybracter_processing = HYBRACTER.out.processing_summary
    
    // Combined assemblies channel
    all_assemblies = SKESA.out.contigs
        .mix(HYBRACTER.out.complete_assemblies)
        .mix(HYBRACTER.out.incomplete_assemblies)
}