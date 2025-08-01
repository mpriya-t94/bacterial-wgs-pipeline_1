// Import assembly modules
include { SKESA } from '../modules/skesa/skesa.nf'
include { HYBRACTER } from '../modules/hybracter/hybracter.nf'

workflow ASSEMBLY {
    take:
    trimmed_short_reads  // Channel from QC workflow FASTP output
    raw_long_reads      // Channel from QC workflow (long reads passed through)
    samplesheet         // Still needed for metadata and determining read types
    
    main:
    
    // Read the samplesheet to get sample metadata and determine read types
    ch_samplesheet = Channel
        .fromPath(samplesheet)
        .splitCsv(header: true)
        .map { row -> 
            [row.sample_id, [
                has_short: row.fastq_1 && row.fastq_2,
                has_long: row.long_reads,
                read_type: (row.fastq_1 && row.fastq_2 && row.long_reads) ? 'hybrid' :
                          (row.fastq_1 && row.fastq_2) ? 'short_only' : 'long_only'
            ]]
        }
    
    // Create a map for easy lookup of sample read types
    ch_sample_info = ch_samplesheet
        .collectFile(name: 'sample_info.tmp') { sample_id, info ->
            "${sample_id}\t${info.read_type}\n"
        }
    
    // Process trimmed short reads and add read type information
    ch_processed_short = trimmed_short_reads
        .combine(ch_samplesheet)
        .filter { meta, reads, sample_id, info -> 
            meta.id == sample_id && info.has_short 
        }
        .map { meta, reads, sample_id, info ->
            [[id: meta.id, read_type: info.read_type], reads]
        }
    
    // Process raw long reads and add read type information  
    ch_processed_long = raw_long_reads
        .combine(ch_samplesheet)
        .filter { meta, reads, sample_id, info ->
            meta.id == sample_id && info.has_long
        }
        .map { meta, reads, sample_id, info ->
            [[id: meta.id, read_type: info.read_type], reads]
        }
    
    // Split processed reads by assembly strategy
    ch_processed_short
        .branch { meta, reads ->
            short_only: meta.read_type == 'short_only'
                return [meta, reads]
            hybrid: meta.read_type == 'hybrid'
                return [meta, reads]
        }
        .set { short_read_channels }
    
    ch_processed_long
        .branch { meta, reads ->
            long_only: meta.read_type == 'long_only'
                return [meta, reads]
            hybrid: meta.read_type == 'hybrid'
                return [meta, reads]
        }
        .set { long_read_channels }
    
    // Run SKESA assembly for short-only samples
    SKESA(short_read_channels.short_only)
    
    // Prepare channels for HYBRACTER
    
    // For long-only samples
    ch_hybracter_long_only = long_read_channels.long_only
        .map { meta, long_reads -> 
            [meta, long_reads.flatten()[0]]  // Ensure single file, not list
        }
    
    // For hybrid samples - combine short and long reads by sample ID
    ch_hybrid_combined = short_read_channels.hybrid
        .join(long_read_channels.hybrid, by: 0)  // Join by meta (sample ID)
        .map { meta, short_reads, long_reads ->
            [meta, short_reads, long_reads.flatten()[0]]  // meta, [short_R1, short_R2], long_reads_file
        }
    
    // Split hybrid channel for HYBRACTER inputs
    ch_hybracter_hybrid_long = ch_hybrid_combined
        .map { meta, short_reads, long_reads ->
            [meta, long_reads]  // Extract long reads
        }
    
    ch_hybracter_hybrid_short = ch_hybrid_combined
        .map { meta, short_reads, long_reads ->
            [meta, short_reads]  // Extract short reads
        }
    
    // Combine all long reads for HYBRACTER
    ch_hybracter_all_long = ch_hybracter_long_only
        .mix(ch_hybracter_hybrid_long)
    
    // Create short reads channel for HYBRACTER (empty list for long-only samples)
    ch_hybracter_all_short = ch_hybracter_hybrid_short
        .mix(
            ch_hybracter_long_only.map { meta, long_reads -> 
                [meta, []]  // Empty short reads for long-only samples
            }
        )
    
    // Run HYBRACTER assembly
    HYBRACTER(
        ch_hybracter_all_long,
        ch_hybracter_all_short
    )
    
    emit:
    // SKESA outputs
    skesa_contigs = SKESA.out.contigs
    
    // HYBRACTER outputs  
    hybracter_complete = HYBRACTER.out.complete_assemblies
    hybracter_incomplete = HYBRACTER.out.incomplete_assemblies
    hybracter_plasmids = HYBRACTER.out.plasmids
    hybracter_summary = HYBRACTER.out.summary
    hybracter_processing = HYBRACTER.out.processing_summary
    
    // Combined assemblies channel (all assembly outputs)
    all_assemblies = SKESA.out.contigs
        .mix(HYBRACTER.out.complete_assemblies)
        .mix(HYBRACTER.out.incomplete_assemblies)
}
