nextflow.enable.dsl = 2

include { HYBRACTER } from './hybracter.nf'

workflow {

    // hybrid mode: multiple samples from sample sheet
    def hybrid_meta = [id: 'hybrid']
    def hybrid_input = tuple(
        hybrid_meta,
        file("/home/mpriya_t94/hybracter/hybrid.txt"),
        [],
        'hybrid',
        'auto'
    )

    // hybrid-single mode: one sample, direct files
    def hybrid_single_meta = [id: 'hybrid-single']
    def hybrid_single_input = tuple(
    hybrid_single_meta,
    file('11-1036.fastq.gz'),
    [file('11-1036_R1.fastq.gz'), file('11-1036_R2.fastq.gz')],
    'hybrid-single',
    'auto'
    )


    // long mode: multiple samples from sample sheet
    def long_meta = [id: 'long']
    def long_input = tuple(
        long_meta,
        file("/home/mpriya_t94/hybracter/long.txt"),
        [],
        'long',
        'auto'
    )

    // long-single mode: one sample, direct file
    def long_single_meta = [id: 'long-single']
    def long_single_input = tuple(
        long_single_meta,
        file('11-1036.fastq.gz'),
        [],
        'long-single',
        'auto'
    )

    Channel
        .of(hybrid_input, hybrid_single_input, long_input, long_single_input)
        .set { hybracter_input }

    HYBRACTER(hybracter_input)
    HYBRACTER.out.complete_assemblies.view()
}