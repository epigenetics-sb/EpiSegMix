include { CLEAN_AND_INDEX_BAM }     from '../../modules/local/clean_and_index_bam'
include { GENERATE_COUNT_MATRIX_BAM } from '../../modules/local/generate_count_matrix_bam'

workflow PROCESS_HISTONES {
    
    take:
    ch_input_bam
    chrom_sizes

    main:

    CLEAN_AND_INDEX_BAM(
        ch_input_bam
    )

    ch_grouped_bams = CLEAN_AND_INDEX_BAM.out.sanitized_bam
        .groupTuple(by: 0)
        .map { meta, bams, bais, reps, marks, mods, pes ->
            [ meta, bams, bais, reps, marks, mods, pes ]
        }

    GENERATE_COUNT_MATRIX_BAM(
        ch_grouped_bams,
        params.genome,
        chrom_sizes
    )

    emit:
    counts = GENERATE_COUNT_MATRIX_BAM.out.counts
}
