include { GENERATE_COUNT_MATRIX_BED } from '../../modules/local/generate_count_matrix_bed'

workflow PROCESS_METHYL {
    
    take:
    ch_input_bed

    main:

    GENERATE_COUNT_MATRIX_BED(
        ch_input_bed
    )

    emit:
    counts = GENERATE_COUNT_MATRIX_BED.out.counts_bed
}
