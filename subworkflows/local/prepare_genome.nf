include { GET_GENOME_CHROM_SIZE } from '../../modules/local/get_genome_chrom_size'
include { GENERATE_GENOME_BINS } from '../../modules/local/generate_genome_bins'

workflow PREPARE_GENOME {
    
    main:
    GET_GENOME_CHROM_SIZE(
        params.genome
    )

    GENERATE_GENOME_BINS(
        GET_GENOME_CHROM_SIZE.out.genome_chrom_size_file,
        params.genome,
        params.binsize
    )

    emit:
    chrom_sizes = GET_GENOME_CHROM_SIZE.out.genome_chrom_size_file
    bins        = GENERATE_GENOME_BINS.out.bins_file
}
