include { GENERATE_METHYLATION_BINS }   from '../../modules/local/genrate_methylation_bins.nf'
workflow GENERATE_BINS {
    
    take:
    ch_methyl_data
    ch_bins

    main:

    ch_methyl_grouped = ch_methyl_data
        .map { meta, rep, mark, file, mod, pe, tab ->
            [ meta, mod, tab ]
        }
        .groupTuple(by: 0)

        
    GENERATE_METHYLATION_BINS(
        ch_methyl_grouped.map { meta, mods, tabs ->
            [ 
              meta, 
              mods, 
              tabs 
            ]
        },
        ch_bins
    )


    emit:
    split_counts = GENERATE_METHYLATION_BINS.out.merged_counts
}
