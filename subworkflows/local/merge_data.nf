include { MERGE_COUNTS }        from '../../modules/local/merge_counts'
include { SPLIT_MERGED_COUNTS } from '../../modules/local/split_merged_counts'

workflow MERGE_DATA {
    
    take:
    ch_histone_data
    ch_methyl_data
    ch_bins
    ref_file

    main:

    ch_methyl_grouped = ch_methyl_data
        .map { meta, rep, mark, file, mod, pe, tab ->
            [ meta, mod, tab ]
        }
        .groupTuple(by: 0)

    ch_merged_inputs = ch_histone_data
        .join(ch_methyl_grouped)
        
    MERGE_COUNTS(
        ch_merged_inputs.map { meta, counts, mods, tabs ->
            [ 
              meta, 
              [], [], [], 
              mods, 
              [], 
              counts, 
              tabs 
            ]
        },
        params.genome,
        ref_file,
        ch_bins
    )

    SPLIT_MERGED_COUNTS(
        MERGE_COUNTS.out.merged_counts
    )

    emit:
    split_counts = SPLIT_MERGED_COUNTS.out.inputs
}
