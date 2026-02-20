include { EPISEGMIX_DISTFIT_HISTONE_TRAIN  } from '../../modules/local/episegmix_distfit_histone_train'
include { EPISEGMIX_DISTFIT_HISTONE_ASSESS } from '../../modules/local/episegmix_distfit_histone_assess'

workflow DISTRIBUTION_FITTING {

    take:
        ch_input_raw  
        val_dists

    main:
    
        // Extract histone data
        ch_histone_data = ch_input_raw
            .map { it -> tuple(it[0], it[1]) }
            .filter { it[1] }

        ch_histone_inputs = ch_histone_data
            .flatMap { meta, h ->
                def header = h.withReader { it.readLine() }
                def marks_list = header.trim().split(/\s+/).toList().drop(3)
                marks_list.collectMany { mark ->
                    val_dists.collect { dist -> tuple(meta, h, mark, dist) }
                }
            }

        ch_histone_models = EPISEGMIX_DISTFIT_HISTONE_TRAIN(ch_histone_inputs)

        // Prepare assess input
        ch_h_assess_input = ch_histone_data
            .map { meta, h -> tuple(meta.id, tuple(meta, h)) }
            .cross(ch_histone_models.model.map { m, f -> tuple(m.id, f) }.groupTuple())
            .map { input, models -> tuple(input[1][0], input[1][1], models[1]) }
            
        ch_histone_assess = EPISEGMIX_DISTFIT_HISTONE_ASSESS(ch_h_assess_input)

}
