include { EPISEGMIX_DM_PREPARE } from '../../modules/local/episegmix_dm_prepare'
include { EPISEGMIX_DM_TRAIN   } from '../../modules/local/episegmix_dm_train'
include { EPISEGMIX_DM_DECODE  } from '../../modules/local/episegmix_dm_decode'
include { EPISEGMIX_DM_REPORT  } from '../../modules/local/episegmix_dm_report'

workflow MODEL_TRAINING_DM {

    take:
    ch_input

    main:

    ch_inputs_split = ch_input.flatMap { meta, histone, meth, states ->

        def state_list = (states instanceof String)
            ? states.split(',')
            : states

        state_list.collect { state ->
            def new_meta = meta.clone()
            new_meta.id    = "${meta.id}_s${state}"
            new_meta.state = state

            [ new_meta, histone, meth, state ]
        }
    }

    EPISEGMIX_DM_PREPARE(ch_inputs_split, params.chr_parameter_estimation)

    ch_train_inputs = EPISEGMIX_DM_PREPARE.out.config
        .join(EPISEGMIX_DM_PREPARE.out.train_counts)
        .join(EPISEGMIX_DM_PREPARE.out.regions)
        .join(EPISEGMIX_DM_PREPARE.out.train_meth)

    EPISEGMIX_DM_TRAIN(ch_train_inputs)

    ch_decode_inputs = EPISEGMIX_DM_PREPARE.out.config
        .join(EPISEGMIX_DM_TRAIN.out.model)

    EPISEGMIX_DM_DECODE(ch_decode_inputs)

    ch_report_inputs = EPISEGMIX_DM_PREPARE.out.config
        .join(EPISEGMIX_DM_TRAIN.out.model)
        .join(EPISEGMIX_DM_DECODE.out.bed)
        .join(EPISEGMIX_DM_DECODE.out.seg_txt)

    EPISEGMIX_DM_REPORT(ch_report_inputs)
}
