include { EPISEGMIX_STD_PREPARE } from '../../modules/local/episegmix_std_prepare'
include { EPISEGMIX_STD_TRAIN   } from '../../modules/local/episegmix_std_train'
include { EPISEGMIX_STD_DECODE  } from '../../modules/local/episegmix_std_decode'
include { EPISEGMIX_STD_REPORT  } from '../../modules/local/episegmix_std_report'

workflow MODEL_TRAINING_STD {

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

    EPISEGMIX_STD_PREPARE(ch_inputs_split, params.chr_parameter_estimation)

    ch_train_inputs = EPISEGMIX_STD_PREPARE.out.config
        .join(EPISEGMIX_STD_PREPARE.out.train_counts)
        .join(EPISEGMIX_STD_PREPARE.out.regions)
        .join(EPISEGMIX_STD_PREPARE.out.train_meth)

    EPISEGMIX_STD_TRAIN(ch_train_inputs)

    ch_decode_inputs = EPISEGMIX_STD_PREPARE.out.config
        .join(EPISEGMIX_STD_TRAIN.out.model)

    EPISEGMIX_STD_DECODE(ch_decode_inputs)

    ch_report_inputs = EPISEGMIX_STD_PREPARE.out.config
        .join(EPISEGMIX_STD_TRAIN.out.model)
        .join(EPISEGMIX_STD_DECODE.out.bed)
        .join(EPISEGMIX_STD_DECODE.out.seg_txt)

    EPISEGMIX_STD_REPORT(ch_report_inputs)
}
