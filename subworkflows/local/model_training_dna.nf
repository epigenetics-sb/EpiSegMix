include { EPISEGMIX_DNA_PREPARE } from '../../modules/local/episegmix_dna_prepare'
include { EPISEGMIX_DNA_TRAIN   } from '../../modules/local/episegmix_dna_train'
include { EPISEGMIX_DNA_DECODE  } from '../../modules/local/episegmix_dna_decode'
include { EPISEGMIX_DNA_REPORT  } from '../../modules/local/episegmix_dna_report'

workflow MODEL_TRAINING_DNA {

    take:
    ch_input

    main:

    ch_inputs_split = ch_input.flatMap { meta, meth, states ->

        def state_list = (states instanceof String)
            ? states.split(',')
            : states

        state_list.collect { state ->
            def new_meta = meta.clone()
            new_meta.id    = "${meta.id}_s${state}"
            new_meta.state = state

            [ new_meta, meth, state ]
        }
    }

    EPISEGMIX_DNA_PREPARE(ch_inputs_split)

    ch_train_inputs = EPISEGMIX_DNA_PREPARE.out.config
        .join(EPISEGMIX_DNA_PREPARE.out.train_counts)
        .join(EPISEGMIX_DNA_PREPARE.out.train_regions)

    EPISEGMIX_DNA_TRAIN(ch_train_inputs)

    ch_decode_inputs = EPISEGMIX_DNA_PREPARE.out.config
        .join(EPISEGMIX_DNA_TRAIN.out.model)
        .join(EPISEGMIX_DNA_PREPARE.out.counts)
        .join(EPISEGMIX_DNA_PREPARE.out.train_regions)
        
    EPISEGMIX_DNA_DECODE(ch_decode_inputs)

    ch_report_inputs = EPISEGMIX_DNA_PREPARE.out.config
        .join(EPISEGMIX_DNA_TRAIN.out.model)
        .join(EPISEGMIX_DNA_DECODE.out.bed)
        .join(EPISEGMIX_DNA_DECODE.out.seg_txt)

     EPISEGMIX_DNA_REPORT(ch_report_inputs)
}
