include { EPISEGMIX_STD_PREPARE } from '../../../modules/local/episegmix_std_prepare/main.nf'
include { EPISEGMIX_STD_TRAIN   } from '../../../modules/local/episegmix_std_train/main.nf'
include { EPISEGMIX_STD_DECODE  } from '../../../modules/local/episegmix_std_decode/main.nf'
include { EPISEGMIX_STD_REPORT  } from '../../../modules/local/episegmix_std_report/main.nf'

workflow MODEL_TRAINING_STD {

    take:
    ch_input

    main:
    ch_versions = Channel.empty()

    // Create unique entries for each requested state number
    ch_std_prepare_in = ch_input
        .flatMap { meta, histone, meth ->
            def state_list = (params.states instanceof String) ? params.states.split(',') : params.states
            state_list.collect { state ->
                def new_meta = meta.clone()
                new_meta.id = "${meta.id}_s${state}"
                [ new_meta, histone, meth, state ]
            }
        }

    // Generate training configuration and count files
    EPISEGMIX_STD_PREPARE (
        ch_std_prepare_in,    
        params.chr_parameter_estimation   
    )
    ch_versions = ch_versions.mix(EPISEGMIX_STD_PREPARE.out.versions)

    // Assemble components for the training phase
    ch_train_inputs = EPISEGMIX_STD_PREPARE.out.config
        .join(EPISEGMIX_STD_PREPARE.out.train_counts)
        .join(EPISEGMIX_STD_PREPARE.out.regions)
        .join(EPISEGMIX_STD_PREPARE.out.train_meth)

    // Execute standard HMM model training
    EPISEGMIX_STD_TRAIN(ch_train_inputs)
    ch_versions = ch_versions.mix(EPISEGMIX_STD_TRAIN.out.versions)

    // Match model and config for Viterbi decoding
    ch_decode_inputs = EPISEGMIX_STD_PREPARE.out.config
        .join(EPISEGMIX_STD_TRAIN.out.model)

    // Generate genome-wide state segmentations
    EPISEGMIX_STD_DECODE(ch_decode_inputs)
    ch_versions = ch_versions.mix(EPISEGMIX_STD_DECODE.out.versions)

    // Group all data for downstream visualization
    ch_report_inputs = EPISEGMIX_STD_PREPARE.out.config
        .join(EPISEGMIX_STD_TRAIN.out.model)
        .join(EPISEGMIX_STD_DECODE.out.bed)
        .join(EPISEGMIX_STD_DECODE.out.seg_txt)

    // Create final analysis report and plots
    EPISEGMIX_STD_REPORT(ch_report_inputs)
    ch_versions = ch_versions.mix(EPISEGMIX_STD_REPORT.out.versions)

    emit:
    versions = ch_versions
}