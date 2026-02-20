#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { INPUT_CHECK } from '../subworkflows/local/input_check'
include { PREPARE_GENOME } from '../subworkflows/local/prepare_genome'
include { PROCESS_HISTONES } from '../subworkflows/local/process_histones'
include { PROCESS_METHYL } from '../subworkflows/local/process_methyl'
include { MERGE_DATA } from '../subworkflows/local/merge_data'
include { GENERATE_BINS } from '../subworkflows/local/generate_bins'

include { MODEL_TRAINING_STD } from '../subworkflows/local/model_training_std'
include { MODEL_TRAINING_DM } from '../subworkflows/local/model_training_dm'
include { MODEL_TRAINING_DNA } from '../subworkflows/local/model_training_dna'
include { DISTRIBUTION_FITTING } from '../subworkflows/local/distribution_fitting'

workflow EPISEGMIX_PIPELINE {

    INPUT_CHECK()
    PREPARE_GENOME()

    if (params.episegmix_mode != 'DNA') {
        PROCESS_HISTONES(INPUT_CHECK.out.bam, PREPARE_GENOME.out.chrom_sizes)
    }

    if (params.episegmix_mode == 'DNA' || params.merge) {
        PROCESS_METHYL(INPUT_CHECK.out.bed)
    }

    ch_inputs_raw = Channel.empty()

    if (params.episegmix_mode == 'DNA') {
        GENERATE_BINS(PROCESS_METHYL.out.counts, PREPARE_GENOME.out.bins)
        ch_inputs_raw = GENERATE_BINS.out.split_counts
    } else if (params.merge) {
        MERGE_DATA(PROCESS_HISTONES.out.counts, PROCESS_METHYL.out.counts, PREPARE_GENOME.out.bins, PREPARE_GENOME.out.chrom_sizes)
        ch_inputs_raw = MERGE_DATA.out.split_counts
    } else {
        ch_inputs_raw = PROCESS_HISTONES.out.counts.map { meta, counts -> [ meta, counts, [] ] }
    }

    if (params.episegmix_mode == 'fitting') {
        ch_fitting_input = ch_inputs_raw.map { tuple ->
            def meta = tuple[0]
            def h = tuple.size() > 1 ? tuple[1] : []
            return [ meta, h ]
        }

        DISTRIBUTION_FITTING(ch_fitting_input, params.distributions)
    } else {
        ch_model_input = Channel.empty()

        if (params.episegmix_mode == 'DNA') {
            ch_model_input = ch_inputs_raw.map { tuple ->
                def meta = tuple[0]
                def meth = tuple.size() > 1 ? tuple[1] : []
                def states = tuple.size() > 2 ? tuple[2] : params.states
                return [ meta, meth, states ]
            }
        } else {
            ch_model_input = ch_inputs_raw.combine(Channel.value(params.states))
        }

        if (params.episegmix_mode == 'duration') {
            MODEL_TRAINING_DM(ch_model_input)
        } else if (params.episegmix_mode == 'DNA') {
            MODEL_TRAINING_DNA(ch_model_input)
        } else {
            MODEL_TRAINING_STD(ch_model_input)
        }
    }
}