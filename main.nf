#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { EPISEGMIX_PIPELINE } from './workflows/episegmix'

workflow {
    EPISEGMIX_PIPELINE ()
}
