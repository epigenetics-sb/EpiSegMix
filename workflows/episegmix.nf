/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_episegmix_pipeline'

/* LOCAL SUBWORKFLOWS */
include { PREPARE_GENOME       } from '../subworkflows/local/prepare_genome/main.nf'
include { PROCESS_METHYL       } from '../subworkflows/local/process_methyl/main.nf'
include { PROCESS_HISTONES     } from '../subworkflows/local/process_histones/main.nf'
include { MERGE_DATA           } from '../subworkflows/local/merge_data/main.nf'
include { GENERATE_BINS        } from '../subworkflows/local/generate_bins/main.nf'
include { MODEL_TRAINING_STD   } from '../subworkflows/local/model_training_std/main.nf'
include { MODEL_TRAINING_DM    } from '../subworkflows/local/model_training_dm/main.nf'
include { MODEL_TRAINING_DNA   } from '../subworkflows/local/model_training_dna/main.nf'
include { DISTRIBUTION_FITTING } from '../subworkflows/local/distribution_fitting/main.nf' 

/* MODULES */
include { MULTIQC              } from '../modules/nf-core/multiqc/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow EPISEGMIX {

    take:
    ch_samplesheet
    ch_versions_init 

    main:
    ch_versions      = ch_versions_init
    ch_multiqc_files = Channel.empty()

    // ---------------------------------------------------------
    // DETERMINE PIPELINE MODE
    // ---------------------------------------------------------
    // Prioritize specific flags over the legacy 'episegmix_mode' string
    def active_mode = params.episegmix_mode
    if (params.fitting)       { active_mode = 'fitting'  }
    else if (params.duration) { active_mode = 'duration' }
    else if (params.dna)      { active_mode = 'DNA'      }
    else if (params.standard) { active_mode = 'standard' }

    // Create a dictionary of distributions for each sample and mark
    ch_distributions = ch_samplesheet
        .map { meta, file -> 
            def dist = meta.distribution ?: meta.distributions
            [ meta.id, meta.epigenetic_mark, dist ] 
        }
        .groupTuple(by: 0)
        .map { id, marks, dists ->
            def dist_map = [marks, dists].transpose().collectEntries { k, v -> [k, v] }
            return [ id, dist_map ]
        }

    // Step 1: Generate genomic bins and chromosome sizes
    PREPARE_GENOME()
    ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)

    // Step 2: Split input files into BAM (Histones) and BED (Methylation)
    ch_samplesheet
        .transpose() 
        .branch {
            meta, file ->
                def fname = file.name.toString().toLowerCase()
                bam: fname.endsWith('.bam') || fname.endsWith('.bam.gz')
                bed: fname.endsWith('.bed') || fname.endsWith('.bed.gz')
                other: true
        }.set { ch_input_branched }

    // Step 3: Process Histone BAM files
    if (active_mode != 'DNA') {
        PROCESS_HISTONES(ch_input_branched.bam, PREPARE_GENOME.out.chrom_sizes)
        ch_versions = ch_versions.mix(PROCESS_HISTONES.out.versions)
    }
    
    // Step 4: Process Methylation BED files
    if (active_mode == 'DNA' || params.merge) {
        PROCESS_METHYL(ch_input_branched.bed)
        ch_versions = ch_versions.mix(PROCESS_METHYL.out.versions)
    }

    // Step 5: Formatting and Merging Logic
    ch_model_input = Channel.empty()

    if (active_mode == 'DNA') {
        GENERATE_BINS(PROCESS_METHYL.out.counts, PREPARE_GENOME.out.bins)
        ch_model_input = GENERATE_BINS.out.split_counts
        ch_versions    = ch_versions.mix(GENERATE_BINS.out.versions)

    } else if (params.merge) {
        MERGE_DATA(PROCESS_HISTONES.out.counts, PROCESS_METHYL.out.counts, PREPARE_GENOME.out.bins, PREPARE_GENOME.out.chrom_sizes)
        ch_model_input = MERGE_DATA.out.split_counts
        ch_versions    = ch_versions.mix(MERGE_DATA.out.versions)

    } else {
        ch_model_input = PROCESS_HISTONES.out.counts.map { meta, counts_file -> [ meta, counts_file, [] ] }
    }

    // Step 6: Metadata Injection
    ch_ready_for_analysis = ch_model_input
        .map { tuple -> [ tuple[0].id, tuple ] } 
        .join(ch_distributions)                  
        .map { id, original_tuple, dist_map ->
            def new_meta = original_tuple[0].clone()
            new_meta.distributions = dist_map    
            return [ new_meta ] + original_tuple.drop(1) 
        }

    // Step 7: Route to Analysis Models
    if (active_mode == 'fitting') {
        def dist_list = params.distributions ? params.distributions.split(',').collect{ it.trim() } : []
        DISTRIBUTION_FITTING(ch_ready_for_analysis.map { [it[0], it[1]] }, dist_list, file(params.input))
        ch_versions = ch_versions.mix(DISTRIBUTION_FITTING.out.versions)

        if (params.best_fit_segmentation) {
            ch_updated_dist = DISTRIBUTION_FITTING.out.updated_samplesheet
                .splitCsv(header: true)
                .map { meta, row -> 
                    [ meta.id, row.epigenetic_mark, row.distribution ] 
                }
                .groupTuple(by: 0)
                .map { id, marks, dists ->
                    def dist_map = [marks, dists].transpose().collectEntries { k, v -> [k, v] }
                    return [ id, dist_map ]
                }

            ch_best_fit_model_input = ch_ready_for_analysis
                .map { tuple -> [ tuple[0].id, tuple ] } 
                .join(ch_updated_dist)                  
                .map { id, original_tuple, dist_map ->
                    def new_meta = original_tuple[0].clone()
                    new_meta.distributions = dist_map    
                    return [ new_meta ] + original_tuple.drop(1) 
                }

            if (params.duration) {
                MODEL_TRAINING_DM(ch_best_fit_model_input)
                ch_versions = ch_versions.mix(MODEL_TRAINING_DM.out.versions)
            } else {
                MODEL_TRAINING_STD(ch_best_fit_model_input)
                ch_versions = ch_versions.mix(MODEL_TRAINING_STD.out.versions)
            }
        }

    } else if (active_mode == 'duration') {
        MODEL_TRAINING_DM(ch_ready_for_analysis)
        ch_versions = ch_versions.mix(MODEL_TRAINING_DM.out.versions)

    } else if (active_mode == 'DNA') {
        ch_dna_model_input = ch_ready_for_analysis.map { [ it[0], it[1], it.size() > 2 ? it[2] : params.states ] }
        MODEL_TRAINING_DNA(ch_dna_model_input)
        ch_versions = ch_versions.mix(MODEL_TRAINING_DNA.out.versions)

    } else {
        MODEL_TRAINING_STD(ch_ready_for_analysis)
        ch_versions = ch_versions.mix(MODEL_TRAINING_STD.out.versions)
    }

    // ---------------------------------------------------------
    // MULTIQC & VERSION COLLATE
    // ---------------------------------------------------------

    // Collate Software Versions
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_episegmix_software_mqc_versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    // MODULE: MultiQC
    summary_params      = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)

    MULTIQC (
        ch_multiqc_files.collect(),
        [], [], [], [], []
    )

    emit:
    versions = ch_versions 
}