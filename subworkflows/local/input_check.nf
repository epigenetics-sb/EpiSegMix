include { samplesheetToList } from 'plugin/nf-schema'

workflow INPUT_CHECK {
    
    main:

    def input_file  = file(params.input, checkIfExists: true)
    def schema_file = file(params.schema, checkIfExists: true)

    def raw_samples = samplesheetToList(input_file, schema_file)

    ch_samples = Channel.fromList(raw_samples)
        .map { sample_id, replicate, epigenetic_mark, file_path, modality, paired_end ->
            def meta = [
                id: sample_id
            ]
            [ 
                meta, 
                replicate, 
                epigenetic_mark, 
                file(file_path, checkIfExists: true), 
                modality, 
                paired_end 
            ]
        }

    ch_samples
        .branch { meta, rep, mark, f, mod, pe ->
            def fname = f.getName().toLowerCase()
            bam: fname.endsWith('.bam') || fname.endsWith('.bam.gz') || fname.endsWith('.sam')
            bed: fname.endsWith('.bed') || fname.endsWith('.bed.gz') || fname.endsWith('.tab') || fname.endsWith('.txt.gz')
            unknown: true 
        }
        .set { ch_split }
        
    ch_split.unknown.view { meta, rep, mark, f, mod, pe -> 
        "WARNING: Ignoring unrecognized file type for sample ${meta.id}: ${f}" 
    }

    emit:
    bam = ch_split.bam
    bed = ch_split.bed
}
