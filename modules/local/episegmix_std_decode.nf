process EPISEGMIX_STD_DECODE {
    tag "Decode: ${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(config), path(model)

    output:
    tuple val(meta), path("segmentation/${meta.id}.bed.gz"), emit: bed
    tuple val(meta), path("segmentation/${meta.id}.txt")   , emit: seg_txt

    script:
    """
    set -euo pipefail
    export MPLCONFIGDIR=\$(pwd)
    ln -s /app/src src
    
    mkdir -p counts_${meta.id} states_${meta.id} segmentation

    python /app/src/get_counts_for_all.py -d "${config}" -o "counts_${meta.id}"

    for file in \$(find "counts_${meta.id}" -name "counts*" -type f -print); do
        
        name=\${file##*/counts_}
        
        /app/HMM/build/HMMChromSeg \\
            -m "${model}" \\
            -v "states_${meta.id}/viterbi_\$name" \\
            -d "states_${meta.id}/posterior_\$name" \\
            -c "\$file" \\
            -x "counts_${meta.id}/meth_\$name" \\
            -r "counts_${meta.id}/regions_\$name" \\
            -p ${task.cpus}
            
    done

    python /app/src/segmentation_to_bed.py -d "${config}" -i "states_${meta.id}/" -o "segmentation/" -c ${params.decoding_algorithm}

    mv segmentation/*.bed.gz segmentation/${meta.id}.bed.gz
    mv segmentation/*.txt    segmentation/${meta.id}.txt
    """
}