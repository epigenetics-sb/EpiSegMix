process EPISEGMIX_DNA_DECODE {
    tag "Decode: ${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(config), path(model), path(counts), path(regions)

    output:
    tuple val(meta), path("segmentation/${meta.id}.bed.gz"), emit: bed
    tuple val(meta), path("segmentation/${meta.id}.txt")   , emit: seg_txt
    tuple val(meta), path("states_${meta.id}")             , emit: states_dir
    tuple val(meta), path("counts_${meta.id}")             , emit: counts_dir

    script:
    """
    set -euo pipefail
    export MPLCONFIGDIR=\$(pwd)
    ln -s /app/src src

    mkdir -p counts_${meta.id} states_${meta.id} segmentation
    
    /app/HMM/build/TopologyHMM \\
        -m "${model}" \\
        -v "states_${meta.id}/states.txt" \\
        -x ${counts} \\
        -r ${regions} \\
        -p ${task.cpus}

  
    python /app/src/segmentation_to_bed.py \\
        -d "${config}" \\
        -i "states_${meta.id}/states.txt" \\
        -o "segmentation/" \\
        -c viterbi

    mv segmentation/*.bed.gz segmentation/${meta.id}.bed.gz
    mv segmentation/*.txt    segmentation/${meta.id}.txt
    """
}