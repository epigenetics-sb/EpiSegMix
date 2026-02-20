process EPISEGMIX_DNA_TRAIN {
    tag "Train: ${meta.id}"
    label 'process_high'

    input:
    tuple val(meta), path(config), path(train_counts), path(train_regions)

    output:
    tuple val(meta), path("*.model.json"), emit: model
    path "*.log"                         , emit: log

    script:
    """
    set -euo pipefail
    export MPLCONFIGDIR=\$(pwd)
    ln -s /app/src src
    
    python /app/src/init_HMM.py \\
        -e "${train_counts}" \\
        -m "${config}" \\
        -j "${meta.id}.preinit.json"
    
    python /app/src/add_topology_to_init.py  \\
        --input "${meta.id}.preinit.json" \\
        --output "${meta.id}.init.json"

    /app/HMM/build/TopologyHMM \\
        -t -n 2 \\
        -m "${meta.id}.init.json" \\
        -o "${meta.id}.model.json" \\
        -x "${train_counts}" \\
        -r "${train_regions}" \\
        -i "${params.iter}" -e "${params.epsilon}" \\
        -p ${task.cpus} \
        &> "${meta.id}.log"
    """
}
