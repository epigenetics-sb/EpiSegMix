process EPISEGMIX_STD_TRAIN {
    tag "Train: ${meta.id}"
    label 'process_high'

    input:
    tuple val(meta), path(config), path(counts), path(regions), path(meth_counts)

    output:
    tuple val(meta), path("${meta.id}.model.json"), emit: model
    path "${meta.id}.log"                         , emit: log

    script:
    def meth_arg = meth_counts ? "-e ${meth_counts}" : ""
    def train_meth_arg = meth_counts ? "-x ${meth_counts}" : ""

    """
    set -euo pipefail
    export MPLCONFIGDIR=\$(pwd)
    ln -s /app/src src
    
    python /app/src/init_HMM.py \\
        -d "${counts}" \\
        ${meth_arg} \\
        -m "${config}" \\
        -j "${meta.id}.init.json"
    
    /app/HMM/build/HMMChromSeg \\
        -t \\
        -m "${meta.id}.init.json" \\
        -o "${meta.id}.model.json" \\
        -c "${counts}" \\
        ${train_meth_arg} \\
        -r "${regions}" \\
        -i "${params.iter}" -e "${params.epsilon}" \\
        -p ${task.cpus} \\
        &> "${meta.id}.log"
    """
}
