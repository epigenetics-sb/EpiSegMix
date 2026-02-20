process EPISEGMIX_DISTFIT_HISTONE_TRAIN {
    tag "${meta.id} | ${mark} | ${dist}"
    label 'process_medium'

    input:
    tuple val(meta), path(histone_data), val(mark), val(dist)

    output:
    tuple val(meta), path("*final-model.json"), emit: model
    tuple val(meta), val(mark), val(dist),      emit: meta_info 

    script:
    def prefix = "${dist}-model-${mark}"

    """
    set -euo pipefail
    ln -sfn /app/src src
    export MPLCONFIGDIR=\$(pwd)

    # 1. Config (Strictly Histone)
    cat <<EOF > ${prefix}.yaml
states: 3
marker: 1
marker_spec:
  - name: ${mark}
    distribution: ${dist}
data: [\$(readlink -f ${histone_data})]
EOF

    # 2. Get Counts
    python src/get_counts.py \\
        -d "${prefix}.yaml" \\
        -m "${prefix}.trainCountsMeth.txt" \\
        -c "${prefix}.trainCounts.txt" \\
        -r "${prefix}.trainRegions.txt"

    # 3. Init
    python src/init_HMM.py \\
        -d "${prefix}.trainCounts.txt" \\
        -e "${prefix}.trainCountsMeth.txt" \\
        -m "${prefix}.yaml" \\
        -j "${prefix}.preinit.json"

    # 4. Topology
    python -c "
import json
with open('${prefix}.preinit.json') as f:
    hmm = json.load(f)
N = int(hmm['states'])
hmm['topology'] = {str(s+1): [s] for s in range(N)}
with open('${prefix}.init.json', 'w') as f:
    json.dump(hmm, f, indent=4)
"

    # 5. Train
    /app/HMM/build/TopologyHMM \\
        -t -n ${params.adjustment} \\
        -m "${prefix}.init.json" \\
        -o "${prefix}.final-model.json" \\
        -c "${prefix}.trainCounts.txt" \\
        -r "${prefix}.trainRegions.txt" \\
        -i ${params.iter} -p ${task.cpus} -e ${params.epsilon}
    """
}