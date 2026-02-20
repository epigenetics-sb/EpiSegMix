process EPISEGMIX_DM_PREPARE {
    tag "Prep: ${meta.id} (State: ${state})"
    label 'process_medium'

    input:
    tuple val(meta), path(histone), path(meth), val(state)
    val(chr_params)

    output:
    tuple val(meta), path("*.yaml"), emit: config
    tuple val(meta), path("*.trainCounts.txt"), emit: train_counts
    tuple val(meta), path("*.regions.txt"), emit: regions
    tuple val(meta), path("*.trainCountsMeth.txt"), emit: train_meth

    script:
    def prefix = "${meta.id}"
    """
    set -euo pipefail

    export MPLCONFIGDIR=\$(pwd)
    ln -s /app/src src
    WF_DIR=\$(pwd)

    HAS_METH=0
    if [ -s "${meth}" ]; then
        HAS_METH=1
    fi

    HEADER=\$(head -n 1 "${histone}")
    FILE_MARKS=(\$(echo "\$HEADER" | cut -f4-))
    MARKER_SPEC=""

    for M in "\${FILE_MARKS[@]}"; do
        MARKER_SPEC="\${MARKER_SPEC}  - name: \${M}\n    distribution: ${params.dist_histone}\n"
    done

    cat <<EOF > "${prefix}.yaml"
states: ${state}
marker: \${#FILE_MARKS[@]}
marker_spec:
\${MARKER_SPEC}data: [\$WF_DIR/${histone}]
chr: ${chr_params}
EOF

    if [ "\$HAS_METH" -eq 1 ]; then
        cat <<EOF >> "${prefix}.yaml"
dna_methylation: ${params.dist_methyl}
meth_data: [\$WF_DIR/${meth}]
EOF
    fi

    python /app/src/get_counts.py \\
        -d "\$WF_DIR/${prefix}.yaml" \\
        -c "\$WF_DIR/${prefix}.trainCounts.txt" \\
        -m "\$WF_DIR/${prefix}.trainCountsMeth.txt" \\
        -r "\$WF_DIR/${prefix}.regions.txt"
    """
}
