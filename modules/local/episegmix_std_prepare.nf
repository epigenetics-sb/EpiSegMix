process EPISEGMIX_STD_PREPARE {
    tag "Prep: ${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(histone), path(meth), val(state)
    val(chr_params)

    output:
    tuple val(meta), path("${meta.id}.yaml"), emit: config
    tuple val(meta), path("${meta.id}.trainCounts.txt"), emit: train_counts
    tuple val(meta), path("${meta.id}.trainRegions.txt"), emit: regions
    tuple val(meta), path("${meta.id}.trainCountsMeth.txt"), emit: train_meth

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
        -d "${prefix}.yaml" \\
        -c "${prefix}.trainCounts.txt" \\
        -m "${prefix}.trainCountsMeth.txt" \\
        -r "${prefix}.trainRegions.txt"
    """
}
