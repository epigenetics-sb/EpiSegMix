process EPISEGMIX_DNA_PREPARE {

    tag "Prep: ${meta.id} (State: ${state})"
    label 'process_medium'

    input:
    tuple val(meta), path(meth), val(state)

    output:
    tuple val(meta), path("${meta.id}.yaml"), emit: config
    tuple val(meta), path("${meta.id}.trainCounts.txt"), emit: train_counts
    tuple val(meta), path("${meta.id}.trainregions.txt"), emit: train_regions
    tuple val(meta), path("${meta.id}.counts.txt"), emit: counts
    tuple val(meta), path("${meta.id}.regions.txt"), emit: regions

    script:
    def prefix = meta.id

    def meth_list = meth instanceof List ? meth : [meth]

    def yaml_data
    if (meth_list.size() == 1) {
        yaml_data = "\"\$WF_DIR/${meth_list[0]}\""
    } else {
        yaml_data = meth_list.collect { "\"\$WF_DIR/${it}\"" }.join(", ")
        yaml_data = "[${yaml_data}]"  // wrap as YAML array
    }

    """
    set -euo pipefail

    # Ensure matplotlib can run inside containers
    export MPLCONFIGDIR=\$(pwd)

    # Symlink source code safely
    ln -sfn /app/src src
    WF_DIR=\$(pwd)

    # Validate methylation files exist
    for f in ${meth_list.collect { "\$WF_DIR/${it}" }.join(" ")}; do
        if [[ ! -f "\$f" ]]; then
            echo "[ERROR] Methylation file not found: \$f"
            exit 1
        fi
    done

    # Create YAML config
    cat <<EOF > "${prefix}.yaml"
distribution: ${params.dist_methyl}
states: ${state}
data: ${yaml_data}
marker: WGBS
EOF

    echo "[INFO] YAML created: ${prefix}.yaml"

    # Run Python script
    python src/get_meth_counts.py \\
        -d "\$WF_DIR/${prefix}.yaml" \\
        -c "\$WF_DIR/${prefix}.trainCounts.txt" \\
        -r "\$WF_DIR/${prefix}.trainregions.txt" \\
        -C "\$WF_DIR/${prefix}.counts.txt" \\
        -R "\$WF_DIR/${prefix}.regions.txt"

    echo "[INFO] Methylation counts generated for ${prefix}"
    """
}
