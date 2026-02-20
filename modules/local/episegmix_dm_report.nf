process EPISEGMIX_DM_REPORT {
    tag "Report: ${meta.id}"
    label 'process_low'

    input:
    tuple val(meta), path(config), path(model), path(bed), path(seg_txt)

    output:
    path "plots/${meta.id}/*.png" , emit: plots
    path "plots/${meta.id}/*.html", emit: html_report

    script:
    """
    set -euo pipefail

    export MPLCONFIGDIR=\$(pwd)
    ln -s /app/src src
    mkdir -p plots/${meta.id}

    python /app/src/plot_statistics.py \\
        -d "${config}" \\
        -p "plots/${meta.id}/${meta.id}-histogram.png" \\
        -c "plots/${meta.id}/${meta.id}-correlation.png" \\
        -m "plots/${meta.id}/${meta.id}-methylation-density.png"

    python /app/src/results.py \\
        -c "${seg_txt}" \\
        -j "${model}" \\
        -e "plots/${meta.id}/${meta.id}-meanEmission.png" \\
        -t "plots/${meta.id}/${meta.id}-transitionMatrix.png" \\
        -m "plots/${meta.id}/${meta.id}-stateMembership.png" \\
        -l "plots/${meta.id}/${meta.id}-stateLength.png" \\
        -s viterbi \\
        -n "plots/${meta.id}/${meta.id}-normEmission.png" \\
        -d "${bed}"

    python /app/src/plot_state_histograms.py \\
        -c "${seg_txt}" \\
        -j "${model}" \\
        -a "plots/${meta.id}/${meta.id}-stateDistribution.png" \\
        -s viterbi

    python /app/src/plot_state_colors.py \\
        -d "${bed}" \\
        -o "plots/${meta.id}/${meta.id}-state-colors.png"

    bash /app/workflowTopology/segmentation_report.sh \\
        -n "${meta.id}" \\
        -o "plots/${meta.id}" \\
        -i true \\
    """
}
