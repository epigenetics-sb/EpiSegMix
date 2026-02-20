process EPISEGMIX_DNA_REPORT {
    tag "Report: ${meta.id}"
    label 'process_low'

    input:
    tuple val(meta), path(config), path(model), path(bed), path(seg_txt)

    output:
    path "plots/${meta.id}/*.png", emit: plots
    path "plots/${meta.id}/*.html", emit: html_report

    script:
    """
    set -euo pipefail

    export MPLCONFIGDIR=\$(pwd)
    ln -s /app/src src
    mkdir -p plots/${meta.id}

    python /app/src/results.py \\
        -c "${seg_txt}" \\
        -j "${model}" \\
        -e "plots/${meta.id}/${meta.id}-meanEmission.png" \\
        -n "plots/${meta.id}/${meta.id}-normEmission.png" \\
        -t "plots/${meta.id}/${meta.id}-transitionMatrix.png" \\
        -m "plots/${meta.id}/${meta.id}-stateMembership.png" \\
        -l "plots/${meta.id}/${meta.id}-stateLength.png" \\
        -s viterbi \\
        -d ${bed}


    python /app/src/plot_state_colors.py \\
        -d "${bed}" \\
        -o "plots/${meta.id}/${meta.id}-stateColors.png"

    bash /app/workflowDNAMethylation/segmentation_report.sh \\
        -n "${meta.id}" \\
        -o "plots/${meta.id}"
    """
}