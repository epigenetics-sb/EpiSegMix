process EPISEGMIX_STD_REPORT {
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

    # Generate statistics plots
    python /app/src/plot_statistics.py \\
        -d "${config}" \\
        -p "plots/${meta.id}/${meta.id}-histogram.png" \\
        -c "plots/${meta.id}/${meta.id}-correlation.png" \\
        -m "plots/${meta.id}/${meta.id}-methylation-density.png"

    # Generate main segmentation plots
    python /app/src/results.py \\
        -c "${seg_txt}" \\
        -j "${model}" \\
        -e "plots/${meta.id}/${meta.id}-meanEmission-${params.decoding_algorithm}.png" \\
        -t "plots/${meta.id}/${meta.id}-transitionMatrix.png" \\
        -m "plots/${meta.id}/${meta.id}-stateMembership-${params.decoding_algorithm}.png" \\
        -l "plots/${meta.id}/${meta.id}-stateLength-${params.decoding_algorithm}.png" \\
        -s ${params.decoding_algorithm} \\
        -n "plots/${meta.id}/${meta.id}-normEmission-${params.decoding_algorithm}.png" \\
        -d "${bed}"

    # Plot state distributions
    python /app/src/plot_state_histograms.py \\
        -c "${seg_txt}" \\
        -j "${model}" \\
        -a "plots/${meta.id}/${meta.id}-stateDistribution-${params.decoding_algorithm}.png" \\
        -s ${params.decoding_algorithm}

    # Plot state colors
    python /app/src/plot_state_colors.py \\
        -d "${bed}" \\
        -o "plots/${meta.id}/${meta.id}-state-colors.png"

    # Generate HTML report
    bash /app/workflow/segmentation_report.sh \\
        -n "${meta.id}" \\
        -o "plots/${meta.id}" \\
        -i true \\
        -c ${params.decoding_algorithm}
    """
}
