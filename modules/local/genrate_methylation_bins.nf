process GENERATE_METHYLATION_BINS {

    tag "Generate methyl bins: ${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta),
          val(meth_modalities),
          path(meth_tabs)
    path  segment_coords

    output:
    tuple val(meta),
          path("${meta.id}_seg_tabs_counts.txt"),
          emit: merged_counts

    script:
    def sample_id = meta.id

    def mods = meth_modalities instanceof List ? meth_modalities : [meth_modalities]
    def tabs = meth_tabs        instanceof List ? meth_tabs        : [meth_tabs]

    def wgbs_idx  = mods.findIndexOf { it == 'WGBS' }
    def wgbs_path = wgbs_idx >= 0 ? tabs[wgbs_idx] : null

    """
    set -euo pipefail

    echo "[GENERATE_METHYLATION_BINS] Processing: ${sample_id}"

    SORTED_SEG="${sample_id}_segments_sorted.bed"
    FINAL_OUT="${sample_id}_seg_tabs_counts.txt"

    bedtools sort -i "${segment_coords}" > "\$SORTED_SEG"

    if [[ -z "${wgbs_path ?: ''}" ]]; then
        echo "[GENERATE_METHYLATION_BINS] ERROR: WGBS modality required but not found"
        exit 1
    fi

    echo "[GENERATE_METHYLATION_BINS] Mapping WGBS: ${wgbs_path}"

    # extract header and body
    head -n 1 "${wgbs_path}" > header.txt
    tail -n +2 "${wgbs_path}" > wgbs_body.bed

    # map + add header in one stream
    bedtools map \\
        -a "\$SORTED_SEG" \\
        -b wgbs_body.bed \\
        -c 4,5 \\
        -o median \\
        -null 0 |
    awk 'NR==FNR{print; next} {print}' header.txt - > "\$FINAL_OUT"
    """
}
