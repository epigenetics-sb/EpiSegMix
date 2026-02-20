process GENERATE_GENOME_BINS {
    tag "Bins: ${genome} (${bin_size}bp)"
    label 'process_single'

    input:
    path chrom_sizes    
    val genome          
    val bin_size        

    output:
    path "${out_filename}", emit: bins_file

    script:
    out_filename = "${genome}_${bin_size}bp_bins.bed"
    """
    set -euo pipefail

    echo "Generating ${bin_size}bp bins for ${genome}..."

    bedtools makewindows \\
        -g "${chrom_sizes}" \\
        -w "${bin_size}" \\
    | awk -v size="${bin_size}" 'BEGIN{OFS="\\t"} (\$3-\$2) == size' \\
    > "${out_filename}"
    """
}
