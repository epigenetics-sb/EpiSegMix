process GET_GENOME_CHROM_SIZE {
    tag "Getting ${genome} chrom sizes"
    label 'process_single'

    input:
    val genome          

    output:
    path "*.chrom.sizes", emit: genome_chrom_size_file

    script:
    """
    set -euo pipefail

    RAW_FILE="${genome}.chrom.sizes"
    TMP_FILE="temp.chrom.sizes"
    URL="http://hgdownload.soe.ucsc.edu/goldenPath/${genome}/bigZips/\${RAW_FILE}"

    echo "Downloading from: \${URL}"
    curl -s -O "\${URL}"

    if [[ "${genome}" == hg* ]]; then
        echo "Processing Human (${genome}): Keeping 1-22, X, Y"
        
        awk 'BEGIN{OFS="\\t"} {
            sub(/^chr/, "", \$1);
            if (\$1 ~ /^(1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16|17|18|19|20|21|22|X|Y)\$/) print
        }' "\${RAW_FILE}" > "\${TMP_FILE}"

    elif [[ "${genome}" == mm* ]]; then
        echo "Processing Mouse (${genome}): Keeping 1-19, X, Y"
        
        awk 'BEGIN{OFS="\\t"} {
            sub(/^chr/, "", \$1);
            if (\$1 ~ /^(1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16|17|18|19|X|Y)\$/) print
        }' "\${RAW_FILE}" > "\${TMP_FILE}"

    else
        echo "Warning: Unrecognized genome '${genome}'. Removing 'chr' only."
        awk 'BEGIN{OFS="\\t"} { 
            sub(/^chr/, "", \$1); 
            print 
        }' "\${RAW_FILE}" > "\${TMP_FILE}"
    fi

    sort -k1,1V "\${TMP_FILE}" > "\${RAW_FILE}"
    
    rm "\${TMP_FILE}"
    """
}
