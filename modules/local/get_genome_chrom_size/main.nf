process GET_GENOME_CHROM_SIZE {
    tag "Getting ${genome} chrom sizes"
    label 'process_single'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://biocontainers/biocontainers:v1.2.0_cv1' :
        'biocontainers/biocontainers:v1.2.0_cv1' }"

    input:
    val genome          

    output:
    path "${genome}.chrom.sizes", emit: genome_chrom_size_file
    path "versions.yml"         , emit: versions

    script:
    """
    set -euo pipefail

    DOWNLOAD_FILE="downloaded.chrom.sizes"
    TMP_FILE="temp.chrom.sizes"
    FINAL_FILE="${genome}.chrom.sizes"
    
    URL="http://hgdownload.soe.ucsc.edu/goldenPath/${genome}/bigZips/\${FINAL_FILE}"

    curl -s -o "\${DOWNLOAD_FILE}" "\${URL}"

    if [[ "${genome}" == hg* ]]; then
        awk 'BEGIN{OFS="\\t"} {
            sub(/^chr/, "", \$1);
            if (\$1 ~ /^(1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16|17|18|19|20|21|22|X|Y)\$/) print
        }' "\${DOWNLOAD_FILE}" > "\${TMP_FILE}"

    elif [[ "${genome}" == mm* ]]; then
        awk 'BEGIN{OFS="\\t"} {
            sub(/^chr/, "", \$1);
            if (\$1 ~ /^(1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16|17|18|19|X|Y)\$/) print
        }' "\${DOWNLOAD_FILE}" > "\${TMP_FILE}"

    else
        awk 'BEGIN{OFS="\\t"} { 
            sub(/^chr/, "", \$1); 
            print 
        }' "\${DOWNLOAD_FILE}" > "\${TMP_FILE}"
    fi

    sort -k1,1V "\${TMP_FILE}" > "\${FINAL_FILE}"
    
    rm "\${DOWNLOAD_FILE}" "\${TMP_FILE}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        curl: \$(curl --version | head -n 1 | awk '{print \$2}')
        awk: \$(awk --version | head -n 1 | awk '{print \$3}')
    END_VERSIONS
    """

    stub:
    """
    touch "${genome}.chrom.sizes"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        curl: \$(curl --version | head -n 1 | awk '{print \$2}')
        awk: \$(awk --version | head -n 1 | awk '{print \$3}')
    END_VERSIONS
    """
}