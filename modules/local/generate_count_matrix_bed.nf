process GENERATE_COUNT_MATRIX_BED {
    tag "Counts BED: ${meta.id} | ${modality}"
    label 'process_medium'

    input:
    tuple val(meta), 
          val(replicate), 
          val(mark), 
          path(file_path), 
          val(modality), 
          val(paired_end)
    
    output:
    tuple val(meta), 
          val(replicate), 
          val(mark), 
          path(file_path), 
          val(modality), 
          val(paired_end), 
          path("${meta.id}_${modality}/${meta.id}_${modality}.tab"), 
          emit: counts_bed

    script:
    def sample_id = meta.id
    """
    set -euo pipefail
    export LC_ALL=C

    echo "[INFO] Processing ${file_path} for modality ${modality}"

    mkdir -p "${sample_id}_${modality}"
    OUTPUT_FILE="${sample_id}_${modality}/${sample_id}_${modality}.tab"

    READ_CMD="cat"
    if [[ "${file_path}" == *.gz ]]; then
        READ_CMD="zcat"
    fi

    {
        echo -e "chr\tstart\tend\tCov\tMeth"
        
        join -1 1 -2 1 -a 1 -a 2 \\
            <( \${READ_CMD} "${file_path}" \\
                | awk -v OFS='\\t' '\$6=="+" {print \$1"_"\$2, \$1, \$2, \$6, \$5, (\$11*\$5)/100}' \\
                | sort -k1,1 ) \\
            <( \${READ_CMD} "${file_path}" \\
                | awk -v OFS='\\t' '\$6=="-" {print \$1"_"(\$2-1), \$1, (\$2-1), \$6, \$5, (\$11*\$5)/100}' \\
                | sort -k1,1 ) \\
        | sort -k2,2 -k3,3n \\
        | awk -v OFS='\\t' '{
            sub(/^chr/, "", \$2);
            printf ("%s\\t%d\\t%d\\t%d\\t%d\\n", \$2, \$3, \$3+1, \$5+\$10, \$6+\$11+0.5)
        }' \\
        | grep -E -v "random|GL|NC|M|hs|hap|Un|J|EBV|ph|L"

    } > "\${OUTPUT_FILE}"
    """
}
