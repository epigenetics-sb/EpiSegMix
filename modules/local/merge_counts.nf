process MERGE_COUNTS {
    tag "Merge: ${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta),            
          val(histone_reps),    
          val(histone_marks),   
          path(bam_files),      
          val(meth_modalities), 
          path(bam_idxs),       
          path(histone_counts), 
          path(meth_tabs)       
    val   genome
    path  ref                   
    path  segment_coords        

    output:
    tuple val(meta), 
          path("${meta.id}_seg_tabs_counts.bed"), 
          emit: merged_counts

    script:
    def sample_id = meta.id
    def wgbs_path = ""
    def nome_path = ""
    
    def mods_list = meth_modalities instanceof List ? meth_modalities : [meth_modalities]
    def tabs_list = meth_tabs instanceof List ? meth_tabs : [meth_tabs]

    mods_list.eachWithIndex { mod, i ->
        if (mod == 'WGBS') { wgbs_path = tabs_list[i] }
        else if (mod == 'NOMe-seq') { nome_path = tabs_list[i] }
    }

    """
    set -euo pipefail

    echo "Processing ${sample_id}..."

    SORTED_SEG="${sample_id}_segments_sorted.bed"
    bedtools sort -i "${segment_coords}" > "\$SORTED_SEG"

    CURRENT_BED="\$SORTED_SEG"
    HAS_WGBS="false"
    HAS_NOME="false"

    if [[ -n "${wgbs_path}" ]]; then
        echo "Mapping WGBS: ${wgbs_path}"
        
        bedtools map \\
            -a "\$CURRENT_BED" \\
            -b <(tail -n +2 "${wgbs_path}") \\
            -c 4,5 -o median -null 0 \\
            > temp_wgbs.bed
        
        mv temp_wgbs.bed "\$CURRENT_BED"
        HAS_WGBS="true"
    fi

    if [[ -n "${nome_path}" ]]; then
        echo "Mapping NOME: ${nome_path}"
        
        bedtools map \\
            -a "\$CURRENT_BED" \\
            -b <(tail -n +2 "${nome_path}") \\
            -c 4,5 -o median -null 0 \\
            > temp_nome.bed
        
        mv temp_nome.bed "\$CURRENT_BED"
        HAS_NOME="true"
    fi

    FINAL_OUT="${sample_id}_seg_tabs_counts.bed"
    INTERSECT_TEMP="intersect_temp.bed"
    
    HISTONE_HEADER=\$(head -n 1 "${histone_counts}")

    bedtools intersect \\
        -a "\$CURRENT_BED" \\
        -b <(tail -n+2 "${histone_counts}") \\
        -wa -wb \\
        > "\$INTERSECT_TEMP"

    if [[ "\$HAS_WGBS" == "true" && "\$HAS_NOME" == "true" ]]; then
        echo "Format: WGBS + NOME"
        echo -e "\${HISTONE_HEADER}\tCov_WGBS\tMeth_WGBS\tCov_NOME\tMeth_NOME" > "\$FINAL_OUT"
        
        awk -v OFS='\\t' '{
            for(i=8; i<=NF; i++) printf "%s%s", \$i, OFS;
            printf "%s\\t%s\\t%s\\t%s\\n", \$4, \$5, \$6, \$7;
        }' "\$INTERSECT_TEMP" >> "\$FINAL_OUT"

    elif [[ "\$HAS_WGBS" == "true" ]]; then
        echo "Format: WGBS Only"
        echo -e "\${HISTONE_HEADER}\tCov_WGBS\tMeth_WGBS" > "\$FINAL_OUT"
        
        awk -v OFS='\\t' '{
            for(i=6; i<=NF; i++) printf "%s%s", \$i, OFS;
            printf "%s\\t%s\\n", \$4, \$5;
        }' "\$INTERSECT_TEMP" >> "\$FINAL_OUT"

    elif [[ "\$HAS_NOME" == "true" ]]; then
        echo "Format: NOME Only"
        echo -e "\${HISTONE_HEADER}\tCov_NOME\tMeth_NOME" > "\$FINAL_OUT"
        
        awk -v OFS='\\t' '{
            for(i=6; i<=NF; i++) printf "%s%s", \$i, OFS;
            printf "%s\\t%s\\n", \$4, \$5;
        }' "\$INTERSECT_TEMP" >> "\$FINAL_OUT"

    else
        echo "Error: No methylation inputs processed."
        exit 1
    fi
    """
}
