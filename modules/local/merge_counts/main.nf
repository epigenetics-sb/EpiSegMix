process MERGE_COUNTS {
    tag "Merge: ${meta.id}"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://community.wave.seqera.io/library/bedtools:2.31.1--7c4ce4cb07c09ee4' :
        'community.wave.seqera.io/library/bedtools:2.31.1--7c4ce4cb07c09ee4' }"

    input:
    tuple val(meta), 
          val(meta_list_meth), 
          path(files_meth), 
          val(meta_list_hist), 
          path(files_hist) 
    val   genome
    path  ref            
    path  segment_coords 

    output:
    tuple val(meta), 
          path("${meta.id}_seg_tabs_counts.bed"), 
          emit: merged_counts
    path "versions.yml", emit: versions

    script:
    def sample_id = meta.id
    def wgbs_path = ""
    def nome_path = ""

    meta_list_meth.eachWithIndex { m, i ->
        def current_file = files_meth instanceof List ? files_meth[i].name : files_meth.name
        if (m.modality == 'WGBS') { wgbs_path = current_file } 
        else if (m.modality == 'NOMe-seq') { nome_path = current_file }
    }

    def histone_path = files_hist instanceof List ? files_hist[0] : files_hist

    """
    set -euo pipefail

    SORTED_SEG="${sample_id}_segments_sorted.bed"
    bedtools sort -i "${segment_coords}" > "\$SORTED_SEG"

    CURRENT_BED="\$SORTED_SEG"
    HAS_WGBS="false"
    HAS_NOME="false"

    if [[ -n "${wgbs_path}" ]]; then
        bedtools map \\
            -a "\$CURRENT_BED" \\
            -b <(tail -n +2 "${wgbs_path}") \\
            -c 4,5 -o median -null 0 \\
            > temp_wgbs.bed
        mv temp_wgbs.bed "\$CURRENT_BED"
        HAS_WGBS="true"
    fi

    if [[ -n "${nome_path}" ]]; then
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
    HISTONE_HEADER=\$(head -n 1 "${histone_path}")

    bedtools intersect \\
        -a "\$CURRENT_BED" \\
        -b <(tail -n+2 "${histone_path}") \\
        -wa -wb \\
        > "\$INTERSECT_TEMP"

    if [[ "\$HAS_WGBS" == "true" && "\$HAS_NOME" == "true" ]]; then
        echo -e "\${HISTONE_HEADER}\tCov_WGBS\tMeth_WGBS\tCov_NOME\tMeth_NOME" > "\$FINAL_OUT"
        awk -v OFS='\\t' '{
            for(i=8; i<=NF; i++) printf "%s%s", \$i, OFS;
            printf "%s\\t%s\\t%s\\t%s\\n", \$4, \$5, \$6, \$7;
        }' "\$INTERSECT_TEMP" >> "\$FINAL_OUT"

    elif [[ "\$HAS_WGBS" == "true" ]]; then
        echo -e "\${HISTONE_HEADER}\tCov_WGBS\tMeth_WGBS" > "\$FINAL_OUT"
        awk -v OFS='\\t' '{
            for(i=6; i<=NF; i++) printf "%s%s", \$i, OFS;
            printf "%s\\t%s\\n", \$4, \$5;
        }' "\$INTERSECT_TEMP" >> "\$FINAL_OUT"

    elif [[ "\$HAS_NOME" == "true" ]]; then
        echo -e "\${HISTONE_HEADER}\tCov_NOME\tMeth_NOME" > "\$FINAL_OUT"
        awk -v OFS='\\t' '{
            for(i=6; i<=NF; i++) printf "%s%s", \$i, OFS;
            printf "%s\\t%s\\n", \$4, \$5;
        }' "\$INTERSECT_TEMP" >> "\$FINAL_OUT"

    else
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
        awk: \$(awk --version | head -n 1 | awk '{print \$3}')
    END_VERSIONS
    """

    stub:
    def sample_id = meta.id

    """
    touch "${sample_id}_seg_tabs_counts.bed"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
        awk: \$(awk --version | head -n 1 | awk '{print \$3}')
    END_VERSIONS
    """
}