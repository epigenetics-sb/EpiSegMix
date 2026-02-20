process GENERATE_COUNT_MATRIX_BAM {
    tag "Counts: ${meta.id}"
    label 'process_high'

    input:
    tuple val(meta),        
          path(bams),       
          path(bais),       
          val(replicates), 
          val(marks), 
          val(modalities), 
          val(paired_ends)  
    val   genome
    path  reference_file 
    
    output:
    tuple val(meta), 
          path("${out_dir}/*_refined_counts.txt"), 
          emit: counts

    script:
    def sample_id = meta.id
    out_dir = "${sample_id}_Histones"
    
    def pe_lines = []
    def se_lines = []

    marks.eachWithIndex { mark, index ->
        def bam_file = bams[index]
        def is_pe    = paired_ends[index]
        
        def line = "${mark}\t${bam_file}"
        if (is_pe) { pe_lines.add(line) } 
        else       { se_lines.add(line) }
    }

    def pe_content = pe_lines.join('\n')
    def se_content = se_lines.join('\n')

    """
    set -euo pipefail

    mkdir -p "${out_dir}"
    
    source /opt/conda/etc/profile.d/conda.sh
    conda activate counts
    
    cp /app/* .

    PE_OUT="${out_dir}/${sample_id}_PE_${genome}_refined_counts.txt"
    SE_OUT="${out_dir}/${sample_id}_SE_${genome}_refined_counts.txt"
    FINAL_OUT="${out_dir}/${sample_id}_${genome}_refined_counts.txt"

    awk -v OFS="\t" '{print \$1, 1, \$2}' "${reference_file}" > new_reference

    if [[ -n "${pe_content}" ]]; then
        echo "Running Paired-End counts..."
        echo "${pe_content}" > "PE_info.txt"
        bash counts.sh \\
            -t "PE_info.txt" \\
            -o "${out_dir}" \\
            -c ${task.cpus} \\
            -p midpoint \\
            -f "${sample_id}_PE" \\
            -g "${genome}" \\
            -r ./new_reference \\
            -b ${params.binsize}
    fi

    awk -v OFS="\t" '{print \$1, 1, \$2}' "${reference_file}" > new_reference

    if [[ -n "${se_content}" ]]; then
        echo "Running Single-End counts..."
        echo "${se_content}" > "SE_info.txt"
        bash counts.sh \\
            -t "SE_info.txt" \\
            -o "${out_dir}" \\
            -c ${task.cpus} \\
            -p ignore \\
            -f "${sample_id}_SE" \\
            -g "${genome}" \\
            -r ./new_reference \\
            -b ${params.binsize}
    fi

    if [[ -f "\$PE_OUT" && -f "\$SE_OUT" ]]; then
        echo "Merging PE and SE outputs..."
        paste "\$PE_OUT" <(cut -f4- "\$SE_OUT") > "\$FINAL_OUT"
        rm "\$PE_OUT" "\$SE_OUT"

    elif [[ -f "\$PE_OUT" ]]; then
        mv "\$PE_OUT" "\$FINAL_OUT"

    elif [[ -f "\$SE_OUT" ]]; then
        mv "\$SE_OUT" "\$FINAL_OUT"
        
    else
        echo "Error: No output files generated for ${sample_id}"
        ls -R
        exit 1
    fi

    rm -f Dockerfile env.yaml counts.* new_reference *_info.txt
    """
}
