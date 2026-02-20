process CLEAN_AND_INDEX_BAM {
    tag "Clean & Index: ${meta.id} | ${modality}"
    label 'process_medium'

    input:
    tuple val(meta), 
          val(replicate), 
          val(mark), 
          path(bam), 
          val(modality), 
          val(paired_end)

    output:
    tuple val(meta), 
          path("*.nochr.bam"),     
          path("*.nochr.bam.bai"), 
          val(replicate),    
          val(mark),         
          val(modality),     
          val(paired_end),   
          emit: sanitized_bam

    script:
    def out_bam = "${bam.baseName}.nochr.bam"
    
    """
    set -euo pipefail

    echo "[INFO] Processing ${bam} (Modality: ${modality})"

    if samtools view -H "${bam}" | grep "SN:chr" > /dev/null; then
        echo "[INFO] 'chr' prefix detected. Reheading..."
        
        samtools view -H "${bam}" \\
            | sed 's/SN:chr/SN:/g' \\
            | samtools reheader - "${bam}" > "${out_bam}"
    else
        echo "[INFO] No 'chr' prefix. Symlinking..."
        ln -s "${bam}" "${out_bam}"
    fi

    samtools index -@ ${task.cpus} "${out_bam}"
    """
}
