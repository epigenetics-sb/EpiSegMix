process SPLIT_MERGED_COUNTS {
    tag "Split: ${meta.id}"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://biocontainers/biocontainers:v1.2.0_cv1' :
        'biocontainers/biocontainers:v1.2.0_cv1' }"

    input:
    tuple val(meta), path(merged_bed)

    output:
    tuple val(meta), 
          path("${meta.id}_Histone_Input.txt"), 
          path("${meta.id}_Meth_Input.txt"), 
          emit: inputs
    path "versions.yml", emit: versions

    script:
    def sample_id = meta.id
    """
    set -euo pipefail

    HIST_OUT="${sample_id}_Histone_Input.txt"
    METH_OUT="${sample_id}_Meth_Input.txt"

    : > "\$METH_OUT"
    : > "\$HIST_OUT"

    awk -v hist_out="\$HIST_OUT" \\
        -v meth_out="\$METH_OUT" \\
        '
    BEGIN { OFS="\\t" }

    NR==1 {
        wgbs_c=0; wgbs_m=0; nome_c=0; nome_m=0;

        for(i=1; i<=NF; i++) {
            if(\$i == "Cov_WGBS")  wgbs_c=i;
            if(\$i == "Meth_WGBS") wgbs_m=i;
            if(\$i == "Cov_NOME")  nome_c=i;
            if(\$i == "Meth_NOME") nome_m=i;
        }

        first_meth_col = 9999;
        if (wgbs_c > 0 && wgbs_c < first_meth_col) first_meth_col = wgbs_c;
        if (nome_c > 0 && nome_c < first_meth_col) first_meth_col = nome_c;

        if (first_meth_col == 9999) first_meth_col = NF + 1;
        histone_end = first_meth_col - 1;

        for(i=1; i<=histone_end; i++)
            printf "%s%s", \$i, (i==histone_end ? "" : OFS) > hist_out;

        if(nome_c > 0)
            printf "%sNOME", OFS > hist_out;

        printf "\\n" > hist_out;

        if(wgbs_c > 0) {
            print "Cov\\tMeth" > meth_out;
        }
    }

    NR > 1 {
        for(i=1; i<=histone_end; i++)
            printf "%s%s", \$i, (i==histone_end ? "" : OFS) > hist_out;

        if(nome_c > 0) {
            cov  = \$(nome_c);
            meth = \$(nome_m);
            pct  = (cov > 0) ? (meth / cov)*100 : 0;
            printf "%s%.4f", OFS, pct > hist_out;
        }

        printf "\\n" > hist_out;

        if (wgbs_c > 0) {
            print \$(wgbs_c), \$(wgbs_m) > meth_out;
        }
    }
    ' "${merged_bed}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version | head -n 1 | awk '{print \$3}')
    END_VERSIONS
    """

    stub:
    def sample_id = meta.id

    """
    touch "${sample_id}_Histone_Input.txt"
    touch "${sample_id}_Meth_Input.txt"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version | head -n 1 | awk '{print \$3}')
    END_VERSIONS
    """
}