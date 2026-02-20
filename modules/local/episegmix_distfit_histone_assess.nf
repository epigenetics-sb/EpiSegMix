process EPISEGMIX_DISTFIT_HISTONE_ASSESS {
    tag "${meta.id} - Evaluate"
    label 'process_medium'

    input:
    tuple val(meta), path(histone_data), path(models)

    output:
    tuple val(meta), path("DISTFIT_${meta.id}.yaml"), emit: final_config

    script:
    """
    set -euo pipefail
    ln -sfn /app/src src

    RAW_HEADER=\$(head -n 1 ${histone_data})
    export HEADER_STR="\$RAW_HEADER"
    export HISTONE_ABS=\$(readlink -f ${histone_data})

    MARKERS=\$(python -c "import os; print(' '.join(os.environ.get('HEADER_STR', '').strip().split()[3:]))")

    mkdir -p temp_counts
    touch log_likelihoods.txt

    # Generate counts per marker
    for MARK in \$MARKERS; do
        cat <<EOF > dummy_\${MARK}.yaml
states: 3
marker: 1
marker_spec:
  - name: \${MARK}
    distribution: NBI
data: [\$HISTONE_ABS]
EOF
        
        mkdir -p temp_counts/\${MARK}
        python src/get_counts_for_all.py -d dummy_\${MARK}.yaml -o temp_counts/\${MARK}
    done

    # Compute log likelihoods
    for model in ${models.join(' ')}; do
        fname=\$(basename "\$model" .final-model.json)
        dist=\${fname%%-*}
        mark=\${fname#*-model-}

        for cfile in temp_counts/\${mark}/counts*.txt; do
            [ -f "\$cfile" ] || continue
            base=\$(basename \$cfile)
            suffix="\${base#counts_}"
            rfile="temp_counts/\${mark}/regions_\$suffix"

            if [ -f "\$rfile" ]; then
                echo -n -e "\$dist\t\$mark\t" >> log_likelihoods.txt
                /app/HMM/build/LogLikelihood \\
                    -m "\$model" \\
                    -c "\$cfile" \\
                    -r "\$rfile" >> log_likelihoods.txt || true
                echo "" >> log_likelihoods.txt
            fi
        done
    done

    # Create final YAML
    python - <<EOF
import pandas as pd
import yaml
import os

data = []

if os.path.exists('log_likelihoods.txt') and os.path.getsize('log_likelihoods.txt') > 0:
    with open('log_likelihoods.txt') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 3:
                try:
                    score = float(parts[-1])
                    data.append({'Mark': parts[1], 'Dist': parts[0], 'Score': score})
                except ValueError:
                    continue

num_states = 3

if data:
    df = pd.DataFrame(data)
    best_models = df.loc[df.groupby('Mark')['Score'].idxmax()]
    
    model = {}
    model['states'] = num_states
    model['marker'] = len(best_models)
    
    marker_spec = []
    for _, row in best_models.iterrows():
        marker_spec.append({'name': row['Mark'], 'distribution': row['Dist']})
        
    model['marker_spec'] = marker_spec
    model['data'] = [os.environ.get('HISTONE_ABS')]
    
    with open("DISTFIT_${meta.id}.yaml", 'w') as file:
        yaml.dump(model, file, sort_keys=False)
else:
    with open("DISTFIT_${meta.id}.yaml", 'w') as file:
        file.write("# Error: No log likelihood data generated\\n")
EOF
    """
}
