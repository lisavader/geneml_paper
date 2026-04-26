#!/bin/bash
set -euo pipefail

# Default to running all steps if no arguments given
if [ $# -eq 0 ]; then
    STEPS="model_training benchmarking alt_transcripts reannotation figures"
else
    STEPS="$@"
fi

for step in $STEPS; do
    case $step in
        model_training)
            echo "=== Step 1: Model training ==="
            for script in scripts/01_model_training/*.sh; do
                echo "Executing $script..."
                bash $script
            done
            ;;
        benchmarking)
            echo "=== Step 2: Benchmarking ==="
            for script in scripts/02_benchmarking/*.sh; do
                echo "Executing $script..."
                bash $script
            done
            ;;
        alt_transcripts)
            echo "=== Step 3: Alternative transcript analysis ==="
            for script in scripts/03_alt_transcripts/*.sh; do
                echo "Executing $script..."
                bash $script
            done
            ;;
        reannotation)
            echo "=== Step 4: Training set reannotation ==="
            for script in scripts/04_reannotation/*.sh; do
                echo "Executing $script..."
                bash $script
            done
            ;;
        figures)
            echo "=== Step 5: Create figures / tables ==="

            mkdir -p figures tables
            python3 -m ipykernel install --sys-prefix --name python3

            for nb in scripts/05_figures/*.ipynb; do
                echo "Executing $nb..."
                jupyter nbconvert --to notebook --execute --inplace $nb
            done
            ;;
        *)
            echo "Unknown step name: $step"
            exit 1
            ;;
    esac
done

echo "=== Pipeline complete ==="
