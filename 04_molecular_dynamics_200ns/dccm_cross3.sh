#!/bin/bash

# --- Configuration ---
GMX_BIN="gmx_mpi"
ANALYSIS_GROUP="C-alpha"

echo "--- Starting Residue-level DCCM Analysis ---"

# Loop over ZINC directories
for ZINC_DIR in ZINC*/; do
    ZINC_ID=$(basename "$ZINC_DIR")
    echo "Processing $ZINC_ID..."

    # Loop over all subdirectories (poses) inside the ZINC folder
    for SUBDIR in "$ZINC_DIR"*/; do
        if [ -d "$SUBDIR/gromacs" ]; then
            GMX_DIR="$SUBDIR/gromacs"
            echo "Found GROMACS directory: $GMX_DIR"

            TPR_FILE="$GMX_DIR/step5_production.tpr"
            TRAJ_FILE="$GMX_DIR/traj_whole_centered.xtc"
            NDX_FILE="$GMX_DIR/new_index.ndx"

            if [ ! -f "$TPR_FILE" ] || [ ! -f "$TRAJ_FILE" ] || [ ! -f "$NDX_FILE" ]; then
                echo "WARNING: Essential files missing in $GMX_DIR. Skipping."
                continue
            fi

            # Output directory per subdirectory
            OUTPUT_DIR="$SUBDIR/DCCM_Analysis"
            mkdir -p "$OUTPUT_DIR"

            # Check if C-alpha exists in the index file
            if ! grep -q "\[ $ANALYSIS_GROUP \]" "$NDX_FILE"; then
                echo "WARNING: $ANALYSIS_GROUP group not found in $NDX_FILE. Skipping."
                continue
            fi

            echo "Running covar for $SUBDIR..."
            echo -e "$ANALYSIS_GROUP\n$ANALYSIS_GROUP" | \
            $GMX_BIN covar \
                -f "$TRAJ_FILE" \
                -s "$TPR_FILE" \
                -n "$NDX_FILE" \
                -xpm "$OUTPUT_DIR/covariance_matrix.xpm" \
                -o "$OUTPUT_DIR/total_variance.xvg" \
                -v "$OUTPUT_DIR/eigenvectors.trr" \
                -l "$OUTPUT_DIR/eigenvalues.log" \
                -ascii "$OUTPUT_DIR/covariance_matrix_ascii.dat"

            if [ $? -eq 0 ]; then
                echo "Covariance calculation finished for $SUBDIR."
            else
                echo "ERROR: Covariance calculation failed for $SUBDIR"
            fi

            # Optional: call Python residue-level DCCM plotting
            python3 gpt_plot4.py "$OUTPUT_DIR/covariance_matrix_ascii.dat" "$TPR_FILE" "$OUTPUT_DIR"
        fi
    done
    echo "----------------------------------------"
done

echo "--- All ZINC compounds processed. Residue-level DCCM analysis complete ---"

