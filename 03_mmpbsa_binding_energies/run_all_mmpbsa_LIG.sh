#!/bin/bash

# --- Configuration ---
# Define the pattern for the directories containing the gromacs folder.
# This assumes the structure: ZINC_ID/POSE_TAG/gromacs
TARGET_DIRS=$(find . -maxdepth 3 -type d -path "./ZINC*/*/gromacs")
TARGET_DIRS+=" $(find . -maxdepth 3 -type d -path "./4zzj*/*/gromacs")"
TARGET_DIRS+=" $(find . -maxdepth 3 -type d -path "./4zzj*/gromacs")" # For 4zzj_nad_p53 which might be one level shorter

# Remove duplicates and clean up the list of directories
TARGET_DIRS=$(echo "$TARGET_DIRS" | tr ' ' '\n' | sort -u)

# --- Execution ---
echo "Starting gmx_MMPBSA calculations across directories..."
echo "----------------------------------------------------"

# Loop through each unique gromacs directory found
for GMX_DIR in $TARGET_DIRS; do
    # Check if the required input files exist in the directory before running
    if [ -f "$GMX_DIR/mmpbsa.in" ] && [ -f "$GMX_DIR/step5_production.tpr" ]; then
        
        echo -e "\nProcessing directory: $GMX_DIR"
        
        # Change into the gromacs directory
        cd "$GMX_DIR" || { echo "Error: Failed to change directory to $GMX_DIR. Skipping."; continue; }

        # --- GMX_MMPBSA COMMAND EXECUTION ---
        # The command is split across multiple lines for readability
        gmx_MMPBSA -O -i mmpbsa.in \
                   -cs step5_production.tpr \
                   -ci new_index.ndx \
                   -cg 1 17 \
                   -ct traj_whole_centered.xtc \
                   -cp topol.top \
                   -o FINAL_RESULTS_MMPBSA.dat
        
        # Capture the exit status of the previous command
        EXIT_STATUS=$?

        # Check for successful execution
        if [ $EXIT_STATUS -eq 0 ]; then
            echo "SUCCESS: gmx_MMPBSA completed successfully in $GMX_DIR"
        else
            echo "FAILURE: gmx_MMPBSA returned error code $EXIT_STATUS in $GMX_DIR"
        fi

        # Change back to the starting directory to prepare for the next loop iteration
        cd - > /dev/null
    else
        echo "Skipping $GMX_DIR: Missing required input files (mmpbsa.in or step5_production.tpr)."
    fi
done

echo -e "\n----------------------------------------------------"
echo "Batch processing complete."
