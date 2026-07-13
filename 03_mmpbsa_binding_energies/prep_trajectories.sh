#!/bin/bash

# Define the base directory (assuming the script is run from here)
BASE_DIR=$(pwd)

# --- Find all ZINC directories, but exclude reference/utility files/folders ---
# The pattern ZINC* covers all your ZINC directories (e.g., ZINC01078623)
ZINC_LIGANDS=$(find "$BASE_DIR" -maxdepth 1 -type d -name "ZINC*" | sort)

# Check if any ZINC directories were found
if [ -z "$ZINC_LIGANDS" ]; then
    echo "Error: No ZINC* directories found in $BASE_DIR."
    exit 1
fi

echo "Starting trajectory preparation for ZINC ligands..."
echo "---"

# Loop through each main ZINC ligand directory
for LIGAND_DIR in $ZINC_LIGANDS; do
    LIGAND_NAME=$(basename "$LIGAND_DIR")
    echo "Processing Ligand: $LIGAND_NAME"

    # Loop through the 4 pose subdirectories inside the ligand directory
    for POSE_DIR in "$LIGAND_DIR"/*; do
        if [ -d "$POSE_DIR" ]; then
            POSE_NAME=$(basename "$POSE_DIR")
            GROMACS_PATH="$POSE_DIR/gromacs"

            # Check for the required gromacs directory
            if [ -d "$GROMACS_PATH" ]; then
                echo "  -> Entering: $POSE_NAME/gromacs"
                
                # Change to the gromacs directory
                cd "$GROMACS_PATH" || { echo "Error: Failed to enter $GROMACS_PATH. Skipping."; continue; }

                # Check for input files before running
                if [ ! -f "step5_production.tpr" ] || [ ! -f "step5_production.xtc" ] || [ ! -f "new_index.ndx" ]; then
                    echo "    !! Missing required input files in $GROMACS_PATH. Skipping pose."
                    cd "$BASE_DIR"
                    continue
                fi

                # --- 1st Command: Remove PBC Jumps ---
                # Input '0' (System) is piped to select the group
                echo "    1. Running gmx trjconv (nojump)..."
                echo "0" | gmx trjconv -s step5_production.tpr -f step5_production.xtc -o traj_nojump.xtc -n new_index.ndx -pbc nojump
                
                if [ $? -ne 0 ]; then
                    echo "    !! Error running first gmx trjconv. Skipping rest of pose."
                    cd "$BASE_DIR"
                    continue
                fi
                
                # --- 2nd Command: Center and Make Whole ---
                # Input '1' (Group for centering) and '0' (Group for output - System) are piped
                echo "    2. Running gmx trjconv (center/whole)..."
                # Using printf to ensure both inputs are piped correctly: '1' for centering, '0' for output
                printf "1\n0\n" | gmx trjconv -s step5_production.tpr -f traj_nojump.xtc -o traj_whole_centered.xtc -n new_index.ndx -pbc mol -ur compact -center
                
                if [ $? -ne 0 ]; then
                    echo "    !! Error running second gmx trjconv. Continuing to next pose."
                else
                    echo "    ✅ New trajectory traj_whole_centered.xtc created."
                fi

                # Change back to the base directory
                cd "$BASE_DIR"
            fi
        fi
    done
    echo "---"
done

echo "Script finished. All target directories processed."

# Next recommended step: You can now use traj_whole_centered.xtc for MM/PBSA
echo "The file 'traj_whole_centered.xtc' is ready to be used as your trajectory input for gmx_MMPBSA."
