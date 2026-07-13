#!/bin/bash

# --- Configuration ---

# List of files to be copied (must exist in the current directory)
FILES_TO_COPY=(
    "mmpbsa.in"
    "step5_production.mdp"
    "run_all.sh"
    "run.sh"
)

# The target directory name (the final destination folder)
TARGET_DIR="gromacs"

# --- Script Logic ---

echo "Starting file copy process..."
echo "Source files: ${FILES_TO_COPY[@]}"
echo "Target pattern: ZINC*/<4 sub dirs>/gromacs/"
echo "---"

# 1. Use 'find' to locate all the target 'gromacs' directories.
#    -path "./ZINC*/*/gromacs": Matches the full path structure: 
#      CurrentDir/ZINC_DIR/SUB_DIR/gromacs
#    -type d: Ensures we are matching only directories.
#    -print0: Safely handles directory names with spaces.

find . -path "./ZINC*/*/${TARGET_DIR}" -type d -print0 | while IFS= read -r -d $'\0' dest_dir; do
    
    echo "Copying to: $dest_dir"
    
    # 2. Loop through each file in the list and copy it.
    for file in "${FILES_TO_COPY[@]}"; do
        
        # Check if the source file exists before trying to copy
        if [ -f "$file" ]; then
            
            # 'cp -f' copies and forces the overwrite if the file already exists.
            cp -f "$file" "$dest_dir/"
            
        else
            echo "⚠️ WARNING: Source file '$file' not found in current directory. Skipping."
        fi
        
    done
    
done

echo "---"
echo "✅ File copy complete. All four files have been copied and overwritten in all targeted 'gromacs' subdirectories."

# Optional: Set executable permissions for the shell scripts
echo "Setting executable permissions for run_all.sh and run.sh in all target directories..."

find . -path "./ZINC*/*/${TARGET_DIR}/run*.sh" -type f -exec chmod +x {} \;

echo "✅ Executable permissions updated."
