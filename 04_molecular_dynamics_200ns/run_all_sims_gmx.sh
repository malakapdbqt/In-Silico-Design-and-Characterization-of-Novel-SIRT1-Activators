#!/bin/bash

# The name of the script to execute in each directory
SCRIPT_TO_RUN="run_all.sh"

# The target directory name (the final destination folder)
TARGET_DIR="gromacs"

echo "Starting sequential execution of $SCRIPT_TO_RUN in all target directories..."
echo "Target pattern: ZINC*/<4 sub dirs>/gromacs/"
echo "---"

# 1. Use 'find' to locate all the target 'gromacs' directories.
#    -path "./ZINC*/*/${TARGET_DIR}": Matches the full path structure: 
#      CurrentDir/ZINC_DIR/SUB_DIR/gromacs
#    -type d: Ensures we are matching only directories.
#    -print0: Safely handles directory names with spaces.

find . -path "./ZINC*/*/${TARGET_DIR}" -type d -print0 | while IFS= read -r -d $'\0' dest_dir; do
    
    echo "▶️ Entering directory: $dest_dir"
    
    # Check if the script exists and is executable before running
    if [ -x "$dest_dir/$SCRIPT_TO_RUN" ]; then
        
        # Change into the target directory
        cd "$dest_dir" || { echo "ERROR: Could not change directory to $dest_dir. Skipping."; continue; }
        
        # Execute the script
        echo "   Executing $SCRIPT_TO_RUN..."
        
        # We use ./$SCRIPT_TO_RUN to ensure the script is run from the current directory
        # The '&>' redirects both standard output and standard error to the log file
        ./"$SCRIPT_TO_RUN" &> run_all.log
        
        # Check the exit status of the script
        if [ $? -eq 0 ]; then
            echo "   ✅ $SCRIPT_TO_RUN completed successfully."
        else
            echo "   ❌ $SCRIPT_TO_RUN failed! Check $dest_dir/run_all.log for details."
        fi
        
        # Change back to the starting directory to continue the loop
        cd - > /dev/null
        
    else
        echo "   ⚠️ WARNING: $SCRIPT_TO_RUN not found or not executable in $dest_dir. Skipping."
    fi
    
    echo "---"

done

echo "🎉 Sequential execution complete for all targeted directories."
