#!/bin/bash

# Define the new value for nsteps
NEW_NSTEPS="10000000"

# Define the target filename
TARGET_FILE="step5_production.mdp"

# 1. Use 'find' to locate all the target files.
#    This pattern is specific to: ./ZINC*/<4 sub dirs>/*/step5_production.mdp
find . -path "./ZINC*/*/gromacs/$TARGET_FILE" -type f -print0 | while IFS= read -r -d $'\0' filepath; do
    
    echo "Processing: $filepath"
    
    # 2. Use 'sed' to perform an in-place replacement.
    #    The regex ensures we find the 'nsteps' line regardless of existing spacing.
    #    The replacement string "nsteps = ${NEW_NSTEPS}" enforces the desired
    #    single space around the equals sign for this specific line.
    
    # Pattern to match: ^nsteps followed by any space, '=', any space, and the number
    # Replacement string: nsteps <single space> = <single space> 10000000
    sed -i "s/^nsteps[[:space:]]*=[[:space:]]*[0-9]*/nsteps = ${NEW_NSTEPS}/" "$filepath"
    
done

echo "---"
echo "✅ All 'nsteps' values in step5_production.mdp files under ZINC*/.../gromacs/ have been successfully updated to the standard format: nsteps = $NEW_NSTEPS."
