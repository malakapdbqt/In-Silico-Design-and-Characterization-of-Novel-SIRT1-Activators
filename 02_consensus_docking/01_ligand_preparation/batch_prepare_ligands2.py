#!/usr/bin/env python
"""
batch_prepare_ligands.py

Iterates through all PDBQT files in an input directory, runs 
prepare_ligand4.py on each by changing the current working directory (CWD) 
to match the tool's requirement, and saves the prepared output to a 
separate output directory with the same base name.

NOTE: This script is fully compatible with Python 2.7.
"""

import os
import sys
import glob
import subprocess
from tqdm import tqdm

# --- Configuration ---
PREPARE_TOOL = 'prepare_ligand4.py' # Name of the AutoDock Vina tool
INPUT_EXTENSION = '.pdbqt'
OUTPUT_EXTENSION = '.pdbqt' 
# IMPORTANT: This script ASSUMES prepare_ligand4.py supports an '-o' flag 
# for specifying a separate output file/directory.

def prepare_ligands_batch(input_dir, output_dir):
    """
    Finds PDBQT files in input_dir and runs prepare_ligand4.py on them.
    Changes CWD to input_dir before processing to satisfy the tool's path requirements.
    """
    
    # Check if the preparation tool is accessible (using Popen for Python 2.7 compatibility)
    try:
        p = subprocess.Popen([PREPARE_TOOL, '-h'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        p.communicate() 
    except OSError: # FileNotFoundError is OSError in Python 2.7
        print "ERROR: The preparation tool '%s' was not found." % PREPARE_TOOL
        print "Please ensure it is in your system PATH or use its full path."
        return

    # 1. Setup paths and check input
    input_dir = os.path.abspath(input_dir)
    output_dir = os.path.abspath(output_dir)

    if not os.path.isdir(input_dir):
        print "ERROR: Input directory not found: %s" % input_dir
        sys.exit(1)

    # os.makedirs(..., exist_ok=True) is Python 3+. Use try/except in Python 2.
    try:
        os.makedirs(output_dir)
    except OSError:
        if not os.path.isdir(output_dir):
            raise
    
    # 2. Find all PDBQT files
    search_path = os.path.join(input_dir, '*%s' % INPUT_EXTENSION)
    ligand_files = glob.glob(search_path)
    
    if not ligand_files:
        print "WARNING: No %s files found in %s" % (INPUT_EXTENSION, input_dir)
        return

    print "Found %d ligands to prepare." % len(ligand_files)
    
    successful_preparations = 0
    errors = []
    
    # 3. Handle directory changes
    original_cwd = os.getcwd()
    
    try:
        # CRITICAL FIX: Change CWD to the input directory
        os.chdir(input_dir)
        print "\nTemporarily changed directory to: %s" % os.getcwd()

        # 4. Process each file
        for input_file_path in tqdm(ligand_files, desc="Preparing Ligands"):
            
            # Get the filename (e.g., ZINC30009884.pdbqt)
            input_filename = os.path.basename(input_file_path)
            
            # Construct the absolute output path in the target directory
            output_file_path = os.path.join(output_dir, input_filename)

            # Check if output already exists (to skip and save time)
            if os.path.exists(output_file_path):
                successful_preparations += 1
                continue

            # 5. Construct and run the command (using Popen for Python 2.7 compatibility)
            command = [
                PREPARE_TOOL, 
                '-l', input_filename, 
                '-o', output_file_path
            ]
            
            try:
                # Use Popen to run the command and capture output
                p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                stdout, stderr = p.communicate()
                
                if os.path.exists(output_file_path):
                    successful_preparations += 1
                else:
                    # Command ran, but output file wasn't created
                    error_msg = stderr.strip() if stderr else 'No error message captured.'
                    errors.append((input_filename, "Output file not created. Error message:\n" + error_msg[:200] + "..."))
                
            except Exception as e:
                errors.append((input_filename, "Script execution failed: %s" % e))

    finally:
        # 6. Restore the original CWD, even if errors occurred
        os.chdir(original_cwd)
        print "\nRestored directory to: %s" % os.getcwd()


    # 7. Final Report
    print "\n" + "="*50
    print "BATCH LIGAND PREPARATION SUMMARY"
    print "="*50
    print "Total files found: %d" % len(ligand_files)
    print "Successfully prepared: %d" % successful_preparations
    print "Files with errors: %d" % len(errors)

    if errors:
        print "\n--- First 5 Errors ---"
        for filename, msg in errors[:5]:
            print "[%s]: %s" % (filename, msg.replace('\n', ' '))

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print "Usage: prepare_ligands.py <input_dir> <output_dir>"
        print "\nExample:"
        print "prepare_ligands.py ligands_raw/ ligands_prepared/"
        sys.exit(1)
        
    input_directory = sys.argv[1]
    output_directory = sys.argv[2]
    
    prepare_ligands_batch(input_directory, output_directory)
