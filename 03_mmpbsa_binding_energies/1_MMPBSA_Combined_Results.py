import os
import re
import csv
from collections import OrderedDict, defaultdict

# --- Configuration ---
# The base directory where the ZINC... and reference folders are located
BASE_DIR = os.path.expanduser("~/malakafahs/gmx/b2_ligs")
OUTPUT_FILE = "MMPBSA_Combined_Results.csv"
POSE_ORDER = ["ad4", "vina", "vinardo", "dock6"]


def parse_final_results(filepath, prefix=""):
    """
    Parses FINAL_RESULTS_MMPBSA.dat to extract all energy components (Average)
    and the final TOTAL.
    """
    data = OrderedDict()
    # Regex to capture the energy component name and its Average value
    # Assumes Average is the second column (after Energy Component)
    # Allows for optional internal spaces in component names (e.g., '1-4 VDW')
    pattern = re.compile(r"^(Δ[\w\s-]+)\s+(-?\d+\.\d+)\s+")

    try:
        with open(filepath, 'r') as f:
            content = f.read()

        # Extract TOTAL first
        total_match = re.search(
            r"^ΔTOTAL\s+(-?\d+\.\d+)", content, re.MULTILINE
        )
        if total_match:
            data[f"{prefix}ΔTOTAL"] = float(total_match.group(1))

        # Extract all other components
        for line in content.splitlines():
            match = pattern.match(line)
            if match and match.group(1) != "ΔTOTAL":
                key = match.group(1).replace(" ", "")  # Clean up key name
                value = float(match.group(2))
                data[f"{prefix}{key}"] = value

    except FileNotFoundError:
        # This is the expected warning if the file is missing
        print(f"Warning: File not found at {filepath}")
    except Exception as e:
        print(f"Error parsing {filepath}: {e}")

    return data


def parse_final_decomp(filepath, prefix=""):
    """
    Parses FINAL_DECOMP_MMPBSA.dat to extract per-residue decomposition data.
    It focuses on the TOTAL (Average) value for each residue.
    The TOTAL Average is typically the 17th column (index 16) when splitting by whitespace.
    """
    residue_data = OrderedDict()
    # The TOTAL (Average) column index based on the standard 19-column output
    TOTAL_AVG_INDEX = 16 

    try:
        with open(filepath, 'r') as f:
            lines = f.readlines()

        # Find the line that starts the residue data (e.g., 'R:A:...')
        start_index = -1
        for i, line in enumerate(lines):
            # Check for lines starting with 'R:' which indicates residue data
            if line.strip().startswith('R:'):
                start_index = i
                break

        if start_index == -1:
            # print(f"Warning: Could not find residue data start in {filepath}")
            return residue_data

        # Process residue lines from the start_index
        for line in lines[start_index:]:
            line = line.strip()
            if not line.startswith('R:'):
                continue
                
            # FIX: Replace commas with spaces to handle comma-separated fields
            # and ensure split() works correctly across mixed delimiters.
            line_cleaned = line.replace(',', ' ') 

            # Use split() without arguments to handle variable whitespace
            parts = line_cleaned.split() 
            
            # Data lines must have enough parts to reach the TOTAL Average
            if len(parts) < TOTAL_AVG_INDEX + 1:
                continue

            # The Residue ID is the first part (e.g., R:A:LEU:202)
            # We strip trailing comma or other non-numeric chars from the ID
            residue_id = parts[0].strip('>,') 

            try:
                # Extract the TOTAL Average value (index 16)
                total_avg = float(parts[TOTAL_AVG_INDEX].strip('>,'))
                
                # Format the key clearly: RES_R:A:THR:260_TOTAL_AVG
                key = f"{prefix}RES_{residue_id}_TOTAL_AVG"
                residue_data[key] = total_avg
            
            except ValueError:
                # Skip lines where the target part is not a valid float 
                continue

    except FileNotFoundError:
        print(f"Warning: File not found at {filepath}")
    except Exception as e:
        print(f"Error parsing {filepath}: {e}")

    return residue_data


def get_pose_tag(pose_dir_name):
    """
    Determines the pose tag (ad4, vina, vinardo, dock6) from the pose directory name.
    It returns the first matching tag.
    """
    for tag in POSE_ORDER:
        if tag in pose_dir_name.lower():
            return tag
    return "unknown"


def process_zinc_directory(zinc_dir, full_path):
    """
    Processes a single ZINC... directory, collecting data from its 4 pose subdirectories.
    """
    molecule_name = os.path.basename(full_path)
    all_pose_data = []

    # Iterate through subdirectories (the 4 docking poses)
    for pose_dir_name in sorted(os.listdir(full_path)):
        pose_path = os.path.join(full_path, pose_dir_name)
        if not os.path.isdir(pose_path):
            continue

        # Check if it contains a 'gromacs' folder
        gromacs_path = os.path.join(pose_path, "gromacs")
        if not os.path.isdir(gromacs_path):
            continue

        # Determine the pose tag
        pose_tag = get_pose_tag(pose_dir_name)

        # --- Data Structure for this Pose ---
        pose_data = OrderedDict()
        pose_data["Molecule_Name"] = molecule_name
        pose_data["Pose_Directory"] = pose_dir_name
        pose_data["Pose_Tag"] = pose_tag
        pose_data["Type"] = "ZINC_LIGAND" # Mark for ZINC ligands

        # --- 1. ZINC Ligand (Main) Data ---
        lig_results_file = os.path.join(gromacs_path, "FINAL_RESULTS_MMPBSA.dat")
        lig_decomp_file = os.path.join(gromacs_path, "FINAL_DECOMP_MMPBSA.dat")

        # Parse main results (no prefix needed for main ligand)
        pose_data.update(parse_final_results(lig_results_file))

        # Parse main decomposition (RESIDUE DATA ADDED HERE)
        pose_data.update(parse_final_decomp(lig_decomp_file))

        # --- 2. NAD Ligand (Secondary) Data ---
        nad_gromacs_path = os.path.join(gromacs_path, "NAD")
        if os.path.isdir(nad_gromacs_path):
            nad_results_file = os.path.join(nad_gromacs_path, "FINAL_RESULTS_MMPBSA.dat")
            nad_decomp_file = os.path.join(nad_gromacs_path, "FINAL_DECOMP_MMPBSA.dat")

            # Parse NAD results (with 'NAD_' prefix)
            pose_data.update(parse_final_results(nad_results_file, prefix="NAD_"))

            # Parse NAD decomposition (with 'NAD_' prefix) (RESIDUE DATA ADDED HERE)
            pose_data.update(parse_final_decomp(nad_decomp_file, prefix="NAD_"))

        all_pose_data.append(pose_data)

    return all_pose_data


def process_reference_directory(ref_dir, full_path, is_nad_only=False):
    """
    Processes the reference directories (4zzj_nad_p53_4tq/ and 4zzj_nad_p53/).
    """
    molecule_name = os.path.basename(full_path)
    ref_data = OrderedDict()
    ref_data["Molecule_Name"] = molecule_name
    ref_data["Pose_Directory"] = "N/A"
    ref_data["Pose_Tag"] = "Reference"
    ref_data["Type"] = "REFERENCE"

    gromacs_path = os.path.join(full_path, "gromacs")

    if not os.path.isdir(gromacs_path):
        print(f"Warning: 'gromacs' not found in reference directory {molecule_name}")
        return []

    if not is_nad_only:
        # --- Reference Ligand (4TQ) Data - Only for 4zzj_nad_p53_4tq ---
        lig_results_file = os.path.join(gromacs_path, "FINAL_RESULTS_MMPBSA.dat")
        lig_decomp_file = os.path.join(gromacs_path, "FINAL_DECOMP_MMPBSA.dat")

        # Parse main results (no prefix needed for reference ligand)
        ref_data.update(parse_final_results(lig_results_file))

        # Parse main decomposition (RESIDUE DATA ADDED HERE)
        ref_data.update(parse_final_decomp(lig_decomp_file))

    # --- NAD Ligand Data (Secondary) ---
    nad_gromacs_path = os.path.join(gromacs_path, "NAD")
    if os.path.isdir(nad_gromacs_path):
        nad_results_file = os.path.join(nad_gromacs_path, "FINAL_RESULTS_MMPBSA.dat")
        nad_decomp_file = os.path.join(nad_gromacs_path, "FINAL_DECOMP_MMPBSA.dat")

        # Prefix is 'NAD_' even if it's the only data (for consistency)
        ref_data.update(parse_final_results(nad_results_file, prefix="NAD_"))
        ref_data.update(parse_final_decomp(nad_decomp_file, prefix="NAD_"))

    # If the reference is NAD-only (4zzj_nad_p53), we need to ensure the NAD data is present
    if is_nad_only and not any(k.startswith("NAD_") for k in ref_data.keys()):
        print(f"Warning: NAD data not found in {molecule_name}")
        return []

    return [ref_data]


def main():
    """
    Main function to orchestrate the directory traversal and data collection.
    """
    all_data_rows = []

    # 1. Get all relevant directories
    try:
        top_level_dirs = sorted([
            d for d in os.listdir(BASE_DIR)
            if os.path.isdir(os.path.join(BASE_DIR, d)) and d != 'others'
        ])
    except FileNotFoundError:
        print(f"Error: Base directory not found at {BASE_DIR}. Please check the path.")
        return

    # 2. Process ZINC ligand directories
    print("Processing ZINC ligand directories...")
    zinc_dirs = [d for d in top_level_dirs if d.startswith("ZINC")]
    for zinc_dir in zinc_dirs:
        full_path = os.path.join(BASE_DIR, zinc_dir)
        all_data_rows.extend(process_zinc_directory(zinc_dir, full_path))
        print(f"  > Processed {zinc_dir} ({len(all_data_rows)} rows collected so far)")

    # 3. Process Reference directories
    print("\nProcessing Reference directories...")
    ref_dirs = [d for d in top_level_dirs if d.startswith("4zzj")]

    # 4zzj_nad_p53_4tq (Reference with 4TQ ligand and NAD data)
    if "4zzj_nad_p53_4tq" in ref_dirs:
        path = os.path.join(BASE_DIR, "4zzj_nad_p53_4tq")
        all_data_rows.extend(process_reference_directory("4zzj_nad_p53_4tq", path, is_nad_only=False))
        print("  > Processed 4zzj_nad_p53_4tq (4TQ and NAD)")

    # 4zzj_nad_p53 (Reference with only NAD data)
    if "4zzj_nad_p53" in ref_dirs:
        path = os.path.join(BASE_DIR, "4zzj_nad_p53")
        all_data_rows.extend(process_reference_directory("4zzj_nad_p53", path, is_nad_only=True))
        print("  > Processed 4zzj_nad_p53 (NAD only)")

    # 4. Final compilation and CSV writing
    if not all_data_rows:
        print("\nNo data rows were collected. Check directory structure and file names.")
        return

    # Collect all possible column headers from all rows
    print("\nPreparing to write CSV...")
    fieldnames = set()
    for row in all_data_rows:
        fieldnames.update(row.keys())

    # Ensure fixed columns are first, then sort the rest alphabetically
    fixed_headers = ["Molecule_Name", "Pose_Directory", "Pose_Tag", "Type"]
    
    # Separate the MMPBSA and Decomposition headers for better control
    mmpbsa_headers = sorted([h for h in fieldnames if h.startswith("Δ") or h.startswith("NAD_Δ")])
    decomp_headers = sorted([h for h in fieldnames if h.startswith("RES_") or h.startswith("NAD_RES_")])
    
    # Remaining headers (should only be the fixed ones, but kept for robustness)
    other_headers = sorted(list(fieldnames - set(fixed_headers) - set(mmpbsa_headers) - set(decomp_headers)))

    # Construct the final sorted fieldnames list
    # Decomposition headers are correctly placed at the end of the file
    sorted_fieldnames = fixed_headers + mmpbsa_headers + other_headers + decomp_headers

    # Write the data to the CSV file
    try:
        with open(OUTPUT_FILE, 'w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=sorted_fieldnames)
            writer.writeheader()
            writer.writerows(all_data_rows)
        print(f"✅ Success! Data collected into **{OUTPUT_FILE}** with {len(all_data_rows)} rows.")
    except Exception as e:
        print(f"An error occurred while writing the CSV file: {e}")

if __name__ == "__main__":
    main()
