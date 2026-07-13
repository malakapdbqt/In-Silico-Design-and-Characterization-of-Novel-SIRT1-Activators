import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
import re
import csv
from operator import itemgetter
from collections import OrderedDict
from matplotlib.gridspec import GridSpec

# --- Configuration ---
BASE_DIR = os.path.expanduser("~/malakafahs/gmx/b2_ligs")
RAW_INPUT_FILE = "MMPBSA_Combined_Results.csv"
ANALYSIS_OUTPUT_FILE = "MMPBSA_Analysis_Ready.csv"

ALLOST_N = 10 
NAD_N = 15
POSES_PER_PAGE = 7  
MAX_Y_LIMIT = 0.1

RESIDUE_COLOR_MAP = {}
RESIDUE_COLORS = list(plt.cm.tab20.colors) + list(plt.cm.tab20b.colors) + list(plt.cm.tab20c.colors)

# NEW: Mapping for shortened pose tags to full descriptions
# You can customize this dictionary based on the actual short tags in your CSV
# and their corresponding full names.
POSE_TAG_MAPPING = {
    "vina": "Autodock Vina",
    "vinardo": "Vinardo",
    "smina": "smina (SMINA)",
    # Add more mappings as needed for other docking programs
}


# --- CSV GENERATION LOGIC ---

def parse_residue_tag(residue_key):
    try:
        core_id = residue_key.replace("NAD_RES_", "").replace("RES_", "").replace("_TOTAL_AVG", "")
        parts = core_id.split(':')
        if len(parts) >= 4:
            return f"{parts[1]}:{parts[2]}:{parts[3]}"
        return core_id
    except:
        return ""

def generate_analysis_csv():
    if not os.path.exists(RAW_INPUT_FILE):
        print(f"Error: Raw input file '{RAW_INPUT_FILE}' not found.")
        return False
    
    try:
        with open(RAW_INPUT_FILE, mode='r', newline='') as csvfile:
            reader = csv.DictReader(csvfile)
            data_rows = list(reader)
    except Exception as e:
        print(f"Error reading raw CSV {RAW_INPUT_FILE}: {e}")
        return False

    if not data_rows:
        print("Raw CSV is empty.")
        return False

    all_analysis_rows = []
    output_headers = ["Molecule_Name", "Pose_Tag", "ΔTOTAL", "NAD_ΔTOTAL"]
    for i in range(1, ALLOST_N + 1):
        output_headers.append(f"ALLOST_RES_E_{i:02d}")
        output_headers.append(f"ALLOST_RES_T_{i:02d}")
    for i in range(1, NAD_N + 1):
        output_headers.append(f"NAD_RES_E_{i:02d}")
        output_headers.append(f"NAD_RES_T_{i:02d}")

    for row in data_rows:
        new_row = OrderedDict()
        new_row["Molecule_Name"] = row.get("Molecule_Name", "")
        new_row["Pose_Tag"] = row.get("Pose_Tag", "")
        try:
            new_row["ΔTOTAL"] = float(row.get("ΔTOTAL"))
        except (ValueError, TypeError):
            new_row["ΔTOTAL"] = np.nan
        try:
            new_row["NAD_ΔTOTAL"] = float(row.get("NAD_ΔTOTAL"))
        except (ValueError, TypeError):
            new_row["NAD_ΔTOTAL"] = np.nan

        def extract_and_rank(prefix, n_count):
            residue_prefix = f"{prefix}RES_"
            residue_energies = []
            for header, value in row.items():
                if header.startswith(residue_prefix):
                    try:
                        energy = float(value)
                        tag = parse_residue_tag(header)
                        if tag:
                            residue_energies.append((tag, energy))
                    except (ValueError, TypeError):
                        continue
            return sorted(residue_energies, key=itemgetter(1))[:n_count]

        allost_ranked = extract_and_rank("", ALLOST_N)
        for i, (tag, energy) in enumerate(allost_ranked):
            new_row[f"ALLOST_RES_E_{i+1:02d}"] = energy
            new_row[f"ALLOST_RES_T_{i+1:02d}"] = tag

        nad_ranked = extract_and_rank("NAD_", NAD_N)
        for i, (tag, energy) in enumerate(nad_ranked):
            new_row[f"NAD_RES_E_{i+1:02d}"] = energy
            new_row[f"NAD_RES_T_{i+1:02d}"] = tag

        all_analysis_rows.append(new_row)

    try:
        with open(ANALYSIS_OUTPUT_FILE, 'w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=output_headers, extrasaction='ignore')
            writer.writeheader()
            writer.writerows(all_analysis_rows)
        print(f"✅ Analysis CSV generated: {ANALYSIS_OUTPUT_FILE}")
        return True
    except Exception as e:
        print(f"Error writing analysis CSV: {e}")
        return False

# --- PLOTTING LOGIC ---

def get_residue_color(residue_id):
    global RESIDUE_COLOR_MAP
    if residue_id not in RESIDUE_COLOR_MAP:
        color_index = len(RESIDUE_COLOR_MAP) % len(RESIDUE_COLORS)
        RESIDUE_COLOR_MAP[residue_id] = RESIDUE_COLORS[color_index]
    return RESIDUE_COLOR_MAP[residue_id]

def get_plot_limits(page_data):
    page_data = page_data.copy()
    page_data[['ΔTOTAL', 'NAD_ΔTOTAL']] = page_data[['ΔTOTAL', 'NAD_ΔTOTAL']].fillna(0)
    min_total_energy = page_data['ΔTOTAL'].min()
    min_nad_total_energy = page_data['NAD_ΔTOTAL'].min()
    allost_res_cols = [c for c in page_data.columns if c.startswith("ALLOST_RES_E_")]
    nad_res_cols = [c for c in page_data.columns if c.startswith("NAD_RES_E_")]
    min_allost_res = page_data[allost_res_cols].min().min() if allost_res_cols else -1.0
    min_nad_res = page_data[nad_res_cols].min().min() if nad_res_cols else -1.0
    abs_min_total = min(min_total_energy, min_nad_total_energy)
    total_limit = np.floor(abs_min_total / 10.0) * 10.0 - 10.0
    if total_limit > -20.0: total_limit = -120.0
    limit_allost = np.floor(min_allost_res / 1.0) * 1.0 - 2.0
    if limit_allost > -3.0: limit_allost = -5.0
    limit_nad = np.floor(min_nad_res / 5.0) * 5.0 - 5.0
    if limit_nad > -5.0: limit_nad = -15.0
    return total_limit, limit_allost, limit_nad

def create_page_plot(page_data, page_num):
    num_poses = len(page_data)
    fig_width = 18 
    fig_height = 4 + num_poses * 6.5
    fig = plt.figure(figsize=(fig_width, fig_height))
    gs = GridSpec(num_poses, 4, figure=fig, width_ratios=[0.2, 1.0, 0.2, 1.0])
    total_limit, allost_res_limit, nad_res_limit = get_plot_limits(page_data)

    pose_axes = []  # collect axes to compute proper title placement

    for plot_idx, row in page_data.iterrows():
        axA = fig.add_subplot(gs[plot_idx, 0])
        axB = fig.add_subplot(gs[plot_idx, 1])
        axC = fig.add_subplot(gs[plot_idx, 2])
        axD = fig.add_subplot(gs[plot_idx, 3])
        pose_axes.append((axA, axD))

        allost_energies = [row[c] for c in row.index if c.startswith("ALLOST_RES_E_") and pd.notna(row[c])]
        allost_tags = [row[c] for c in row.index if c.startswith("ALLOST_RES_T_") and pd.notna(row[c])]
        nad_energies = [row[c] for c in row.index if c.startswith("NAD_RES_E_") and pd.notna(row[c])]
        nad_tags = [row[c] for c in row.index if c.startswith("NAD_RES_T_") and pd.notna(row[c])]

        # --- A --- (Total ΔG Allost)
        gA = row.get("ΔTOTAL")
        if pd.notna(gA):
            axA.bar(["Total"], [gA], color='darkred', alpha=0.9)
            axA.set_ylim(total_limit, MAX_Y_LIMIT)
            axA.axhline(0, color='black', lw=0.5)
            axA.set_title("Total ΔG Allost", fontsize=10, color='darkred')
            axA.text(0, gA, f"{gA:.2f}", ha='center', va='bottom' if gA>0 else 'top', fontsize=8, color='darkred')
        else:
            axA.text(0.5,0.5,"N/A",ha='center',va='center',color='red')
        axA.grid(axis='y', ls='--', alpha=0.5)

        # --- B --- (Allosteric Residues)
        if allost_tags:
            bars_B = axB.bar(allost_tags, allost_energies, color=[get_residue_color(t) for t in allost_tags], alpha=0.8)
            axB.set_ylim(allost_res_limit, MAX_Y_LIMIT)
            axB.axhline(0, color='black', lw=0.5)
            axB.set_title(f"Allosteric Residues (Top {ALLOST_N})", fontsize=10, color='darkgreen')
            axB.tick_params(axis='x', rotation=90, labelsize=9)
            
            # Add values to allost residue bars
            for bar in bars_B:
                yval = bar.get_height()
                axB.text(bar.get_x() + bar.get_width()/2., yval, 
                         f"{yval:.2f}", 
                         ha='center', va='bottom' if yval > 0 else 'top', 
                         fontsize=7, color='black')
        else:
            axB.text(0.5,0.5,"No Allosteric Residue Data",ha='center',va='center',color='red')
        axB.grid(axis='y', ls='--', alpha=0.5)

        # --- C --- (NAD Total ΔG)
        gN = row.get("NAD_ΔTOTAL")
        if pd.notna(gN):
            axC.bar(["Total"], [gN], color='darkmagenta', alpha=0.9)
            axC.set_ylim(total_limit, MAX_Y_LIMIT)
            axC.axhline(0, color='black', lw=0.5)
            axC.set_title("NAD Total ΔG", fontsize=10, color='darkmagenta')
            axC.text(0, gN, f"{gN:.2f}", ha='center', va='bottom' if gN>0 else 'top', fontsize=8, color='darkmagenta')
        else:
            axC.text(0.5,0.5,"N/A",ha='center',va='center',color='red')
        axC.grid(axis='y', ls='--', alpha=0.5)

        # --- D --- (NAD Residues)
        if nad_tags:
            bars_D = axD.bar(nad_tags, nad_energies, color=[get_residue_color(t) for t in nad_tags], alpha=0.8)
            axD.set_ylim(nad_res_limit, MAX_Y_LIMIT)
            axD.axhline(0, color='black', lw=0.5)
            axD.set_title(f"NAD Residues (Top {NAD_N})", fontsize=10, color='darkblue')
            axD.tick_params(axis='x', rotation=90, labelsize=9)
            axD.yaxis.set_label_position("right")
            axD.yaxis.tick_right()
            # Add values to NAD residue bars
            for bar in bars_D:
                yval = bar.get_height()
                axD.text(bar.get_x() + bar.get_width()/2., yval, 
                         f"{yval:.2f}", 
                         ha='center', va='bottom' if yval > 0 else 'top', 
                         fontsize=7, color='black')
        else:
            axD.text(0.5,0.5,"N/A",ha='center',va='center',color='red')
        axD.grid(axis='y', ls='--', alpha=0.5)

    fig.tight_layout(rect=[0, 0.05, 1, 0.98], h_pad=1.5)

    # --- FIXED TITLE POSITION & KEY CHANGE FOR POSE TAG EXPANSION ---
    fig.canvas.draw()  # ensures layout is computed
    for plot_idx, (axA, axD) in enumerate(pose_axes):
        bbox = axA.get_position()
        y_top = bbox.y1 + 0.01
        
        molecule_name = page_data.loc[plot_idx, 'Molecule_Name']
        raw_pose_tag = page_data.loc[plot_idx, 'Pose_Tag']
        
        # KEY CHANGE: Apply mapping to expand the pose tag for the title
        display_pose_tag = POSE_TAG_MAPPING.get(raw_pose_tag.lower(), raw_pose_tag)
        
        title = f"--- {molecule_name} | Pose: {display_pose_tag} ---"
        
        fig.text(0.01, y_top, title, fontsize=12, fontweight='bold', ha='left', va='bottom')

    fig.suptitle(f"MMPBSA Residue Decomposition Comparison (Page {page_num})", fontsize=16, y=0.995)
    plt.savefig(f"MMPBSA_Residue_Comparison_Page_{page_num}.png", dpi=300)
    plt.close(fig)
    print(f"Generated MMPBSA_Residue_Comparison_Page_{page_num}.png")

def plot_mmpbsa_data():
    try:
        df = pd.read_csv(ANALYSIS_OUTPUT_FILE)
    except Exception as e:
        print(f"Error reading analysis CSV: {e}")
        return
    all_tags = set()
    for c in df.columns:
        if c.startswith("ALLOST_RES_T_") or c.startswith("NAD_RES_T_"):
            all_tags.update(df[c].dropna().unique())
    for t in sorted(all_tags):
        get_residue_color(t)
    total_poses = len(df)
    print(f"Total poses found: {total_poses}")
    num_pages = int(np.ceil(total_poses / POSES_PER_PAGE))
    print(f"Generating {num_pages} pages with {POSES_PER_PAGE} poses per page...")
    for i in range(num_pages):
        page_data = df.iloc[i*POSES_PER_PAGE:(i+1)*POSES_PER_PAGE].reset_index(drop=True)
        if not page_data.empty:
            create_page_plot(page_data, i+1)
    print("\n✅ Plot generation complete.")

def main_workflow():
    print("--- Starting MMPBSA Analysis and Plotting Workflow ---")
    if generate_analysis_csv():
        plot_mmpbsa_data()
    else:
        print("Workflow aborted due to CSV generation error.")

if __name__ == "__main__":
    try:
        main_workflow()
    except ModuleNotFoundError as e:
        print(f"Error: {e}. Please ensure pandas, matplotlib, and numpy are installed.")
