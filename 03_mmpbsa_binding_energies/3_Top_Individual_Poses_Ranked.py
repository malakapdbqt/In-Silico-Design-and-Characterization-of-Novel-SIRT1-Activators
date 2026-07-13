import pandas as pd
import io
import numpy as np

# --- Configuration ---
INPUT_FILE = "MMPBSA_Analysis_Ready.csv"
OUTPUT_FILE_BEST_POSE = "Top_Individual_Poses_Ranked.csv"

def find_best_individual_pose(input_file, output_file):
    """
    Calculates ΔΔG_NAD for every pose and ranks them to find the single best pose.
    """
    try:
        df = pd.read_csv(input_file)
    except FileNotFoundError:
        print(f"Error: Input file '{input_file}' not found.")
        return

    # --- 1. Define Baseline and Calculate Activation Metric (ΔΔG_NAD) ---

    # Baseline NAD binding energy (ΔG_NAD, apo) from 4zzj_nad_p53
    ref_nad_row = df[df['Molecule_Name'] == '4zzj_nad_p53']
    if ref_nad_row.empty or ref_nad_row['NAD_ΔTOTAL'].isnull().all():
        print("\nERROR: Could not find valid NAD binding reference (4zzj_nad_p53).")
        return
    
    # Use the specific value from the baseline reference
    delta_g_nad_apo = ref_nad_row['NAD_ΔTOTAL'].iloc[0]
    print(f"Baseline ΔG_NAD, apo: {delta_g_nad_apo:.2f} kcal/mol")
    
    # Calculate ΔΔG_NAD for all poses
    df['ΔΔG_NAD (Activation Efficacy)'] = df['NAD_ΔTOTAL'] - delta_g_nad_apo

    # --- 2. Filter, Clean, and Rank Poses ---
    
    # Filter out the NAD-only reference (since it has no allosteric data)
    best_pose_df = df[df['Molecule_Name'] != '4zzj_nad_p53'].copy()
    
    # Rename columns for clarity
    best_pose_df.rename(columns={'ΔTOTAL': 'ΔG_Allost_Binding (kcal/mol)', 
                                 'NAD_ΔTOTAL': 'ΔG_NAD_Holo (kcal/mol)'}, 
                        inplace=True)
    
    # Sort: Prioritize by the most favorable (most negative) activation metric
    best_pose_df.sort_values(by='ΔΔG_NAD (Activation Efficacy)', ascending=True, inplace=True)

    # --- 3. Report the Top Poses ---
    
    # Select key columns for the final report
    report_cols = ['Molecule_Name', 'Pose_Tag', 'ΔG_Allost_Binding (kcal/mol)', 
                   'ΔG_NAD_Holo (kcal/mol)', 'ΔΔG_NAD (Activation Efficacy)']
    
    final_report_df = best_pose_df[report_cols]
    
    # Save the report
    final_report_df.to_csv(output_file, index=False, float_format='%.2f')
    
    print(f"\n✅ Top poses ranked successfully and saved to: **{output_file}**")
    print("\n--- Top 5 Individual Poses (Highest Activation Efficacy) ---")
    print(final_report_df.head(5).to_markdown(index=False, floatfmt=".2f"))
    
    # --- 4. Identify the BEST Pose overall ---
    best_overall_pose = final_report_df.iloc[0]
    
    print("\n------------------------------------------------------------")
    print("🥇 The Single Best Pose for Allosteric Activation:")
    print(f"Molecule: {best_overall_pose['Molecule_Name']}")
    print(f"Pose: {best_overall_pose['Pose_Tag']}")
    print(f"ΔΔG_NAD (Efficacy): {best_overall_pose['ΔΔG_NAD (Activation Efficacy)']:.2f} kcal/mol")
    print(f"ΔG_Allost (Affinity): {best_overall_pose['ΔG_Allost_Binding (kcal/mol)']:.2f} kcal/mol")
    print("------------------------------------------------------------")

if __name__ == "__main__":
    find_best_individual_pose(INPUT_FILE, OUTPUT_FILE_BEST_POSE)
