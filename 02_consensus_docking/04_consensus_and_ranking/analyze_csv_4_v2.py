import pandas as pd
import numpy as np
from collections import defaultdict

# --- Configuration ---
FILE_LIST = [
    '/home/kumara/malakafahs/dock/b2_ad4_final_results_2.csv',
    '/home/kumara/malakafahs/dock/b2_dock6_final_results_2.csv',
    '/home/kumara/malakafahs/dock/b2_vina_final_results_2.csv',
    '/home/kumara/malakafahs/dock/b2_vinardo_final_results_2.csv',
]
# Columns needed for the pose identification and filtering
REQUIRED_COLS = ['Full_Name', 'Binding_Energy', 'com_distance', 'rmsd']
TOP_N_CUTOFF = 500  # A reasonable cutoff to analyze the top candidates for consensus
COM_DISTANCE_CUTOFF = 10.0  # Å
RMSD_CUTOFF = 3.0  # Å

def get_ranked_molecules_efficiently(file_path: str) -> pd.DataFrame:
    """
    Reads a large docking results file, filters poses by COM distance and RMSD,
    finds the best valid pose per molecule, and returns the full ranked list, 
    including the Full_Name (pose ID) of the best pose.
    """
    algorithm_name = file_path.split('_')[1]

    print(f"-> Processing {file_path} (Algorithm: {algorithm_name})...")
    
    # 1. Memory-efficient read (Only load essential columns)
    try:
        df = pd.read_csv(file_path, usecols=REQUIRED_COLS)
    except Exception as e:
        print(f"FATAL ERROR: Could not read {file_path} or columns are missing. Details: {e}")
        return pd.DataFrame()

    # 2. Extract the unique molecule identifier (ZINCxxxx) from Full_Name
    df['Molecule_ID'] = df['Full_Name'].str.split('_').str[0]

    # 3a. Apply COM distance filter
    df_filtered = df[df['com_distance'] <= COM_DISTANCE_CUTOFF].copy()
    # 3b. Apply RMSD filter
    df_filtered = df_filtered[df_filtered['rmsd'] <= RMSD_CUTOFF].copy()
    
    if df_filtered.empty:
        print(f"    WARNING: No poses for {algorithm_name} met COM <= {COM_DISTANCE_CUTOFF} Å and RMSD <= {RMSD_CUTOFF} Å.")
        return pd.DataFrame()

    # 4. Find the best valid pose (lowest Binding_Energy) for each unique molecule
    df_unique_best_pose = (
        df_filtered.sort_values(by='Binding_Energy', ascending=True)
                   .drop_duplicates(subset=['Molecule_ID'], keep='first')
    )
    
    # 5. Apply final ranking and create algorithm-specific columns
    df_ranked = df_unique_best_pose.reset_index(drop=True).copy()
    
    # Store these full rank reports separately
    df_full_report = df_ranked[['Molecule_ID', 'Full_Name', 'Binding_Energy', 'com_distance', 'rmsd']].copy()
    df_full_report.columns = ['Molecule_ID', 'Best_Pose_Full_Name', 'Best_Binding_Energy', 'Best_COM_Distance', 'Best_RMSD']
    df_full_report['Algorithm_Rank'] = df_full_report.index + 1
    
    # Save the individual report
    individual_output_file = f'{algorithm_name}_full_ranked_report.csv'
    df_full_report.to_csv(individual_output_file, index=False)
    print(f"    Saved full ranked report to: {individual_output_file}")

    # Create columns needed for the consensus merge
    df_ranked[f'{algorithm_name}_Rank'] = df_ranked.index + 1
    df_ranked[f'{algorithm_name}_Energy'] = df_ranked['Binding_Energy']
    df_ranked[f'{algorithm_name}_COM_Dist'] = df_ranked['com_distance']
    df_ranked[f'{algorithm_name}_RMSD'] = df_ranked['rmsd']
    df_ranked[f'{algorithm_name}_Best_Pose'] = df_ranked['Full_Name']
    df_ranked['Algorithm'] = algorithm_name 

    # 6. Return only the essential columns for consensus analysis (Top N)
    print(f"    Done. {len(df_ranked)} unique molecules ranked.")
    
    return df_ranked[['Molecule_ID', 'Algorithm', f'{algorithm_name}_Rank', f'{algorithm_name}_Energy',
                      f'{algorithm_name}_COM_Dist', f'{algorithm_name}_RMSD', f'{algorithm_name}_Best_Pose']].head(TOP_N_CUTOFF)


# --- Main Execution ---

# 1. Process all files and store the top N results
all_ranked_dfs = []
for file in FILE_LIST:
    df_top_n = get_ranked_molecules_efficiently(file)
    if not df_top_n.empty:
        all_ranked_dfs.append(df_top_n)

if len(all_ranked_dfs) == 0:
    print("\nNo data could be processed successfully. Exiting.")
else:
    # 2. Consolidate all top N lists into a single consensus DataFrame
    df_consensus = all_ranked_dfs[0].drop(columns=['Algorithm'])
    
    for i in range(1, len(all_ranked_dfs)):
        df_consensus = pd.merge(
            df_consensus, 
            all_ranked_dfs[i].drop(columns=['Algorithm']), 
            on='Molecule_ID', 
            how='outer'
        )

    # 3. Determine Consensus Level for each molecule
    rank_cols = [col for col in df_consensus.columns if '_Rank' in col]
    num_algorithms = len(rank_cols)
    
    df_consensus['Consensus_Count'] = df_consensus[rank_cols].notna().sum(axis=1)

    # 4. Include molecules with at least 2-way consensus
    df_consensus_report = df_consensus[df_consensus['Consensus_Count'] >= 2].copy()
    
    # Calculate Average Rank for all remaining molecules
    df_consensus_report['Average_Rank'] = df_consensus_report[rank_cols].mean(axis=1)
    
    # Sort first by Consensus Count (highest first), then by best Average Rank
    df_final_report = df_consensus_report.sort_values(
        by=['Consensus_Count', 'Average_Rank'], 
        ascending=[False, True]
    )

    print("\n--- Consensus Analysis Summary ---")
    print(f"4-Way Consensus: {len(df_final_report[df_final_report['Consensus_Count'] == 4])} molecules.")
    print(f"3-Way Consensus: {len(df_final_report[df_final_report['Consensus_Count'] == 3])} molecules.")
    print(f"2-Way Consensus: {len(df_final_report[df_final_report['Consensus_Count'] == 2])} molecules.")
    
    # 5. Generate and Save Final Consolidated CSV Report
    column_order_base = ['Molecule_ID', 'Consensus_Count', 'Average_Rank']
    available_algorithms = [col.split('_')[0] for col in rank_cols]
    
    # Order for the report: Rank, Energy, COM_Dist, RMSD, Best_Pose (Full_Name)
    metrics = ['Rank', 'Energy', 'COM_Dist', 'RMSD', 'Best_Pose']
    algo_specific_cols = [f'{alg}_{metric}' for alg in available_algorithms for metric in metrics]
    
    final_column_order = column_order_base + algo_specific_cols
    
    # Fill NaN values for cleaner output
    fill_values = {col: np.nan for col in algo_specific_cols if '_Pose' not in col}
    fill_values.update({col: 'N/A' for col in algo_specific_cols if '_Pose' in col})
    fill_values['Average_Rank'] = np.nan
    
    df_final_report_csv = df_final_report[final_column_order].fillna(fill_values)

    # Save the consolidated CSV
    output_csv_file = f'sirtuin_consolidated_consensus_report_top{TOP_N_CUTOFF}.csv'
    df_final_report_csv.to_csv(output_csv_file, index=False)
    
    print(f"\n✅ Consolidated Consensus Report (2, 3, and 4-Way) saved to: {output_csv_file}")
    print("The individual algorithm rank reports were saved in the processing step above.")
