# Physics-Based Discovery and Mechanistic Characterization of Novel SIRT1 Allosteric Activators

![Python](https://img.shields.io/badge/Python-3.12%2B-blue?logo=python&logoColor=white)
![GROMACS](https://img.shields.io/badge/GROMACS-2025.1-orange?logo=gromacs&logoColor=white)
![ORCA](https://img.shields.io/badge/ORCA-6.1.0-purple)
![Multiwfn](https://img.shields.io/badge/Multiwfn-3.8-red)
![License](https://img.shields.io/badge/License-MIT-green)

An end-to-end, reproducible, multi-scale computational chemistry and molecular modeling pipeline tracking the discovery of high-efficacy sirtuin 1 (SIRT1) allosteric regulators. This suite integrates high-throughput cheminformatics curation, multi-engine consensus molecular docking, explicit solvent Molecular Dynamics (MD) trajectories, endpoint binding free energy therapeutics ($MM/PBSA$), reaction coordinate kinetics via Umbrella Sampling (US), and quantum electronic pre-organization analysis using Density Functional Theory (DFT).

---

## 🧬 Scientific Overview & Workflow Funnel

Pharmacological activation of Sirtuin 1 (SIRT1) requires precise modulation of an allosteric network linking the cofactor ($NAD^+$) pocket, catalytic residues, and distal ligand-binding sites. This pipeline utilizes structural coordinates derived from the crystallographic template of the human **SIRT1 enzyme complex (PDB ID: 4ZZJ)**. The protocol defines six comparative operational frameworks:
* **System A (`4zzj_nad_p53` / `4zzj_nad_pep`):** Apo-reference baseline containing strictly the cofactor $NAD^+$ and acetylated p53 peptide substrate without an allosteric modifier.
* **System B (`4zzj_nad_p53_4tq` / `4zzj_nad_pep_indol`):** Positive control benchmark containing the native reference allosteric activator (4TQ / Indol derivatives).
* **Systems C–F:** High-confidence discovery activator scaffolds nominated via our multi-scale triage selection rules.

             [ ZINC15 Database (21.5M Compounds) ]
                              │
                              ▼
       Stage 1: Pharmacophore-Guided Screening (5,338 Mols)
                              │
                              ▼
       Stage 2: Rigid ADMET & Complexity Filters (3,557 Mols)
                              │
                              ▼
       Stage 3: Cross-Platform Consensus Docking (13 Leads)
                              │
                              ▼
       Stage 4: MM/PBSA Validation & ΔΔG Binding Triage
                              │
                              ▼
       Stage 5: 200 ns Explicit Solvent Production MD
                              │
                              ▼
       Stage 6: Biased Pulling & Umbrella Sampling PMF
                              │
                              ▼
       Stage 7: DFT Single-Point Electronic Pre-Organization

---

## 📂 Repository Architecture & Directory Mapping

```text
.
├── 01_ligand_library_curation          # Stage 1: HT virtual database curation & filtration
│   ├── rules                           # Computational workflow diagrams and PDF funnels
│   └── logs                            # High-throughput molecular filtration datasets
├── 02_consensus_docking                # Stage 2: Cross-platform multi-engine virtual screening
│   ├── 01_ligand_preparation           # United-atom/all-atom coordinate preprocessing scripts
│   ├── 02_docking_execution_data_arrange # Native execution scripts and receptor map parameters
│   │   ├── ad4                         # AutoDock4-GPU Lamarckian Genetic Algorithm suite
│   │   ├── dock6                       # UCSF DOCK6 gridbox energy minimizations & spheres
│   │   ├── vina                        # AutoDock Vina multi-threaded docking pathways
│   │   └── vinardo                     # Steric-optimized Vinardo execution configurations
│   ├── 03_Generate_Fingerprints        # Structural pose fingerprinting wrappers
│   └── 04_consensus_and_ranking        # Intersecting candidate consensus matrix processors
├── 03_mmpbsa_binding_energies          # Stage 3: Thermodynamic binding free energy triage
├── 04_molecular_dynamics_200ns         # Stage 4: 200 ns production trajectory parameters
├── 05_umbrella_sampling_pmf            # Stage 5: Reaction coordinate potential landscape tools
├── 06_quantum_mechanics_descriptors    # Stage 6: DFT active site electronic pre-organization
└── 07_admet_predictions                # Stage 7: Druggability indices and safety profiles
🏃 Step-by-Step Workflow Execution & Core CodeModule 1: Ligand Library CurationFilters massive chemical space down to drug-like leads using recursive filtration rules.  Input: Multi-record .sdf libraries from chemical databases.  Output: Standardized three-dimensional all-atom .pdbqt and .mol2 records.  Plaintext  Molecules (5,338) ➔ [PAINS Filter] ➔ 4,948 ➔ [Lipinski Ro5] ➔ 3,703 ➔ [Veber Constraints] ➔ 3,557 Passed Hits
The filtration sequence systematically applies PAINS substructure rules, Lipinski molecular weights ($MW < 500\text{ Da}$), $\log P < 5$, Veber rotational constraints ($\le 10$), polar surface properties ($TPSA \le 140\text{ \AA}^2$), and synthetic accessibility limits ($SA < 6.0$) to avoid hit inflation.  Module 2: Cross-Platform Consensus DockingLeverages four orthogonal grid-based search and scoring paradigms to eliminate software-specific scoring artifacts.  Input: 3,557 drug-like coordinates.  Output: sirtuin_consolidated_consensus_report_top500.csv.  Pythonimport pandas as pd
import numpy as np
from collections import defaultdict

FILE_LIST = [
    'b2_ad4_final_results_2.csv',
    'b2_dock6_final_results_2.csv',
    'b2_vina_final_results_2.csv',
    'b2_vinardo_final_results_2.csv',
]
REQUIRED_COLS = ['Full_Name', 'Binding_Energy', 'com_distance', 'rmsd']
TOP_N_CUTOFF = 500  
COM_DISTANCE_CUTOFF = 10.0  # Å
RMSD_CUTOFF = 3.0  # Å

def get_ranked_molecules_efficiently(file_path: str) -> pd.DataFrame:
    algorithm_name = file_path.split('_')[1]
    try:
        df = pd.read_csv(file_path, usecols=REQUIRED_COLS)
    except Exception as e:
        print(f"ERROR: Could not read {file_path}. Details: {e}")
        return pd.DataFrame()

    df['Molecule_ID'] = df['Full_Name'].str.split('_').str[0]
    df_filtered = df[(df['com_distance'] <= COM_DISTANCE_CUTOFF) & (df['rmsd'] <= RMSD_CUTOFF)].copy()
    
    if df_filtered.empty:
        return pd.DataFrame()

    df_unique_best_pose = (
        df_filtered.sort_values(by='Binding_Energy', ascending=True)
                   .drop_duplicates(subset=['Molecule_ID'], keep='first')
    )
    
    df_ranked = df_unique_best_pose.reset_index(drop=True).copy()
    df_full_report = df_ranked[['Molecule_ID', 'Full_Name', 'Binding_Energy', 'com_distance', 'rmsd']].copy()
    df_full_report.columns = ['Molecule_ID', 'Best_Pose_Full_Name', 'Best_Binding_Energy', 'Best_COM_Distance', 'Best_RMSD']
    df_full_report['Algorithm_Rank'] = df_full_report.index + 1
    df_full_report.to_csv(f'{algorithm_name}_full_ranked_report.csv', index=False)

    df_ranked[f'{algorithm_name}_Rank'] = df_ranked.index + 1
    df_ranked[f'{algorithm_name}_Energy'] = df_ranked['Binding_Energy']
    df_ranked[f'{algorithm_name}_COM_Dist'] = df_ranked['com_distance']
    df_ranked[f'{algorithm_name}_RMSD'] = df_ranked['rmsd']
    df_ranked[f'{algorithm_name}_Best_Pose'] = df_ranked['Full_Name']
    
    return df_ranked[['Molecule_ID', f'{algorithm_name}_Rank', f'{algorithm_name}_Energy',
                      f'{algorithm_name}_COM_Dist', f'{algorithm_name}_RMSD', f'{algorithm_name}_Best_Pose']].head(TOP_N_CUTOFF)

all_ranked_dfs = []
for file in FILE_LIST:
    df_top_n = get_ranked_molecules_efficiently(file)
    if not df_top_n.empty:
        all_ranked_dfs.append(df_top_n)

if len(all_ranked_dfs) > 0:
    df_consensus = all_ranked_dfs[0]
    for i in range(1, len(all_ranked_dfs)):
        df_consensus = pd.merge(df_consensus, all_ranked_dfs[i], on='Molecule_ID', how='outer')

    rank_cols = [col for col in df_consensus.columns if '_Rank' in col]
    df_consensus['Consensus_Count'] = df_consensus[rank_cols].notna().sum(axis=1)
    df_consensus_report = df_consensus[df_consensus['Consensus_Count'] >= 2].copy()
    df_consensus_report['Average_Rank'] = df_consensus_report[rank_cols].mean(axis=1)
    df_final_report = df_consensus_report.sort_values(by=['Consensus_Count', 'Average_Rank'], ascending=[False, True])
    df_final_report.to_csv(f'sirtuin_consolidated_consensus_report_top{TOP_N_CUTOFF}.csv', index=False)
    print("✅ Successfully generated combined consensus report matrix.")
Module 3: Binding Free Energy Thermodynamics ($MM/PBSA$)Calculates ensemble-averaged free energies ($\Delta G$) across short trajectories to isolate true thermodynamic affinity profiles[cite: 10].Input: GROMACS equilibration frame snapshots.  Output: Top_Individual_Poses_Ranked.csv evaluating relative cofactor activation potency.  $$\Delta\Delta G_{\text{NAD Activation}} = \Delta G_{\text{NAD (Holo)}} - \Delta G_{\text{NAD (Apo Baseline Ref)}}$$Pythonimport pandas as pd
import numpy as np

INPUT_FILE = "MMPBSA_Analysis_Ready.csv"
OUTPUT_FILE_BEST_POSE = "Top_Individual_Poses_Ranked.csv"

def find_best_individual_pose(input_file, output_file):
    try:
        df = pd.read_csv(input_file)
    except FileNotFoundError:
        print(f"Error: {input_file} not found.")
        return

    ref_nad_row = df[df['Molecule_Name'] == '4zzj_nad_p53']
    if ref_nad_row.empty or ref_nad_row['NAD_ΔTOTAL'].isnull().all():
        print("ERROR: Missing valid baseline reference (4zzj_nad_p53).")
        return
    
    delta_g_nad_apo = ref_nad_row['NAD_ΔTOTAL'].iloc[0]
    df['ΔΔG_NAD (Activation Efficacy)'] = df['NAD_ΔTOTAL'] - delta_g_nad_apo

    best_pose_df = df[df['Molecule_Name'] != '4zzj_nad_p53'].copy()
    best_pose_df.rename(columns={'ΔTOTAL': 'ΔG_Allost_Binding (kcal/mol)', 
                                 'NAD_ΔTOTAL': 'ΔG_NAD_Holo (kcal/mol)'}, inplace=True)
    
    best_pose_df.sort_values(by='ΔΔG_NAD (Activation Efficacy)', ascending=True, inplace=True)
    report_cols = ['Molecule_Name', 'Pose_Tag', 'ΔG_Allost_Binding (kcal/mol)', 
                   'ΔG_NAD_Holo (kcal/mol)', 'ΔΔG_NAD (Activation Efficacy)']
    
    final_report_df = best_pose_df[report_cols]
    final_report_df.to_csv(output_file, index=False, float_format='%.2f')
    print(f"🥇 Single Top Activator: {final_report_df.iloc[0]['Molecule_Name']} ({final_report_df.iloc[0]['Pose_Tag']})")

if __name__ == "__main__":
    find_best_individual_pose(INPUT_FILE, OUTPUT_FILE_BEST_POSE)
Module 4: 200 ns Explicit Solvent Molecular DynamicsSimulates equilibrium trajectories to characterize complex stability, global compactness ($R_g$), and correlated mechanical profiles[cite: 10].Configuration Specifications: $NSTEPS = 100,000,000$ (200 ns), $dt = 0.002\text{ ps}$, $T = 310\text{ K}$, $P = 1.0\text{ bar}$.  Forcefield Framework: AMBER ff19SB (protein) + CGenFF/C4 parameters (ligands) inside an explicit OPC water sphere.  Plaintext| System | Protein Ca-RMSD (nm) | Radius of Gyration Rg (nm) | Protein-NAD H-Bonds |
| :--- | :---: | :---: | :---: |
| 4ZZJ (Apo Baseline) | 0.7774 ± 0.1687 | 2.2817 | N/A |
| 4ZZJ + NAD+ Holo | 0.7545 ± 0.1716 | 2.2638 | 10.16 |
| Human 4TQ Reference | 0.5433 ± 0.1315 | 2.2902 | 11.26 |
| ZINC12660894 (Lead) | 0.4987 ± 0.1349 | 2.3085 | 9.58 |
| ZINC01078623 (Lead) | 0.5660 ± 0.1465 | 2.3468 | 10.83 |
Observations confirm that our leading discovery scaffolds restrict structural thermal fluctuations, matching or exceeding the stability profiles of co-crystallized industrial references.  Module 5: Potential of Mean Force Kinetics (Umbrella Sampling)Performs biased pulling simulations along a center-of-mass (COM) separation axis to determine the free energy dissociation barrier ($\Delta G^\ddagger$) of the cofactor.  Windows Alignment: ~38 distinct production windows spaced at $0.1\text{ nm}$ intervals, tracking displacements from $0.25\text{ nm}$ to $5.5\text{ nm}$.  Integration Method: Weighted Histogram Analysis Method (WHAM) parsing.  Plaintext  - Apo Complex NAD+ Escape Barrier (Baseline):  19.70 kJ/mol
  - 4TQ Reference Activator Barrier:             22.80 kJ/mol
  - ZINC01078623 (Pose 1 - DOCK6 Lead):         33.45 kJ/mol (+70% Potency Enhancement)
  - ZINC62004166 (Pose 3 - Vina Lead):          33.70 kJ/mol (+71% Potency Enhancement)
Module 6: Active Site Electronic Pre-Organization (DFT Analysis)Probes long-range quantum electronic shifts in charge distributions and orbital structures across local active site architectures[cite: 10].Theory Level Configuration: ! B3LYP def2-TZVP TightSCF SP RIJCOSX AutoAux D4 CPCM(water)[cite: 10, 11].Reactivity Descriptors Suite: Mayer Bond Orders, Intrabond Strength Index Weighted (IBSIW), Mulliken charges, and condensed nucleophilic Fukui ($f^-$) metrics parsed from Multiwfn wavefunctions[cite: 10, 11].Plaintext| Operational Complex System | Scissile C-N Mayer BO | Substrate Oxygen Charge (q) | C-O IBSIW Index |
| :--- | :---: | :---: | :---: |
| System A (Unmodulated Apo) | 0.8852 | -0.49647 | 66.38 |
| System B (4TQ Control Ref) | 0.8974 | -0.47769 | 0.41 |
| System E (ZINC12660894 Lead) | 0.8413 | -0.56672 | 55.89 |
💡 Core Electronic Discovery:
While industrial standards like 4TQ strengthen the scissile C-N bond, the prioritized ZINC12660894 platform drives a significant 5.0% reduction in glycosidic C-N bond order[cite: 10, 11]. This charge polarization destabilizes the ground state, lowering the intrinsic activation barrier for nicotinamide cleavage and accelerating catalytic turnover ($k_{\text{cat}}$)[cite: 10]. Simultaneously, it restores active site coordination, increasing the C-O IBSIW index to 55.89 to optimize the nucleophilic trajectory of the substrate oxygen[cite: 10, 11].🛠️ Instructions to Reproduce the Pipeline1. Environment Setup & DependenciesEnsure your cluster nodes or local execution workstations contain the following software configurations added to your system environment variables:Cheminformatics Core: Python 3.12+ with RDKit, Pandas, NumPy, and Matplotlib.  Parallel Virtual Docking: AutoDock-GPU (v1.6.7), DOCK6 (v6.13), and AutoDock Vina (v1.2.7).  Ensemble Trajectories: GROMACS 2025.1 compiled with CUDA acceleration.  Thermodynamics & Solvation: gmx_MMPBSA v1.6.4 and AmberTools 25.  Electronic Structure: ORCA 6.1.0 and Multiwfn v3.8.  2. Quick-Start Pipeline CommandsTo execute or reperform a specific component of the research pipeline, navigate to the respective directory path and run the target script modules:Bash# Step 1: Run virtual library drug-likeness filtration
cd 01_ligand_library_curation
python3 02_4_layer_filtering.py

# Step 2: Run multi-engine consensus docking analysis
cd ../02_consensus_docking/04_consensus_and_ranking
python3 analyze_csv_4_v2.py

# Step 3: Extract thermodynamic binding free energies
cd ../../03_mmpbsa_binding_energies
bash run_all_mmpbsa_LIG.sh
python3 3_Top_Individual_Poses_Ranked.py

# Step 4: Submit 200 ns GROMACS production runs
cd ../04_molecular_dynamics_200ns
gmx mdrun -v -deffnm step5_production -nb gpu

# Step 5: Run Umbrella Sampling PMF Extraction
cd ../05_umbrella_sampling_pmf
python3 umbrella_window_analyzer.py
python3 fel.py

# Step 6: Trigger ORCA single-point electronic evaluation calculations
cd ../06_quantum_mechanics_descriptors
orca active_site_cluster.inp > active_site_cluster.out