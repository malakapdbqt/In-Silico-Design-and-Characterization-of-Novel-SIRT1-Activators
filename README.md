# 🧬 In Silico Design and Characterization of Novel SIRT1 Activators

## 📌 Overview
This repository contains a **complete, end-to-end computational drug discovery pipeline** developed for the identification and characterization of **novel SIRT1 allosteric activators**.

The workflow integrates:
- Ligand library curation
- Multi-engine consensus docking
- Binding free energy calculations (MM/PBSA)
- Molecular dynamics simulations (200 ns)
- Umbrella sampling (PMF)
- Quantum mechanical descriptor analysis
- ADMET profiling

This is not just a script collection — it is a **fully reproducible research pipeline**, designed for scalability (HPC-ready), traceability, and scientific rigor.

---

## 🎯 Research Objective
To identify **high-confidence SIRT1 activators** targeting the **NAD⁺ allosteric binding site**, using a **multi-layer computational screening and validation strategy**.

---

## ⚙️ Workflow Architecture


01_ligand_library_curation
02_consensus_docking
03_mmpbsa_binding_energies
04_molecular_dynamics_200ns
05_umbrella_sampling_pmf
06_quantum_mechanics_descriptors
07_admet_predictions


Each module represents a **logically independent but sequentially connected stage**.

---

# 🔬 1. Ligand Library Curation

### 📂 Path:
`01_ligand_library_curation/`

### 🧠 Purpose:
To construct a **high-quality, drug-like ligand dataset** from raw chemical libraries.

### ⚙️ Pipeline Steps:
1. Merge and deduplicate SDF files
2. Apply **multi-layer filtering**
3. Convert structures to docking-ready formats (PDBQT)
4. Prepare ligands with proper protonation and charges

### 📜 Filtering Criteria:
| Filter | Condition |
|------|--------|
| PAINS | Removal of pan-assay interference compounds |
| Lipinski | MW < 500, logP < 5, HBD < 5, HBA < 10 |
| Veber | Rotatable bonds < 10, TPSA < 140 |
| SA Score | < 6.0 |

### 📊 Results:
- **Total molecules processed:** 5338  
- **Final curated library:** 3557  
- **Success rate:** 66.64%  

### 🧾 Key Scripts:
- `01.merge_unique_sdf.py`
- `02.4_layer_filtering.py`
- `03.sdf_to_pdbqt_by_id.py`
- `04.prepare_all_ligands.py`

---

# ⚡ 2. Consensus Docking

### 📂 Path:
`02_consensus_docking/`

### 🧠 Purpose:
To improve prediction reliability using **multiple docking engines** and eliminate algorithm-specific bias.

### 🧪 Docking Engines Used:
- AutoDock4 (GPU)
- AutoDock Vina
- Vinardo
- DOCK6

### 🔁 Workflow:
1. Ligand preparation
2. Docking execution (engine-specific)
3. Pose extraction and normalization
4. Feature extraction (COM distance, RMSD)
5. Consensus ranking

---

## 📊 Filtering Strategy (Critical Step)

Only poses satisfying:
- **COM distance ≤ 10 Å**
- **RMSD ≤ 3 Å**

➡️ This ensures **binding site relevance + pose stability**

---

## 🧠 Consensus Algorithm

- Top **500 ligands per engine** selected
- Molecules merged across engines
- **Consensus Count** computed (2–4 engines agreement)
- Final ranking based on:
  - Highest consensus level
  - Lowest average rank

---

## 🏆 Key Outcome

- All **3557 ligands docked in 4 engines**
- High-confidence candidates selected via:
  - **4-way consensus**
  - **3-way consensus**
  - **2-way consensus**

Example (Top Candidate):

ZINC14632781 → 4-way consensus
Average Rank: 89


---

## 📜 Key Scripts:
- `analyze_csv_4_v2.py` ⭐ (Core consensus logic)
- `build_results_summary.py`
- `collect_consensus_poses.py`

---

# 💧 3. MM/PBSA Binding Energy Calculations

### 📂 Path:
`03_mmpbsa_binding_energies/`

### 🧠 Purpose:
To estimate **binding free energies (ΔG binding)** using post-docking trajectories.

### ⚙️ Method:
- AMBER-based MM/PBSA
- Frame sampling:
  - Start: 20
  - End: 101
  - Interval: 2

### 📈 Outputs:
- ΔG binding decomposition
- Ligand ranking refinement
- Prioritization maps

### 📜 Key Scripts:
- `1_MMPBSA_Combined_Results.py`
- `2_MMPBSA_Analysis_Ready.py`
- `prioratize2.py`

---

# 🌊 4. Molecular Dynamics Simulation (200 ns)

### 📂 Path:
`04_molecular_dynamics_200ns/`

### 🧠 Purpose:
To evaluate **dynamic stability and conformational behavior** of ligand–protein complexes.

### ⚙️ Engine:
- GROMACS

### 🔬 Analyses Performed:
- RMSD / RMSF
- DCCM (Dynamic Cross-Correlation)
- Allosteric network analysis

### 📜 Key Scripts:
- `run_all_sims_gmx.sh`
- `allosteric_network_analysis.py`
- `plot_dccm_residue.py`

---

# ⛰️ 5. Umbrella Sampling & PMF

### 📂 Path:
`05_umbrella_sampling_pmf/`

### 🧠 Purpose:
To calculate **Potential of Mean Force (PMF)** for ligand binding/unbinding.

### ⚙️ Method:
- Pulling simulations
- Window-based sampling
- WHAM analysis

### 📁 Outputs:
- PMF profiles
- Binding energy landscapes

---

# ⚛️ 6. Quantum Mechanical Descriptor Analysis

### 📂 Path:
`06_quantum_mechanics_descriptors/`

### 🧠 Purpose:
To understand **electronic-level interactions** influencing binding.

### ⚙️ Tools:
- Multiwfn

### 🔬 Descriptors:
- Fukui functions
- Mayer bond order
- Charge distribution
- Orbital interactions

---

# 🧪 7. ADMET Prediction

### 📂 Path:
`07_admet_predictions/`

### 🧠 Purpose:
To evaluate **drug-likeness and pharmacokinetic feasibility**

### 📊 Includes:
- Absorption
- Distribution
- Metabolism
- Excretion
- Toxicity

---

# 🔗 Full Workflow Summary


Library → Filtering → Docking (4 Engines)
→ Consensus Selection
→ MM/PBSA
→ Molecular Dynamics (200 ns)
→ Umbrella Sampling (PMF)
→ QM Descriptors
→ ADMET


---

# 🚀 How to Reuse This Workflow

## Step 1: Ligand Preparation
Run:
```bash
python 01.merge_unique_sdf.py
python 02.4_layer_filtering.py
python 04.prepare_all_ligands.py
Step 2: Docking

Run engine-specific scripts:

bash 01.run_autodock_gpu.sh
bash 01.dockvina.sh
python 01.run_rigid_dock.py
Step 3: Consensus Ranking
python analyze_csv_4_v2.py
Step 4: Post-Processing
MM/PBSA
MD simulations
PMF analysis
🧠 Key Scientific Strengths

✔ Multi-engine consensus docking (reduces bias)
✔ Strict structural filtering (COM + RMSD)
✔ Multi-scale validation (MM → MD → QM)
✔ End-to-end reproducibility
✔ HPC-compatible architecture

⚠️ Notes
Designed for Linux/HPC environments
Requires:
GROMACS
AutoDock tools
DOCK6
Python (RDKit, pandas, numpy)
Large datasets are not fully uploaded (storage limitations)
📜 License

See LICENSE file.

👤 Author

Malaka Sandaruwan
B.Pharm (Hons) – Computational Drug Discovery

⭐ Final Remark

This repository represents a complete computational drug discovery pipeline, moving from raw chemical space to mechanistic-level validation.

If you understand this workflow — you are not just running docking,
you are doing real drug discovery.