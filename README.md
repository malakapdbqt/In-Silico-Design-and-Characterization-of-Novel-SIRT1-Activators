<h1 align="center">🧬 SIRT1 Activator Discovery Pipeline</h1>

<p align="center">
<b>End-to-End Computational Drug Discovery Workflow</b><br>
Docking → Free Energy → Molecular Dynamics → QM → ADMET
</p>

<p align="center">
<img src="https://img.shields.io/badge/Python-3.10-blue">
<img src="https://img.shields.io/badge/GROMACS-2023-orange">
<img src="https://img.shields.io/badge/AutoDock-4%20Engines-green">
<img src="https://img.shields.io/badge/Stage-Complete-success">
<img src="https://img.shields.io/badge/Type-In%20Silico-blueviolet">
</p>

---

## 📌 Overview

This project presents a **fully integrated computational drug discovery pipeline** for identifying **novel SIRT1 activators** targeting the **NAD⁺ allosteric binding site**.

Unlike isolated docking studies, this workflow applies **multi-scale validation**, combining:

- Multi-engine consensus docking  
- Thermodynamic validation (MM/PBSA)  
- Long-timescale molecular dynamics (200 ns)  
- Free energy profiling (PMF)  
- Quantum mechanical descriptor analysis  
- ADMET prediction  

> ⚡ Designed for **scientific rigor, reproducibility, and HPC scalability**

---

## 🧠 Scientific Strategy

- **Target:** SIRT1 NAD⁺ allosteric site  
- **Challenge:** Docking bias & false positives  
- **Solution:** Multi-layer validation pipeline  

### 🔬 Validation Stack

| Layer | Purpose |
|------|--------|
| Docking | Pose generation & ranking |
| Consensus | Reduce algorithm bias |
| MM/PBSA | Binding free energy estimation |
| MD Simulation | Stability & conformational dynamics |
| PMF | Binding/unbinding energetics |
| QM Analysis | Electronic interaction insights |
| ADMET | Drug-likeness evaluation |

---

## 🔁 Workflow Architecture

```mermaid
graph LR

A[📦 Ligand Library] --> B[🧹 Filtering]
B --> C[⚡ Docking - 4 Engines]
C --> D[🧠 Consensus Selection]
D --> E[💧 MM/PBSA]
E --> F[🌊 MD Simulation - 200 ns]
F --> G[⛰️ Umbrella Sampling]
G --> H[⚛️ QM Analysis]
H --> I[🧪 ADMET]
📊 Key Metrics
🧪 5338 compounds screened
🎯 3557 retained after filtering
⚡ 4 docking engines used
🏆 Top hit: ZINC14632781 (4-way consensus)
⚙️ Project Structure
01_ligand_library_curation/
02_consensus_docking/
03_mmpbsa_binding_energies/
04_molecular_dynamics_200ns/
05_umbrella_sampling_pmf/
06_quantum_mechanics_descriptors/
07_admet_predictions/
🔬 1. Ligand Library Curation
<details> <summary>View Details</summary>
Purpose

To generate a high-quality, drug-like ligand dataset

Filtering Criteria
Filter	Condition
PAINS	Removed
Lipinski	MW < 500, logP < 5
Veber	TPSA < 140
SA Score	< 6
Output
Initial: 5338
Final: 3557
Key Scripts
01.merge_unique_sdf.py
02.4_layer_filtering.py
04.prepare_all_ligands.py
</details>
⚡ 2. Consensus Docking
<details> <summary>View Details</summary>
Engines Used
AutoDock4 (GPU)
AutoDock Vina
Vinardo
DOCK6
Filtering Criteria
COM distance ≤ 10 Å
RMSD ≤ 3 Å
Strategy
Top 500 ligands per engine
Merge results
Compute consensus score (2–4 engines)
Outcome

High-confidence ligand selection based on agreement across engines.

</details>
💧 3. MM/PBSA Binding Energies
<details> <summary>View Details</summary>
AMBER-based calculations
Frame sampling: 20–101 (interval 2)
Output: ΔG binding & decomposition
</details>
🌊 4. Molecular Dynamics (200 ns)
<details> <summary>View Details</summary>
Engine
GROMACS
Analyses
RMSD / RMSF
DCCM
Allosteric network analysis
</details>
⛰️ 5. Umbrella Sampling & PMF
<details> <summary>View Details</summary>
Pulling simulations
Window sampling
WHAM analysis
Output: Free energy landscape
</details>
⚛️ 6. Quantum Mechanical Analysis
<details> <summary>View Details</summary>
Tool
Multiwfn
Descriptors
Fukui functions
Mayer bond order
Charge distribution
</details>
🧪 7. ADMET Prediction
<details> <summary>View Details</summary>

Evaluates:

Absorption
Distribution
Metabolism
Toxicity
</details>
🚀 Usage Guide
1. Ligand Preparation
python 01.merge_unique_sdf.py
python 02.4_layer_filtering.py
python 04.prepare_all_ligands.py
2. Docking
bash 01.run_autodock_gpu.sh
bash 01.dockvina.sh
python 01.run_rigid_dock.py
3. Consensus Ranking
python analyze_csv_4_v2.py
4. Downstream Analysis
MM/PBSA
MD simulations
PMF
QM
🧠 Key Strengths

✔ Multi-engine docking (bias reduction)
✔ Strict pose validation (COM + RMSD)
✔ Multi-scale validation (MM → MD → QM)
✔ HPC-ready architecture
✔ Reproducible pipeline

⚠️ Requirements
Linux / HPC environment
Python (RDKit, numpy, pandas)
GROMACS
AutoDock tools
DOCK6
👤 Author

Malaka Sandaruwan
B.Pharm (Hons) – Computational Drug Discovery

⭐ Final Note

This project represents a complete computational drug discovery pipeline, progressing from chemical space exploration to mechanistic validation.

If you understand this workflow, you are not just running docking —
you are performing real drug discovery.