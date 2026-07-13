# In Silico Design and Characterization of Novel SIRT1 Activators

A complete, reproducible computational pipeline for the discovery and characterization of novel small-molecule SIRT1 activators — spanning ligand library curation, consensus molecular docking, binding free-energy estimation, molecular dynamics, umbrella sampling, quantum mechanical descriptor generation, and ADMET prediction.

This repository documents every stage of the workflow used in the underlying thesis project, with scripts organized by pipeline step so the full process can be understood, audited, and reused.

---

## Table of Contents

- [Overview](#overview)
- [Pipeline Summary](#pipeline-summary)
- [Repository Structure](#repository-structure)
- [Requirements](#requirements)
- [Workflow Details](#workflow-details)
  - [1. Ligand Library Curation](#1-ligand-library-curation)
  - [2. Consensus Docking](#2-consensus-docking)
  - [3. MM-PBSA Binding Energies](#3-mm-pbsa-binding-energies)
  - [4. Molecular Dynamics (200 ns)](#4-molecular-dynamics-200-ns)
  - [5. Umbrella Sampling / PMF](#5-umbrella-sampling--pmf)
  - [6. Quantum Mechanics Descriptors](#6-quantum-mechanics-descriptors)
  - [7. ADMET Predictions](#7-admet-predictions)
- [Reproducing the Pipeline](#reproducing-the-pipeline)
- [Citation](#citation)
- [License](#license)

---

## Overview

SIRT1 (Sirtuin 1) is an NAD⁺-dependent deacetylase implicated in aging, metabolic regulation, and several disease pathways, making it a high-value target for small-molecule activator design. This project implements a **multi-stage, consensus-driven computational funnel** that progressively filters a large virtual ligand library down to a small set of high-confidence candidate activators, validated through increasingly rigorous and computationally expensive methods:
Library curation → Consensus docking (4 engines) → MM-PBSA rescoring
→ 200 ns MD → Umbrella sampling (PMF) → QM descriptors → ADMET filtering

Each stage acts as a filter, so only the most promising ligands are carried forward to the next, more expensive stage — a standard funnel strategy in structure-based virtual screening.

---

## Pipeline Summary

| Stage | Folder | Purpose | Key Tools |
|---|---|---|---|
| 1 | `01_ligand_library_curation` | Merge, filter, and prepare a clean ligand library | RDKit / OpenBabel, custom Python |
| 2 | `02_consensus_docking` | Dock ligands using 4 independent engines and rank by consensus | AutoDock-GPU, AutoDock Vina, Vinardo, DOCK6 |
| 3 | `03_mmpbsa_binding_energies` | Rescore top poses with MM-PBSA free energy calculations | AMBER / gmx_MMPBSA |
| 4 | `04_molecular_dynamics_200ns` | Validate binding stability over long-timescale MD | GROMACS |
| 5 | `05_umbrella_sampling_pmf` | Compute potential of mean force (PMF) for ligand (un)binding | GROMACS (pull code) + WHAM |
| 6 | `06_quantum_mechanics_descriptors` | Generate electronic/QM descriptors for top hits | Multiwfn, Gaussian/ORCA outputs |
| 7 | `07_admet_predictions` | Predict pharmacokinetic and toxicity profiles | ADMET prediction tools |

---

## Repository Structure
.
├── 01_ligand_library_curation/
│   ├── logs/
│   └── rules/
│
├── 02_consensus_docking/
│   ├── 01_ligand_preparation/
│   ├── 02_docking_execution_data_arrange/
│   │   ├── ad4/
│   │   │   ├── dpf_gpf_files/
│   │   │   └── ligands_prepared_hs/
│   │   ├── dock6/
│   │   │   ├── doc6_2/
│   │   │   │   ├── 001_structure/
│   │   │   │   ├── 002_spheres/
│   │   │   │   ├── 003_gridbox/
│   │   │   │   ├── 004_energy_min/
│   │   │   │   ├── 005_rigid_dock/
│   │   │   │   ├── 006_flex_dock/
│   │   │   │   ├── 007_footprint/
│   │   │   │   ├── 008_virtual_screen/
│   │   │   │   ├── 009_gen_alg/
│   │   │   │   └── 010_de_novo/
│   │   │   └── ligands_mol2/
│   │   ├── vina/
│   │   │   └── all_prepared_ligand/
│   │   └── vinardo/
│   │       └── all_prepared_ligand/
│   ├── 03_Generate_Fingerprints/
│   ├── 03_raw_data_collection/
│   └── 04_consensus_and_ranking/
│
├── 03_mmpbsa_binding_energies/
├── 04_molecular_dynamics_200ns/
├── 05_umbrella_sampling_pmf/
│   ├── pmf_profiles/
│   └── pulling_mdp/
├── 06_quantum_mechanics_descriptors/
│   ├── Multiwfn_parses/
│   └── scripts/
├── 07_admet_predictions/
│   └── profiles/
│
├── LICENSE
└── README.md

**Folder-by-folder purpose:**

| Path | Contents |
|---|---|
| `01_ligand_library_curation/` | Library merge/filter/prep scripts · `logs/` filtering reports · `rules/` filtering flowchart & rule docs |
| `02_consensus_docking/01_ligand_preparation/` | Batch ligand prep shared across all 4 docking engines |
| `02_consensus_docking/02_docking_execution_data_arrange/ad4/` | AutoDock4-GPU: grid files, prepared ligands, run/rename scripts |
| `02_consensus_docking/02_docking_execution_data_arrange/dock6/` | DOCK6: staged pipeline (structure → spheres → gridbox → energy min → rigid/flex dock → footprint/screen/de novo) |
| `02_consensus_docking/02_docking_execution_data_arrange/vina/` | AutoDock Vina docking + prepared ligands |
| `02_consensus_docking/02_docking_execution_data_arrange/vinardo/` | Vinardo scoring function docking + prepared ligands |
| `02_consensus_docking/03_Generate_Fingerprints/` | Per-engine rerun + interaction fingerprint scripts |
| `02_consensus_docking/03_raw_data_collection/` | Aggregated raw docking outputs |
| `02_consensus_docking/04_consensus_and_ranking/` | Cross-engine consensus scoring, ranking, top-pose collection |
| `03_mmpbsa_binding_energies/` | MM-PBSA rescoring, prioritization mapping |
| `04_molecular_dynamics_200ns/` | 200 ns production MD, DCCM/allosteric network analysis |
| `05_umbrella_sampling_pmf/` | Umbrella sampling windows (`pulling_mdp/`) and PMF profiles (`pmf_profiles/`) |
| `06_quantum_mechanics_descriptors/` | QM descriptor scripts + Multiwfn parses |
| `07_admet_predictions/` | Per-compound ADMET prediction profiles |

> Large intermediate outputs (trajectories, docking maps, raw pose files) are excluded from version control — only scripts, configs, and summary data are tracked.
> **Note:** Large intermediate outputs (trajectories, docking maps, raw pose files) are excluded from version control via `.gitattributes`/`.gitignore` conventions — only scripts, configs, and summary data are tracked. Regenerate bulk outputs by rerunning the corresponding stage script(s).

---

## Requirements

| Tool | Used In | Notes |
|---|---|---|
| Python 3.9+ | All stages | RDKit, pandas, numpy, matplotlib |
| OpenBabel | Ligand prep | Format conversion (SDF ↔ PDBQT ↔ MOL2) |
| AutoDock-GPU | Stage 2 | GPU-accelerated AutoDock4 scoring |
| AutoDock Vina | Stage 2 | Vina & Vinardo scoring functions |
| DOCK6 (UCSF) | Stage 2 | Sphere generation, grid, rigid/flexible docking |
| AMBER / gmx_MMPBSA | Stage 3 | MM-PBSA binding free energy |
| GROMACS | Stages 4–5 | MD production runs, umbrella sampling pull code |
| WHAM | Stage 5 | PMF profile reconstruction |
| Multiwfn | Stage 6 | Wavefunction analysis, QM descriptor extraction |
| ADMET prediction tool(s) | Stage 7 | e.g., SwissADME / pkCSM / ADMETlab (specify which you used) |

Install Python dependencies:

```bash
pip install -r requirements.txt
```

*(Add a `requirements.txt` listing exact package versions used, for full reproducibility.)*

---

## Workflow Details

### 1. Ligand Library Curation
`01_ligand_library_curation/`

Builds a clean, non-redundant ligand library from raw SDF sources.

1. **`01.merge_unique_sdf.py`** — Merges multiple source SDF files and removes duplicate structures.
2. **`02.4_layer_filtering.py`** — Applies a four-layer filtering scheme (see `rules/thesis_filtering_pipeline_v2.pdf` for the full flowchart) to remove non-drug-like or unsuitable molecules.
3. **`03.sdf_to_pdbqt_by_id.py`** — Converts filtered ligands from SDF to PDBQT format, indexed by compound ID.
4. **`04.prepare_all_ligands.py`** / **`prepare_ligand4.py`** — Final ligand preparation (protonation, torsion tree setup) ahead of docking.

Logs and QC reports from the filtering run are saved to `logs/filtering_data_report.dat`.

**Run:**
```bash
cd 01_ligand_library_curation
python 01.merge_unique_sdf.py
python 02.4_layer_filtering.py
python 03.sdf_to_pdbqt_by_id.py
python 04.prepare_all_ligands.py
```

---

### 2. Consensus Docking
`02_consensus_docking/`

The curated ligand library is docked against the SIRT1 receptor structure using **four independent docking engines**, and results are combined into a consensus ranking to reduce single-algorithm bias.

#### 2.1 Ligand Preparation
`01_ligand_preparation/batch_prepare_ligands2.py` — batch-prepares ligands for all four docking engines simultaneously.

#### 2.2 Docking Execution
`02_docking_execution_data_arrange/`

- **`ad4/`** — AutoDock4-GPU pipeline:
  1. `01.run_autodock_gpu.sh` — runs GPU-accelerated docking using precomputed grid maps (`dpf_gpf_files/`)
  2. `02.extract_poses_from_dlg.sh` — extracts docked poses from `.dlg` output
  3. `03.pdbqt_to_mol2.sh` — converts poses to MOL2
  4. `04.rename_logs.sh` / `05.rename_mol2.sh` — standardizes output filenames for downstream parsing

- **`dock6/`** — DOCK6 pipeline, run in structured stages under `doc6_2/`:
  1. `001_structure` — receptor/ligand structure prep
  2. `002_spheres` — DOCK sphere generation (`sphgen`)
  3. `003_gridbox` — energy grid box generation
  4. `004_energy_min` — pre-dock energy minimization
  5. `005_rigid_dock` — rigid-body docking
  6. `006_flex_dock` — flexible ligand docking
  7. `007_footprint` – `010_de_novo` — reserved for footprint similarity, virtual screening, genetic algorithm, and de novo design extensions
  
  Orchestrated by `01.run_rigid_dock.py` → `02.collect_dock6_outputs.sh` → `03.seperate_mol2_2.sh` → `04.rename_logs.sh` / `05.rename_mol2.sh`.

- **`vina/`** — AutoDock Vina docking: `01.dockvina.sh` → `02.split_vina_outputs.sh` → `03.pdbqt_to_mol2.sh` → rename scripts. `progress_estimator.sh` tracks batch docking progress.

- **`vinardo/`** — Identical structure to `vina/`, using the Vinardo scoring function (`01.dock_all.sh` → `02.split_vinardo_outputs.sh` → …).

Each engine folder is self-contained: receptor file (`REC.pdbqt`), prepared ligands, run scripts, and output-renaming utilities, so any single engine can be rerun independently.

#### 2.3 Fingerprint Generation
`03_Generate_Fingerprints/` — per-engine rerun + interaction fingerprint scripts (`final_b2_*.py`) that recompute/verify scores and generate protein–ligand interaction fingerprints for each docking method, orchestrated by `run_all_reruns_b2.sh`.

#### 2.4 Consensus & Ranking
`04_consensus_and_ranking/`
1. `analyze_csv_4_v2.py` — parses and merges per-engine result CSVs
2. `build_results_summary.py` — builds a unified summary table across all four engines
3. `collect_consensus_mol2.py` / `collect_consensus_poses.py` — gathers top consensus poses for downstream MM-PBSA/MD stages

**Run (per engine, example — Vina):**
```bash
cd 02_consensus_docking/02_docking_execution_data_arrange/vina
./01.dockvina.sh
./02.split_vina_outputs.sh
./03.pdbqt_to_mol2.sh
```

---

### 3. MM-PBSA Binding Energies
`03_mmpbsa_binding_energies/`

Top consensus hits are rescored using MM-PBSA to obtain more physically rigorous binding free energy estimates.

1. `prep_trajectories.sh` — prepares MD trajectories/topologies for MM-PBSA input
2. `run_all_mmpbsa_LIG.sh` — batch-runs MM-PBSA (`mmpbsa.in` config) across all ligand complexes
3. `1_MMPBSA_Combined_Results.py` — aggregates raw MM-PBSA output
4. `2_MMPBSA_Analysis_Ready.py` — cleans/formats results for analysis
5. `3_Top_Individual_Poses_Ranked.py` — ranks poses by binding free energy
6. `consolidate_mmpbsa_step.py` — merges intermediate results
7. `prioratize2.py` / `plot_prioritization_map.py` — generates a prioritization map to select final candidates for MD

See stage-specific `README.md` inside this folder for parameter details.

---

### 4. Molecular Dynamics (200 ns)
`04_molecular_dynamics_200ns/`

Selected candidates undergo 200 ns all-atom MD simulations in GROMACS to assess binding stability and dynamic behavior.

1. `copy_gmx_files.sh` — stages simulation input files per ligand
2. `change_nstep.sh` — configures simulation length (`step5_production.mdp`)
3. `run_all_sims_gmx.sh` — batch-launches production MD runs
4. `dccm_cross3.sh` / `plot_dccm_residue.py` — dynamic cross-correlation matrix (DCCM) analysis
5. `allosteric_network_analysis.py` — identifies allosteric communication networks from MD trajectories

See stage-specific `README.md` inside this folder for GROMACS version and MDP parameter details.

---

### 5. Umbrella Sampling / PMF
`05_umbrella_sampling_pmf/`

Computes the potential of mean force (PMF) for ligand binding/unbinding along a reaction coordinate, using GROMACS pull code and WHAM reconstruction.

- `pulling_mdp/` — MDP configuration files for steered MD and umbrella sampling windows
- `pmf_profiles/` — WHAM-reconstructed PMF profiles and free-energy barrier plots

---

### 6. Quantum Mechanics Descriptors
`06_quantum_mechanics_descriptors/`

Generates electronic-structure descriptors (e.g., HOMO/LUMO gap, electrostatic potential, Fukui indices) for top candidates to support structure-activity rationalization.

- `scripts/` — QM job setup and descriptor extraction scripts
- `Multiwfn_parses/` — Multiwfn wavefunction analysis outputs

---

### 7. ADMET Predictions
`07_admet_predictions/`

Final pharmacokinetic/toxicity screening of the top candidate list.

- `profiles/` — per-compound ADMET prediction outputs (absorption, distribution, metabolism, excretion, toxicity)

---

## Reproducing the Pipeline

To rerun the full workflow end-to-end on a new receptor/ligand set:

1. Place raw ligand SDF files in the expected input location for `01_ligand_library_curation`.
2. Run stage 1 scripts sequentially (numbered order).
3. Prepare the receptor and run all four docking engines under `02_consensus_docking` (each is independent — can be parallelized).
4. Run fingerprinting and consensus ranking scripts to shortlist top hits.
5. Run MM-PBSA rescoring (`03_mmpbsa_binding_energies`) on shortlisted complexes.
6. Launch 200 ns MD (`04_molecular_dynamics_200ns`) for MM-PBSA-validated hits.
7. Run umbrella sampling (`05_umbrella_sampling_pmf`) on final candidates for PMF/free-energy profiles.
8. Generate QM descriptors (`06_quantum_mechanics_descriptors`) and ADMET profiles (`07_admet_predictions`) for the final candidate shortlist.

> Each numbered folder is designed to be runnable independently given its expected inputs — see per-stage `README.md` files (where present) for exact input/output specifications.

---

## Citation

If you use this workflow or any part of it in your own research, please cite:
[Author Name]. "In Silico Design and Characterization of Novel SIRT1 Activators."
[Thesis/Institution, Year]. GitHub repository: [repo URL]

*(Replace with your actual citation details.)*

---

## License

This project is licensed under the terms described in [`LICENSE`](./LICENSE).