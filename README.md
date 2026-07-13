# Physics-Based Discovery and Mechanistic Characterization of Novel SIRT1 Allosteric Activators

An end-to-end reproducible, multi-scale computational screening and molecular modeling pipeline tracking the discovery of high-efficacy sirtuin 1 (SIRT1) allosteric regulators. This suite integrates high-throughput cheminformatics curation, multi-engine consensus docking, explicit solvent Molecular Dynamics (MD) trajectories, binding free energy thermodynamics ($G$), reaction coordinate kinetics via Umbrella Sampling (US), and quantum electronic pre-organization analysis.

## 🧬 Scientific Pipeline Overview

The active site architecture and structural configurations are derived from the crystallographic coordinate profiles of the human **SIRT1 enzyme complex (PDB ID: 4ZZJ)**. The pipeline models the complex in three operational frameworks:
* **System A (Apo-Ref Baseline):** Human SIRT1 modeled exclusively with native substrate fragments ($NAD^+$ and acetylated p53 peptide) without an allosteric modifier present.
* **System B (Positive Control Baseline):** Human SIRT1 complexed with the native reference activator (Indol derivatives) to establish comparative target milestones.
* **Systems C–F (Discovery Binders):** Lead allosteric candidate systems selected via multi-parameter sorting matrices for advanced validation testing.

---

## 📂 Repository Architecture & Directory Mapping

```text
.
├── 01_ligand_library_curation          # Stage 1: Initial virtual compound screening library prep
│   ├── rules                           # Pipeline operational screening rule logic
│   └── logs                            # High-throughput molecular filtration datalog matrices
├── 02_consensus_docking                # Stage 2: Cross-platform multi-engine structural docking execution
│   ├── 01_ligand_preparation           # Frontend preparation protocols and coordinate format conversions
│   ├── 02_docking_execution_data_arrange # Operational execution scripts and runtime configurations
│   │   ├── ad4                         # AutoDock4-GPU runtime engines and forcefield parameters
│   │   ├── dock6                       # UCSF DOCK6 structural gridbox, sphere generation, and scoring
│   │   ├── vina                        # AutoDock Vina multi-threaded screening paths
│   │   └── vinardo                     # Vinardo scoring matrix execution directories
│   ├── 03_Generate_Fingerprints        # Structural fragment fingerprints and raw CSV extraction
│   └── 04_consensus_and_ranking        # Intersecting candidate hits evaluation and consensus filters
├── 03_mmpbsa_binding_energies          # Stage 3: Thermodynamic binding free energy triage evaluations
├── 04_molecular_dynamics_200ns         # Stage 4: 200ns explicit solvent trajectory production setups
├── 05_umbrella_sampling_pmf            # Stage 5: Pulling paths and Potential of Mean Force (PMF) kinetics
├── 06_quantum_mechanics_descriptors    # Stage 6: active site electronic pre-organization (ORCA)
└── 07_admet_predictions                # Stage 7: In silico ADMET profiles and druggability matrices

![Python](https://img.shields.io/badge/Python-3.8%2B-blue?logo=python&logoColor=white)
![GROMACS](https://img.shields.io/badge/GROMACS-2022+-orange?logo=gromacs&logoColor=white)
![ORCA](https://img.shields.io/badge/ORCA-QM%20Theory-purple)
![License](https://img.shields.io/badge/License-MIT-green)
