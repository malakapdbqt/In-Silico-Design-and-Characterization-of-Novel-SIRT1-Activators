# SIRT1 Allosteric Hits - Binding Free Energy (MMPBSA) Pipeline

This directory handles the processing, aggregation, and ranking of binding free energies from simulation frames.

## 🏃 Execution Workflow Sequence
1. Ensure `mmpbsa.in` contains your correct structural parameters and grid choices.
2. Execute `./prep_trajectories.sh` to unwrap, fix PBC issues, and extract specific frames.
3. Launch `./run_all_mmpbsa_LIG.sh` to trigger the analytical engine calculation loops.
4. Run `python3 1_MMPBSA_Combined_Results.py` followed by `python3 2_MMPBSA_Analysis_Ready.py` to structure the energy matrices.
5. Apply `python3 3_Top_Individual_Poses_Ranked.py` and `prioratize2.py` to calculate thermodynamic contributions and isolate high-affinity leads.
