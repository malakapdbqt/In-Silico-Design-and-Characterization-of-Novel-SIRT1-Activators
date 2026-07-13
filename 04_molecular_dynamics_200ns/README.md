# SIRT1 Filtered Candidates - 200ns Long Trajectory Production Pipeline

This directory contains the production configurations, runtime setup scripts, and dynamic correlation tools used to extend target allosteric compound trajectories to 200ns.

## 🏃 Execution Workflow Sequence
1. Execute `./copy_gmx_files.sh` to generate standardized directory trees for your selected target compounds.
2. Run `./change_nstep.sh` to update coordinate interval targets to hit the 200ns parameter criteria inside your input parameters.
3. The `step5_production.mdp` file handles explicit system calculations (thermostats, steps, barostat).
4. Launch `./run_all_sims_gmx.sh` to submit jobs across cluster nodes.
5. Post-process trajectory dynamics by evaluating cross-correlation profiles with `./dccm_cross3.sh` and `allosteric_network_analysis.py`.
