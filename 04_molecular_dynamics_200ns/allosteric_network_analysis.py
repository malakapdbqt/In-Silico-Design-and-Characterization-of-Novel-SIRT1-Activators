import os
import mdtraj as md
import numpy as np
import pandas as pd

BASE_DIR = "/home/kumara/malakafahs/gmx/b2_ligs"
def load_md(lig_dir):
    xtc = os.path.join(lig_dir, "gromacs", "traj_whole_centered.xtc")
    tpr = os.path.join(lig_dir, "gromacs", "step5_production.tpr")

    traj = md.load(xtc, top=tpr)
    return traj
def compute_dccm(traj):
    # Use only Cα atoms
    ca_indices = traj.top.select("name CA")
    ca_traj = traj.atom_slice(ca_indices)

    xyz = ca_traj.xyz - np.mean(ca_traj.xyz, axis=0)

    n_res = xyz.shape[1]
    dccm = np.zeros((n_res, n_res))

    for i in range(n_res):
        for j in range(n_res):
            num = np.sum(xyz[:,i,:] * xyz[:,j,:])
            den = np.sqrt(np.sum(xyz[:,i,:]**2) * np.sum(xyz[:,j,:]**2))
            dccm[i,j] = num / den

    return dccm
def extract_residue_correlations(dccm, traj):
    residues = [r for r in traj.top.residues]

    # get Cα mapping
    ca_res_index = []
    for r in residues:
        for atom in r.atoms:
            if atom.name == "CA":
                ca_res_index.append(r.index)

    ca_res_index = np.array(ca_res_index)

    # NAD Cα index
    nad_idx = np.where(ca_res_index == 15)[0][0]
    lig_idx = np.where(ca_res_index == 17)[0][0]

    df = pd.DataFrame({
        "Residue": ca_res_index,
        "Corr_with_NAD": dccm[:, nad_idx],
        "Corr_with_LIG": dccm[:, lig_idx]
    })

    return df
def run_all():
    results = []

    for folder in os.listdir(BASE_DIR):
        if folder.startswith("ZINC"):
            lig_path = os.path.join(BASE_DIR, folder)

            # find subdir containing gromacs
            for root, dirs, files in os.walk(lig_path):
                if "gromacs" in dirs:
                    gmx_dir = os.path.join(root, "gromacs")
                    break

            print("Processing:", folder)

            traj = load_md(root)
            dccm = compute_dccm(traj)
            df = extract_residue_correlations(dccm, traj)

            df.to_csv(os.path.join(lig_path, f"{folder}_DCCM.csv"), index=False)

