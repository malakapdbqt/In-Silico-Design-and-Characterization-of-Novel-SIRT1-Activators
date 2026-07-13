#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.coordinates import TPR

# -----------------------------
# User parameters
# -----------------------------
cov_file = "covariance_matrix_ascii.dat"  # your 3N x 3N GROMACS covar output
tpr_file = "../ZINC01078623_pose1_dock6_b2/gromacs/step5_production.tpr"        # GROMACS TPR file for residue mapping
output_prefix = "DCCM_residue"            # prefix for output files

# -----------------------------
# Load covariance matrix
# -----------------------------
print("Loading covariance matrix...")
cov_3N = np.loadtxt(cov_file)
print(f"Original matrix shape: {cov_3N.shape}")

# -----------------------------
# Load TPR and map residues
# -----------------------------
print("Loading TPR to extract residue mapping...")
u = mda.Universe(tpr_file)
# select C-alpha atoms (or main atoms per residue)
ca_atoms = u.select_atoms("name CA")
N_residues = len(ca_atoms)
print(f"Number of residues (C-alpha): {N_residues}")

# -----------------------------
# Convert to residue-level covariance
# -----------------------------
cov_res = np.zeros((N_residues, N_residues))

for i, res_i in enumerate(ca_atoms.residues):
    for j, res_j in enumerate(ca_atoms.residues):
        # take the 3x3 block corresponding to x,y,z of both residues
        block_i = slice(i*3, i*3+3)
        block_j = slice(j*3, j*3+3)
        block = cov_3N[block_i, block_j]
        cov_res[i, j] = np.mean(block)

print(f"Residue-level covariance matrix shape: {cov_res.shape}")

# -----------------------------
# Convert to cross-correlation (DCCM)
# -----------------------------
# Cij = cov_ij / sqrt(cov_ii * cov_jj)
dccm = np.zeros_like(cov_res)
for i in range(N_residues):
    for j in range(N_residues):
        dccm[i, j] = cov_res[i, j] / np.sqrt(cov_res[i, i] * cov_res[j, j])

# -----------------------------
# Plot DCCM
# -----------------------------
plt.figure(figsize=(10, 8))
im = plt.imshow(dccm, cmap='bwr', vmin=-1, vmax=1)
plt.colorbar(im, label='Correlation coefficient')
plt.xlabel("Residue index")
plt.ylabel("Residue index")
plt.title("Residue-level Dynamical Cross-Correlation Matrix (DCCM)")
plt.tight_layout()
plt.savefig(f"{output_prefix}.png", dpi=300)
plt.savefig(f"{output_prefix}.eps")  # for vector graphics
plt.close()
print(f"Saved DCCM plot to {output_prefix}.png / {output_prefix}.eps")

# -----------------------------
# Save matrix as npy and txt
# -----------------------------
np.savetxt(f"{output_prefix}.txt", dccm, fmt="%.5f")
np.save(f"{output_prefix}.npy", dccm)
print(f"Saved DCCM matrix to {output_prefix}.txt and {output_prefix}.npy")

# -----------------------------
# Eigenvalue analysis (optional)
# -----------------------------
eigenvalues = np.linalg.eigvalsh(cov_res)
plt.figure(figsize=(8, 5))
plt.plot(np.arange(1, len(eigenvalues)+1), eigenvalues, 'o-', markersize=4)
plt.xlabel("Mode index")
plt.ylabel("Eigenvalue (nm^2)")
plt.title("Residue-level covariance eigenvalues")
plt.tight_layout()
plt.savefig(f"{output_prefix}_eigenvalues.png", dpi=300)
plt.savefig(f"{output_prefix}_eigenvalues.eps")
plt.close()
print(f"Saved eigenvalue plot to {output_prefix}_eigenvalues.png / .eps")

