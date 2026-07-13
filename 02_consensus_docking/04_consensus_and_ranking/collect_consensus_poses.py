import os
import shutil
import pandas as pd

# === Paths ===
base_dir = "/home/kumara/malakafahs/dock/allosteric"
output_dir = "/home/kumara/malakafahs/ML_models/NEW_model/b2_results/consensus_poses"

# Create output folder if not exists
os.makedirs(output_dir, exist_ok=True)

# Docking source directories
dock_paths = {
    "ad4": os.path.join(base_dir, "ad4_gpu_dock", "poses_mol2"),
    "dock6": os.path.join(base_dir, "dock6", "separated_mol2"),
    "vina": os.path.join(base_dir, "vina_gpu", "poses_mol2"),
    "vinardo": os.path.join(base_dir, "vinardo", "poses_mol2"),
}

# === Read CSV ===
df = pd.read_csv("sirtuin_consolidated_consensus_report_top500.csv")

# Filter only rows with Consensus_Count = 4
df_filtered = df[df["Consensus_Count"] == 4]

# === Loop through and copy files ===
for idx, row in df_filtered.iterrows():
    molecule = row["Molecule_ID"]

    # Each best pose filename (without extension)
    poses = {
        "ad4": row["ad4_Best_Pose"],
        "dock6": row["dock6_Best_Pose"],
        "vina": row["vina_Best_Pose"],
        "vinardo": row["vinardo_Best_Pose"],
    }

    for engine, pose_name in poses.items():
        src_file = os.path.join(dock_paths[engine], f"{pose_name}.mol2")
        dest_file = os.path.join(output_dir, f"{pose_name}.mol2")

        if os.path.exists(src_file):
            shutil.copy2(src_file, dest_file)
            print(f"✅ Copied {pose_name}.mol2 from {engine}")
        else:
            print(f"⚠️ Missing: {pose_name}.mol2 in {engine}")

print("\n✅ Copying complete! Files saved in:", output_dir)

