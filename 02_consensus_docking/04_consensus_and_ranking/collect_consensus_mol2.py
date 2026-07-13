#!/usr/bin/env python3
"""
collect_consensus_mol2.py
Copy .mol2 best-pose files for rows with Consensus_Count == 4 from engine pose folders
into a single destination directory. Produces a log of copied and missing files.

Usage:
  python3 collect_consensus_mol2.py       # real mode (will copy)
  python3 collect_consensus_mol2.py --dry # dry-run (no copying, just report)
"""

import csv
import os
import shutil
import argparse
from pathlib import Path

# ===== CONFIG - edit if your paths are different =====
csv_path = "sirtuin_consolidated_consensus_report_top500.csv"

# Source directories (as you described)
src_dirs = {
    "ad4": Path("/home/kumara/malakafahs/dock/allosteric/ad4_gpu_dock_2/poses_mol2"),
    "dock6": Path("/home/kumara/malakafahs/dock/allosteric/dock6/separated_mol2"),
    "vina": Path("/home/kumara/malakafahs/dock/allosteric/vina_gpu_2/poses_mol2"),
    "vinardo": Path("/home/kumara/malakafahs/dock/allosteric/vinardo_2/poses_mol2"),
}

# Destination directory (will create)
dest_base = Path("/home/kumara/malakafahs/ML_models/NEW_model/consensus_poses")

# Map CSV column -> engine key used above
col_to_engine = {
    "ad4_Best_Pose": "ad4",
    "dock6_Best_Pose": "dock6",
    "vina_Best_Pose": "vina",
    "vinardo_Best_Pose": "vinardo",
}

# ====================================================


def norm_pose_name(p):
    """Clean pose string from whitespace / surrounding quotes"""
    if p is None:
        return ""
    return str(p).strip()


def main(dry_run=False):
    if not Path(csv_path).exists():
        print(f"ERROR: CSV file not found: {csv_path}")
        return 1

    dest_base.mkdir(parents=True, exist_ok=True)

    # create per-engine subdirs to avoid name collisions (optional)
    per_engine_dirs = {eng: dest_base / eng for eng in src_dirs}
    for d in per_engine_dirs.values():
        d.mkdir(exist_ok=True)

    copied = []
    missing = []
    total_to_process = 0

    with open(csv_path, newline='', encoding='utf-8') as fh:
        reader = csv.DictReader(fh)
        # Ensure required columns exist
        required = ["Consensus_Count"] + list(col_to_engine.keys())
        for r in required:
            if r not in reader.fieldnames:
                print(f"ERROR: CSV missing required column: {r}")
                return 2

        for row in reader:
            try:
                cc = row.get("Consensus_Count", "").strip()
            except Exception:
                cc = row.get("Consensus_Count", "")
            # only process rows where Consensus_Count == 4 (string or int)
            if cc == "" or int(float(cc)) != 4:
                continue

            # for each engine column, attempt to copy
            for col, engine in col_to_engine.items():
                pose = norm_pose_name(row.get(col, ""))
                if not pose or pose.upper() == "N/A":
                    continue

                total_to_process += 1

                # pose in csv likely has no extension; append .mol2 if missing
                # but also allow if CSV already includes .mol2
                if pose.lower().endswith(".mol2"):
                    filename = pose
                else:
                    filename = f"{pose}.mol2"

                src_path = src_dirs[engine] / filename
                dest_path = per_engine_dirs[engine] / filename

                if src_path.exists():
                    if dry_run:
                        print(f"[DRY] Would copy: {src_path} -> {dest_path}")
                        copied.append((src_path, dest_path))
                    else:
                        try:
                            shutil.copy2(src_path, dest_path)
                            print(f"Copied: {src_path.name} -> {dest_path}")
                            copied.append((src_path, dest_path))
                        except Exception as e:
                            print(f"ERROR copying {src_path} -> {dest_path}: {e}")
                            missing.append((src_path, "copy-failed", str(e)))
                else:
                    print(f"Missing: {src_path}")
                    missing.append((src_path, "not-found"))

    # summary
    print("\n=== Summary ===")
    print(f"Total pose files targeted     : {total_to_process}")
    print(f"Files successfully (or dry-listed): {len(copied)}")
    print(f"Missing / failed               : {len(missing)}")
    if missing:
        miss_log = dest_base / "missing_consensus_poses.log"
        with open(miss_log, "w", encoding="utf-8") as mfh:
            for item in missing:
                if len(item) == 2:
                    mfh.write(f"{item[0]}    {item[1]}\n")
                else:
                    mfh.write(f"{item[0]}    {item[1]}    {item[2]}\n")
        print(f"Missing list written to: {miss_log}")

    print(f"Destination root: {dest_base}")
    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Collect consensus .mol2 poses (Consensus_Count==4)")
    parser.add_argument("--dry", action="store_true", help="Dry run (no files copied)")
    args = parser.parse_args()
    exit(main(dry_run=args.dry))

