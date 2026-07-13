#!/bin/bash
# This script runs all four rerun scripts sequentially

echo "Starting final_b2_vinardo_rerun.py..."
python3 final_b2_vinardo_rerun.py
echo "Completed final_b2_vinardo_rerun.py"

echo "Starting final_b2_vina_gpu_rerun.py..."
python3 final_b2_vina_gpu_rerun.py
echo "Completed final_b2_vina_gpu_rerun.py"

echo "Starting final_b2_ad4_gpu_rerun.py..."
python3 final_b2_ad4_gpu_rerun.py
echo "Completed final_b2_ad4_gpu_rerun.py"

echo "Starting final_b2_dock6.py..."
python3 final_b2_dock6.py
echo "Completed final_b2_dock6.py"

echo "All rerun scripts completed ✅"

