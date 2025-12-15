#!/bin/bash
#SBATCH --job-name=analyze_md
#SBATCH --output=logs/analyze_%j.out
#SBATCH --error=logs/analyze_%j.err
#SBATCH --partition=short
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --mem=32G

# Analyze a single simulation directory
# Usage: sbatch submit_analyze_single.sh /path/to/simulation_dir

module load apptainer/1.4.1-gcc-13.4.0

SIM_DIR="$1"
if [ -z "$SIM_DIR" ]; then
    echo "ERROR: No simulation directory provided"
    exit 1
fi

cd /home/kuhfeldr/hpc.cqls/md
./analyze_md.sh "$SIM_DIR"
