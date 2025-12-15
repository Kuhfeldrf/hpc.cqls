#!/bin/bash
#SBATCH --job-name=analyze_md
#SBATCH --output=analyze_md_%j.out
#SBATCH --error=analyze_md_%j.err
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:a30:1
#SBATCH --mem=16G
#SBATCH --time=2:00:00
# Note: Container requires GPU node for initialization, but analysis runs on CPU

echo "=========================================="
echo "GROMACS MD Analysis Job"
echo "=========================================="
echo "Job started at: $(date)"
echo "Running on node: $(hostname)"
echo "SLURM Job ID: $SLURM_JOB_ID"
echo ""

# Load apptainer module
module load apptainer/1.4.1-gcc-13.4.0

# Change to md directory
cd /home/kuhfeldr/hpc.cqls/md

# Run analysis script
# Pass all command line arguments to the analysis script
./analyze_md.sh "$@"

echo ""
echo "=========================================="
echo "Analysis job completed at: $(date)"
echo "=========================================="

