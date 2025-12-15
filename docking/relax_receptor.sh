#!/bin/bash
#SBATCH --job-name=relax_receptor
#SBATCH --output=logs/relax_receptor_%j.out
#SBATCH --error=logs/relax_receptor_%j.err
#SBATCH --partition=short
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G

# Receptor Relaxation using Rosetta Container
# ============================================
# Relaxes the receptor structure before docking

set -e

module load apptainer/1.4.1-gcc-13.4.0

WORK_DIR="/home/kuhfeldr/hpc.cqls/docking"
RECEPTOR_DIR="${WORK_DIR}/receptor"
ROSETTA_SIF="/home/kuhfeldr/rosetta.sif"

cd "${WORK_DIR}"

# Find receptor - must be exactly one PDB file in receptor folder
RECEPTOR_FILES=($(ls ${RECEPTOR_DIR}/*.pdb 2>/dev/null))
if [ ${#RECEPTOR_FILES[@]} -eq 0 ]; then
    echo "ERROR: No receptor PDB file found in ${RECEPTOR_DIR}/"
    exit 1
elif [ ${#RECEPTOR_FILES[@]} -gt 1 ]; then
    echo "ERROR: Multiple receptor PDB files found in ${RECEPTOR_DIR}/"
    echo "  Found: ${RECEPTOR_FILES[*]}"
    echo "  Please keep only one receptor file in the receptor folder."
    exit 1
fi
RECEPTOR="${RECEPTOR_FILES[0]}"
RECEPTOR_NAME=$(basename "${RECEPTOR}" .pdb)

echo "=========================================="
echo "Receptor Relaxation"
echo "=========================================="
echo "Started: $(date)"
echo "Receptor: ${RECEPTOR}"
echo ""

apptainer exec --bind "${WORK_DIR}:${WORK_DIR}" "${ROSETTA_SIF}" \
    relax.default.linuxgccrelease \
    -s "${RECEPTOR}" \
    -relax:quick \
    -nstruct 1 \
    -out:path:pdb "${RECEPTOR_DIR}" \
    -out:prefix "relaxed_"

echo ""
echo "=========================================="
echo "Relaxation complete!"
echo "=========================================="
echo "Completed: $(date)"
ls -la "${RECEPTOR_DIR}/relaxed_"*.pdb 2>/dev/null || echo "No relaxed files found"
