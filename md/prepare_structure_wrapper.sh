#!/bin/bash
# Wrapper script to submit prepare_structure.sh with outputs in structure directory

if [ -z "$1" ]; then
    echo "Error: Please provide a PDB filename as an argument"
    echo "Usage: sbatch prepare_structure_wrapper.sh <structure_name.pdb>"
    exit 1
fi

INPUT_NAME=$1
STRUCTURE_NAME="${INPUT_NAME%.*}"
BASE_DIR="/fs1/home/pi/kuhfeldr"
MD_DIR="${BASE_DIR}/md"
PREP_DIR="${MD_DIR}/${STRUCTURE_NAME}"

# Create preparation directory
mkdir -p "${PREP_DIR}"

# Copy the preparation script to the directory
cp "${MD_DIR}/prepare_structure.sh" "${PREP_DIR}/"

# Submit from within the structure directory so outputs go there
cd "${PREP_DIR}"
JOB_OUTPUT=$(sbatch --job-name="prep_${STRUCTURE_NAME}" \
       --output="prep_${STRUCTURE_NAME}_%j.out" \
       --error="prep_${STRUCTURE_NAME}_%j.err" \
       --partition=cqls_gpu-1080 \
       --account=cqls \
       --nodes=1 \
       --ntasks=4 \
       --cpus-per-task=4 \
       --gres=gpu:1 \
       --mem=32G \
       --time=12:00:00 \
       prepare_structure.sh "${MD_DIR}/${INPUT_NAME}" 2>&1)

cd "${MD_DIR}"

# Extract and display job ID
echo "$JOB_OUTPUT" | grep -E "[0-9]+" | head -1

