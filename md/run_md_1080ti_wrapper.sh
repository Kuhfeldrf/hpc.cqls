#!/bin/bash
# Wrapper script to submit run_md_1080ti.sh with outputs in structure directory

if [ -z "$1" ]; then
    echo "Error: Please provide a PDB filename as an argument"
    echo "Usage: sbatch run_md_1080ti_wrapper.sh <structure_name.pdb>"
    exit 1
fi

INPUT_NAME=$1
STRUCTURE_NAME="${INPUT_NAME%.*}"
BASE_DIR="/fs1/home/pi/kuhfeldr"
MD_DIR="${BASE_DIR}/md"
SIM_DIR="${MD_DIR}/${STRUCTURE_NAME}"

# Create simulation directory
mkdir -p "${SIM_DIR}"

# Copy the MD script to the directory
cp "${MD_DIR}/run_md_1080ti.sh" "${SIM_DIR}/"

# Submit from within the structure directory so outputs go there
cd "${SIM_DIR}"
sbatch --job-name="md_${STRUCTURE_NAME}" \
       --output="md_${STRUCTURE_NAME}_%j.out" \
       --error="md_${STRUCTURE_NAME}_%j.err" \
       --partition=cqls_gpu-1080 \
       --account=cqls \
       --nodes=1 \
       --ntasks=1 \
       --cpus-per-task=4 \
       --gres=gpu:1 \
       --mem=32G \
       --time=288:00:00 \
       run_md_1080ti.sh "${MD_DIR}/${INPUT_NAME}"

cd "${MD_DIR}"

