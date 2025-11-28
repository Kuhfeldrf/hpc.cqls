#!/bin/bash
#SBATCH --job-name=md_1080ti_100ns
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --partition=cqls_gpu-1080
#SBATCH --account=cqls
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:1
#SBATCH --mem=32G
#SBATCH --time=288:00:00

# Check if filename argument is provided
if [ -z "$1" ]; then
    echo "Error: Please provide a filename/structure name as an argument"
    echo "Usage: sbatch run_md_1080ti.sh <structure_name>"
    exit 1
fi

# Strip file extension (e.g., .pdb, .gro) from the input filename
# Handle both relative and absolute paths
if [[ "$1" == /* ]]; then
    FULL_PATH="$1"
    INPUT_NAME=$(basename "$1")
else
    FULL_PATH="${MD_DIR}/$1"
    INPUT_NAME="$1"
fi
STRUCTURE_NAME="${INPUT_NAME%.*}"
BASE_DIR="/fs1/home/pi/kuhfeldr"
MD_DIR="${BASE_DIR}/md"
SIM_DIR="${MD_DIR}/${STRUCTURE_NAME}"

# Create simulation directory first
mkdir -p "${SIM_DIR}"

# Change to simulation directory so SLURM outputs go here
cd "${SIM_DIR}"

# Copy the script itself to the simulation directory
SCRIPT_NAME=$(basename "$0")
SCRIPT_PATH="${MD_DIR}/${SCRIPT_NAME}"
if [ -f "${SCRIPT_PATH}" ]; then
    cp "${SCRIPT_PATH}" "${SIM_DIR}/"
fi

# OpenCL environment (we KNOW this works on GTX 1080 Ti)
export GPU_DEVICE_ORDINAL=0
export OCL_ICD_VENDORS=/etc/OpenCL/vendors
export OPENCL_VENDOR_PATH=/etc/OpenCL/vendors

echo "=========================================="
echo "GTX 1080 Ti Production MD - 100 ns"
echo "Structure: ${STRUCTURE_NAME}"
echo "Job ID: ${SLURM_JOB_ID}"
echo "=========================================="
echo "Job started: $(date)"
echo "Node: $(hostname)"
echo ""

nvidia-smi -L
echo ""

# Copy all input files (.mdp settings files)
echo "Copying input files to ${SIM_DIR}..."
cp "${MD_DIR}"/*.mdp "${SIM_DIR}/" 2>/dev/null || true
echo "Input files copied."

# Copy topology and structure files if they exist in md directory
if [ -f "${MD_DIR}/topol.top" ]; then
    cp "${MD_DIR}/topol.top" "${SIM_DIR}/"
fi
if [ -f "${MD_DIR}/topol_Protein_chain_"*.itp ]; then
    cp "${MD_DIR}/topol_Protein_chain_"*.itp "${SIM_DIR}/" 2>/dev/null || true
fi
if [ -f "${MD_DIR}/posre_Protein_chain_"*.itp ]; then
    cp "${MD_DIR}/posre_Protein_chain_"*.itp "${SIM_DIR}/" 2>/dev/null || true
fi

# Copy the input .tpr file (check structure directory first, then md directory)
if [ -f "${SIM_DIR}/md.tpr" ]; then
    echo "Found md.tpr in structure directory"
elif [ -f "${MD_DIR}/${STRUCTURE_NAME}/md.tpr" ]; then
    cp "${MD_DIR}/${STRUCTURE_NAME}/md.tpr" "${SIM_DIR}/"
    echo "Copied md.tpr from ${STRUCTURE_NAME}/ directory"
elif [ -f "${MD_DIR}/md.tpr" ]; then
    cp "${MD_DIR}/md.tpr" "${SIM_DIR}/"
    echo "Copied md.tpr from md directory"
elif [ -f "${MD_DIR}/${STRUCTURE_NAME}.tpr" ]; then
    cp "${MD_DIR}/${STRUCTURE_NAME}.tpr" "${SIM_DIR}/md.tpr"
    echo "Copied ${STRUCTURE_NAME}.tpr as md.tpr"
else
    echo "ERROR: No .tpr file found!"
    echo "Please prepare the structure first with:"
    echo "  sbatch prepare_structure.sh ${INPUT_NAME}"
    echo ""
    echo "Or ensure md.tpr exists in one of these locations:"
    echo "  - ${SIM_DIR}/md.tpr"
    echo "  - ${MD_DIR}/${STRUCTURE_NAME}/md.tpr"
    echo "  - ${MD_DIR}/md.tpr"
    exit 1
fi

# Copy the checkpoint file if it exists (for continuation)
if [ -f "${MD_DIR}/npt.cpt" ]; then
    cp "${MD_DIR}/npt.cpt" "${SIM_DIR}/" 2>/dev/null || true
fi

echo "Input files copied successfully."
echo ""

# Already in simulation directory (changed above)

echo "Starting 100 ns MD with GTX 1080 Ti GPU..."
echo "Structure: ${STRUCTURE_NAME}"
echo "Expected: ~9 ns/day (11 days total)"
echo "Working directory: ${SIM_DIR}"
echo ""

# This configuration WORKS - we tested it!
# Use structure name in output filename
gmx mdrun -deffnm "${STRUCTURE_NAME}_100ns" -s md.tpr -ntmpi 1 -ntomp 4 -v

echo ""
echo "Simulation completed: $(date)"
echo "Results saved in: ${SIM_DIR}"
