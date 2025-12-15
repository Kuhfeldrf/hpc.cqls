#!/bin/bash
#SBATCH --job-name=flexpep_parallel
#SBATCH --output=/home/kuhfeldr/hpc.cqls/docking/logs/flexpep_par_%A_%a.out
#SBATCH --error=/home/kuhfeldr/hpc.cqls/docking/logs/flexpep_par_%A_%a.err
#SBATCH --partition=normal
#SBATCH --time=8:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G

# FlexPepDocking - PARALLEL REFINEMENT VERSION
# =============================================
# Refines global docking poses with full peptide flexibility.
# Each job produces 1 structure for maximum parallelism.
# Uses best global docking pose as input.
#
# Environment variables (set by run_docking_pipeline.sh):
#   PEPTIDE_LIST_FILE - Path to peptide list file (optional)
#   PEPTIDE_DIR       - Directory containing peptide PDB files
#   OUTPUT_DIR        - Output directory for refined structures
#   GLOBAL_DOCK_DIR   - Directory containing global docking output
#   NSTRUCT           - Number of structures per peptide

set -e

module load apptainer/1.4.1-gcc-13.4.0

# Paths - use environment variables if set, otherwise use defaults
WORK_DIR="/home/kuhfeldr/hpc.cqls/docking"
RECEPTOR_DIR="${WORK_DIR}/receptor"
LIGAND_DIR="${PEPTIDE_DIR:-${WORK_DIR}/peptides}"
OUTPUT_DIR="${OUTPUT_DIR:-${WORK_DIR}/docked_complexes}"
GLOBAL_DOCK_DIR="${GLOBAL_DOCK_DIR:-${WORK_DIR}/global_docked}"
ROSETTA_SIF="/home/kuhfeldr/rosetta.sif"

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

# Configuration - use environment variable if set
STRUCTURES_PER_PEPTIDE="${NSTRUCT:-10}"

mkdir -p "${OUTPUT_DIR}"
mkdir -p "${WORK_DIR}/logs"
cd "${WORK_DIR}"

# Get list of ligands - either from peptide list file or directory
if [ -n "${PEPTIDE_LIST_FILE}" ] && [ -f "${PEPTIDE_LIST_FILE}" ]; then
    # Read peptides from list file
    LIGANDS=()
    while IFS= read -r line || [ -n "$line" ]; do
        line=$(echo "${line}" | tr -d '\r' | xargs)
        if [ -n "${line}" ] && [[ ! "${line}" =~ ^# ]]; then
            pdb_file="${LIGAND_DIR}/${line}.pdb"
            if [ -f "${pdb_file}" ]; then
                LIGANDS+=("${pdb_file}")
            fi
        fi
    done < "${PEPTIDE_LIST_FILE}"
else
    # Get list from directory (sorted for reproducibility)
    LIGANDS=($(ls ${LIGAND_DIR}/*.pdb 2>/dev/null | sort))
fi

NUM_PEPTIDES=${#LIGANDS[@]}

if [ ${NUM_PEPTIDES} -eq 0 ]; then
    echo "ERROR: No ligand PDB files found"
    exit 1
fi

# Calculate which peptide and structure number this job handles
# Array ID 0-(nstruct-1) = peptide 0
# Array ID nstruct-(2*nstruct-1) = peptide 1
# etc.
PEPTIDE_INDEX=$((SLURM_ARRAY_TASK_ID / STRUCTURES_PER_PEPTIDE))
STRUCT_NUM=$((SLURM_ARRAY_TASK_ID % STRUCTURES_PER_PEPTIDE))

if [ ${PEPTIDE_INDEX} -ge ${NUM_PEPTIDES} ]; then
    echo "No work for this task (peptide index ${PEPTIDE_INDEX} >= ${NUM_PEPTIDES})"
    exit 0
fi

LIGAND="${LIGANDS[${PEPTIDE_INDEX}]}"
LIGAND_NAME=$(basename "${LIGAND}" .pdb)

echo "=========================================="
echo "FlexPepDocking REFINEMENT"
echo "=========================================="
echo "Started: $(date)"
echo "Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Peptide: ${LIGAND_NAME} (index ${PEPTIDE_INDEX})"
echo "Refinement structure: ${STRUCT_NUM} of ${STRUCTURES_PER_PEPTIDE}"
echo "Global dock dir: ${GLOBAL_DOCK_DIR}"
echo ""

# Function to run Rosetta
run_rosetta() {
    apptainer exec --bind "${WORK_DIR}:${WORK_DIR}" "${ROSETTA_SIF}" "$@"
}

# Get receptor chains
RECEPTOR_CHAINS=$(grep "^ATOM" "${RECEPTOR}" | awk '{print substr($0,22,1)}' | sort -u | tr -d '\n')

# Setup complex - look for global docking output
COMPLEX_NAME="${LIGAND_NAME}_${RECEPTOR_NAME}"

# Priority order for input:
# 1. Best global docking pose (preferred)
# 2. Top 1 global pose
# 3. Input complex from global docking
GLOBAL_BEST="${GLOBAL_DOCK_DIR}/${COMPLEX_NAME}_best_global.pdb"
GLOBAL_TOP1="${GLOBAL_DOCK_DIR}/${COMPLEX_NAME}_top1_global.pdb"
GLOBAL_INPUT="${GLOBAL_DOCK_DIR}/${COMPLEX_NAME}_input.pdb"

# Determine which input file to use
INPUT_FOR_DOCKING=""

if [ -f "${GLOBAL_BEST}" ]; then
    echo "✓ Found best global docking pose: ${GLOBAL_BEST}"
    INPUT_FOR_DOCKING="${GLOBAL_BEST}"
elif [ -f "${GLOBAL_TOP1}" ]; then
    echo "✓ Found top 1 global docking pose: ${GLOBAL_TOP1}"
    INPUT_FOR_DOCKING="${GLOBAL_TOP1}"
elif [ -f "${GLOBAL_INPUT}" ]; then
    echo "⚠ No global dock output, using input complex: ${GLOBAL_INPUT}"
    INPUT_FOR_DOCKING="${GLOBAL_INPUT}"
else
    echo "ERROR: No input structure found!"
    echo "  Looked for: ${GLOBAL_BEST}"
    echo "  And: ${GLOBAL_TOP1}"
    echo "  And: ${GLOBAL_INPUT}"
    echo ""
    echo "Please run global docking first:"
    echo "  sbatch global_dock_parallel.sh"
    exit 1
fi

echo "Using input: ${INPUT_FOR_DOCKING}"

echo ""
echo "Running FlexPepDocking refinement structure ${STRUCT_NUM}..."
echo "  Input: Global docking pose"
echo "  Output: Refined pose with full peptide flexibility"
echo ""

# Run single structure with unique random seed based on task ID
run_rosetta FlexPepDocking.default.linuxgccrelease \
    -s "${INPUT_FOR_DOCKING}" \
    -flexPepDocking:receptor_chain "${RECEPTOR_CHAINS}" \
    -flexPepDocking:peptide_chain "P" \
    -nstruct 1 \
    -out:path:pdb "${OUTPUT_DIR}" \
    -out:file:scorefile "${OUTPUT_DIR}/${COMPLEX_NAME}_scores_${STRUCT_NUM}.sc" \
    -out:prefix "${COMPLEX_NAME}_s${STRUCT_NUM}_" \
    -out:suffix "_docked" \
    -lowres_preoptimize \
    -pep_refine \
    -flexpep_score_only false \
    -ex1 \
    -ex2aro \
    -use_input_sc \
    -unboundrot "${LIGAND}" \
    -score:weights ref2015 \
    -min_receptor_bb \
    -renumber_pdb false \
    -run:constant_seed \
    -run:jran $((SLURM_ARRAY_TASK_ID * 12345 + 98765)) \
    -overwrite 2>&1

echo ""
echo "=========================================="
echo "Structure ${STRUCT_NUM} complete for ${LIGAND_NAME}"
echo "=========================================="
echo "Completed: $(date)"

