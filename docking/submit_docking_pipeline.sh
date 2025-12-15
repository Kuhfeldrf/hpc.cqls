#!/bin/bash
#
# Docking Pipeline - SLURM SUBMISSION WRAPPER
# ============================================
# Submits all pipeline steps as Slurm jobs with proper dependencies.
# This script submits jobs and exits - NO waiting required.
# Job dependencies ensure proper execution order.
# Automatically skips steps if output already exists.
#
# Usage:
#   ./submit_docking_pipeline.sh <peptide_list.txt> [nstruct_refine] [nstruct_global]
#
# Example:
#   ./submit_docking_pipeline.sh my_peptides.txt 10 50
#

set -e

# =============================================================================
# CONFIGURATION - ALL ABSOLUTE PATHS
# =============================================================================

WORK_DIR="/home/kuhfeldr/hpc.cqls/docking"
MD_DIR="/home/kuhfeldr/hpc.cqls/md"
RECEPTOR_DIR="${WORK_DIR}/receptor"
PEPTIDE_DIR="${WORK_DIR}/peptides"
GLOBAL_DOCK_DIR="${WORK_DIR}/global_docked"
OUTPUT_DIR="${WORK_DIR}/docked_complexes"
LOG_DIR="${WORK_DIR}/logs"

# Default values
DEFAULT_NSTRUCT_REFINE=10
DEFAULT_NSTRUCT_GLOBAL=50

# =============================================================================
# ARGUMENT PARSING
# =============================================================================

print_usage() {
    echo "Usage: $0 <peptide_list.txt> [nstruct_refine] [nstruct_global]"
    echo ""
    echo "Arguments:"
    echo "  peptide_list.txt  - Text file with peptide sequences (one per line)"
    echo "  nstruct_refine    - FlexPepDocking structures per peptide (default: ${DEFAULT_NSTRUCT_REFINE})"
    echo "  nstruct_global    - Global docking structures per peptide (default: ${DEFAULT_NSTRUCT_GLOBAL})"
    echo ""
    echo "Example:"
    echo "  $0 my_peptides.txt 10 50"
    echo ""
}

if [ $# -lt 1 ]; then
    print_usage
    exit 1
fi

# Get absolute path to peptide list
PEPTIDE_LIST="$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
NSTRUCT_REFINE="${2:-$DEFAULT_NSTRUCT_REFINE}"
NSTRUCT_GLOBAL="${3:-$DEFAULT_NSTRUCT_GLOBAL}"

# Validate inputs
if [ ! -f "${PEPTIDE_LIST}" ]; then
    echo "ERROR: Peptide list file not found: ${PEPTIDE_LIST}"
    exit 1
fi

if ! [[ "${NSTRUCT_REFINE}" =~ ^[0-9]+$ ]] || [ "${NSTRUCT_REFINE}" -lt 1 ]; then
    echo "ERROR: nstruct_refine must be a positive integer"
    exit 1
fi

if ! [[ "${NSTRUCT_GLOBAL}" =~ ^[0-9]+$ ]] || [ "${NSTRUCT_GLOBAL}" -lt 1 ]; then
    echo "ERROR: nstruct_global must be a positive integer"
    exit 1
fi

# =============================================================================
# INITIALIZATION
# =============================================================================

echo "=========================================="
echo "   DOCKING PIPELINE - SLURM SUBMISSION"
echo "=========================================="
echo "Submission Time: $(date)"
echo ""
echo "Configuration:"
echo "  Peptide list: ${PEPTIDE_LIST}"
echo "  Global docking structures: ${NSTRUCT_GLOBAL}"
echo "  Refinement structures: ${NSTRUCT_REFINE}"
echo "  Work directory: ${WORK_DIR}"
echo "  MD output directory: ${MD_DIR}"
echo ""

# Create directories
mkdir -p "${PEPTIDE_DIR}"
mkdir -p "${GLOBAL_DOCK_DIR}"
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${LOG_DIR}"
mkdir -p "${MD_DIR}"

# Validate receptor
RECEPTOR_FILES=($(ls ${RECEPTOR_DIR}/*.pdb 2>/dev/null))
if [ ${#RECEPTOR_FILES[@]} -eq 0 ]; then
    echo "ERROR: No receptor PDB file found in ${RECEPTOR_DIR}/"
    exit 1
elif [ ${#RECEPTOR_FILES[@]} -gt 1 ]; then
    echo "ERROR: Multiple receptor PDB files found in ${RECEPTOR_DIR}/"
    echo "  Found: ${RECEPTOR_FILES[*]}"
    exit 1
fi
RECEPTOR="${RECEPTOR_FILES[0]}"
RECEPTOR_NAME=$(basename "${RECEPTOR}" .pdb)
echo "Receptor: ${RECEPTOR_NAME}"

# Read peptide sequences
PEPTIDES=()
while IFS= read -r line || [ -n "$line" ]; do
    line=$(echo "${line}" | tr -d '\r' | xargs)
    if [ -n "${line}" ] && [[ ! "${line}" =~ ^# ]]; then
        PEPTIDES+=("${line}")
    fi
done < "${PEPTIDE_LIST}"

NUM_PEPTIDES=${#PEPTIDES[@]}
if [ ${NUM_PEPTIDES} -eq 0 ]; then
    echo "ERROR: No peptide sequences found in ${PEPTIDE_LIST}"
    exit 1
fi

echo "Peptides to process: ${NUM_PEPTIDES}"
for pep in "${PEPTIDES[@]}"; do
    echo "  - ${pep}"
done
echo ""

# Track submitted jobs for scancel command
SUBMITTED_JOBS=()

# =============================================================================
# STEP 1: PREDICT PEPTIDE STRUCTURES (runs inline - fast)
# =============================================================================

echo "=========================================="
echo "STEP 1: Check/Predict Peptide Structures"
echo "=========================================="

# Check which peptides need prediction
PEPTIDES_TO_PREDICT=()
PEPTIDES_EXISTING=0
for pep in "${PEPTIDES[@]}"; do
    if [ ! -f "${PEPTIDE_DIR}/${pep}.pdb" ]; then
        PEPTIDES_TO_PREDICT+=("${pep}")
        echo "  ○ ${pep}.pdb: needs prediction"
    else
        echo "  ✓ ${pep}.pdb: already exists"
        ((PEPTIDES_EXISTING++)) || true
    fi
done

if [ ${PEPTIDES_EXISTING} -eq ${NUM_PEPTIDES} ]; then
    echo ""
    echo "  → SKIPPING peptide prediction: All ${NUM_PEPTIDES} structures already exist."
elif [ ${#PEPTIDES_TO_PREDICT[@]} -gt 0 ]; then
    echo ""
    echo "Predicting ${#PEPTIDES_TO_PREDICT[@]} peptide structure(s)..."
    
    # Load Python module
    module load intel-python/24.0.0 2>/dev/null || module load python/3.10 2>/dev/null || true
    
    for pep in "${PEPTIDES_TO_PREDICT[@]}"; do
        echo "  Predicting: ${pep}"
        python3 "${WORK_DIR}/predict_peptide_structure.py" --sequence "${pep}" --output-dir "${PEPTIDE_DIR}" || {
            echo "  ⚠ Warning: Failed to predict ${pep}"
        }
    done
fi

# Verify all peptides have structures
MISSING_PEPTIDES=()
for pep in "${PEPTIDES[@]}"; do
    if [ ! -f "${PEPTIDE_DIR}/${pep}.pdb" ]; then
        MISSING_PEPTIDES+=("${pep}")
    fi
done

if [ ${#MISSING_PEPTIDES[@]} -gt 0 ]; then
    echo ""
    echo "ERROR: Missing peptide structures:"
    for pep in "${MISSING_PEPTIDES[@]}"; do
        echo "  - ${pep}"
    done
    exit 1
fi

echo "✓ Step 1 complete: All peptide structures available"
echo ""

# =============================================================================
# CHECK EXISTING OUTPUT - DETERMINE WHAT NEEDS TO RUN
# =============================================================================

echo "=========================================="
echo "Checking Existing Output"
echo "=========================================="

# Check global docking output
echo ""
echo "Global Docking:"
GLOBAL_COMPLETE=0
for pep in "${PEPTIDES[@]}"; do
    COMPLEX_NAME="${pep}_${RECEPTOR_NAME}"
    if [ -f "${GLOBAL_DOCK_DIR}/${COMPLEX_NAME}_best_global.pdb" ]; then
        echo "  ✓ ${pep}: global docking complete"
        ((GLOBAL_COMPLETE++)) || true
    else
        echo "  ○ ${pep}: global docking needed"
    fi
done

# Check FlexPepDocking output
echo ""
echo "FlexPepDocking:"
FLEXPEP_COMPLETE=0
for pep in "${PEPTIDES[@]}"; do
    COMPLEX_NAME="${pep}_${RECEPTOR_NAME}"
    SCORE_COUNT=$(ls "${OUTPUT_DIR}/${COMPLEX_NAME}_scores"*.sc 2>/dev/null | wc -l)
    if [ "${SCORE_COUNT}" -ge "${NSTRUCT_REFINE}" ]; then
        echo "  ✓ ${pep}: FlexPepDocking complete (${SCORE_COUNT} score files)"
        ((FLEXPEP_COMPLETE++)) || true
    else
        echo "  ○ ${pep}: FlexPepDocking needed (found ${SCORE_COUNT}/${NSTRUCT_REFINE})"
    fi
done

# Check MD directory
echo ""
echo "Final Output (MD directory):"
MD_COMPLETE=0
for pep in "${PEPTIDES[@]}"; do
    DEST_FILE="${MD_DIR}/${pep}_${RECEPTOR_NAME}.pdb"
    if [ -f "${DEST_FILE}" ]; then
        SIZE=$(ls -lh "${DEST_FILE}" | awk '{print $5}')
        echo "  ✓ ${pep}: exists (${SIZE})"
        ((MD_COMPLETE++)) || true
    else
        echo "  ○ ${pep}: needs to be copied"
    fi
done

echo ""

# =============================================================================
# STEP 2: SUBMIT GLOBAL DOCKING JOBS (if needed)
# =============================================================================

GLOBAL_DOCK_JOB=""
ARRAY_MAX=$((NUM_PEPTIDES - 1))

if [ ${GLOBAL_COMPLETE} -eq ${NUM_PEPTIDES} ]; then
    echo "=========================================="
    echo "STEP 2: Global Docking - SKIPPING"
    echo "=========================================="
    echo "  All ${NUM_PEPTIDES} peptides already have global docking output."
    echo ""
else
    echo "=========================================="
    echo "STEP 2: Submit Global Docking Jobs"
    echo "=========================================="
    
    echo "Submitting global docking job array (0-${ARRAY_MAX})..."
    echo "  ${NUM_PEPTIDES} peptides × ${NSTRUCT_GLOBAL} global structures each"
    
    GLOBAL_DOCK_JOB=$(sbatch --parsable \
        --array=0-${ARRAY_MAX} \
        --export=ALL,PEPTIDE_LIST_FILE="${PEPTIDE_LIST}",PEPTIDE_DIR="${PEPTIDE_DIR}",OUTPUT_DIR="${GLOBAL_DOCK_DIR}",NSTRUCT_GLOBAL="${NSTRUCT_GLOBAL}" \
        "${WORK_DIR}/global_dock_parallel.sh")
    
    echo "  ✓ Global docking job submitted: ${GLOBAL_DOCK_JOB}"
    SUBMITTED_JOBS+=("${GLOBAL_DOCK_JOB}")
    echo ""
fi

# =============================================================================
# STEP 3: SUBMIT FLEXPEP DOCKING JOBS (if needed)
# =============================================================================

FLEXPEP_JOB=""
TOTAL_JOBS=$((NUM_PEPTIDES * NSTRUCT_REFINE))
DOCK_ARRAY_MAX=$((TOTAL_JOBS - 1))

if [ ${FLEXPEP_COMPLETE} -eq ${NUM_PEPTIDES} ]; then
    echo "=========================================="
    echo "STEP 3: FlexPepDocking - SKIPPING"
    echo "=========================================="
    echo "  All ${NUM_PEPTIDES} peptides already have FlexPepDocking output."
    echo ""
else
    echo "=========================================="
    echo "STEP 3: Submit FlexPepDocking Jobs"
    echo "=========================================="
    
    echo "Submitting refinement job array (0-${DOCK_ARRAY_MAX})..."
    echo "  ${NUM_PEPTIDES} peptides × ${NSTRUCT_REFINE} refinement structures = ${TOTAL_JOBS} jobs"
    
    # Set dependency if global dock job was submitted
    DEPENDENCY_FLAG=""
    if [ -n "${GLOBAL_DOCK_JOB}" ]; then
        DEPENDENCY_FLAG="--dependency=afterok:${GLOBAL_DOCK_JOB}"
        echo "  Dependency: afterok:${GLOBAL_DOCK_JOB}"
    fi
    
    FLEXPEP_JOB=$(sbatch --parsable \
        --array=0-${DOCK_ARRAY_MAX} \
        ${DEPENDENCY_FLAG} \
        --export=ALL,PEPTIDE_LIST_FILE="${PEPTIDE_LIST}",NSTRUCT="${NSTRUCT_REFINE}",PEPTIDE_DIR="${PEPTIDE_DIR}",OUTPUT_DIR="${OUTPUT_DIR}",GLOBAL_DOCK_DIR="${GLOBAL_DOCK_DIR}" \
        "${WORK_DIR}/flexpep_dock_parallel.sh")
    
    echo "  ✓ FlexPepDocking job submitted: ${FLEXPEP_JOB}"
    SUBMITTED_JOBS+=("${FLEXPEP_JOB}")
    echo ""
fi

# =============================================================================
# STEP 4: SUBMIT FINAL SCORING JOB (if needed)
# =============================================================================

SCORING_JOB=""

if [ ${MD_COMPLETE} -eq ${NUM_PEPTIDES} ]; then
    echo "=========================================="
    echo "STEP 4: Scoring/Selection - SKIPPING"
    echo "=========================================="
    echo "  All ${NUM_PEPTIDES} structures already in ${MD_DIR}/"
    echo ""
else
    echo "=========================================="
    echo "STEP 4: Submit Final Scoring Job"
    echo "=========================================="
    
    echo "Submitting scoring/selection job..."
    
    # Set dependency based on what was submitted
    DEPENDENCY_FLAG=""
    if [ -n "${FLEXPEP_JOB}" ]; then
        DEPENDENCY_FLAG="--dependency=afterok:${FLEXPEP_JOB}"
        echo "  Dependency: afterok:${FLEXPEP_JOB}"
    elif [ -n "${GLOBAL_DOCK_JOB}" ]; then
        DEPENDENCY_FLAG="--dependency=afterok:${GLOBAL_DOCK_JOB}"
        echo "  Dependency: afterok:${GLOBAL_DOCK_JOB}"
    fi
    
    SCORING_JOB=$(sbatch --parsable \
        ${DEPENDENCY_FLAG} \
        --job-name=dock_scoring \
        --output="${LOG_DIR}/scoring_%j.out" \
        --error="${LOG_DIR}/scoring_%j.err" \
        --partition=normal \
        --time=1:00:00 \
        --nodes=1 \
        --ntasks=1 \
        --cpus-per-task=1 \
        --mem=4G \
        --export=ALL,PEPTIDE_LIST_FILE="${PEPTIDE_LIST}",OUTPUT_DIR="${OUTPUT_DIR}",MD_DIR="${MD_DIR}",RECEPTOR_NAME="${RECEPTOR_NAME}",GLOBAL_DOCK_DIR="${GLOBAL_DOCK_DIR}" \
        "${WORK_DIR}/compile_scores_and_select.sh")
    
    echo "  ✓ Scoring job submitted: ${SCORING_JOB}"
    SUBMITTED_JOBS+=("${SCORING_JOB}")
    echo ""
fi

# =============================================================================
# SUMMARY
# =============================================================================

echo "=========================================="
echo "       SUMMARY"
echo "=========================================="
echo ""

if [ ${#SUBMITTED_JOBS[@]} -eq 0 ]; then
    echo "No jobs submitted - all steps already complete!"
    echo ""
    echo "Final structures are in: ${MD_DIR}/"
    for pep in "${PEPTIDES[@]}"; do
        echo "  - ${pep}_${RECEPTOR_NAME}.pdb"
    done
else
    echo "Jobs Submitted:"
    STEP_NUM=1
    if [ -n "${GLOBAL_DOCK_JOB}" ]; then
        echo "  ${STEP_NUM}. Global Docking: ${GLOBAL_DOCK_JOB} (array 0-${ARRAY_MAX})"
        ((STEP_NUM++))
    fi
    if [ -n "${FLEXPEP_JOB}" ]; then
        DEP_MSG=""
        [ -n "${GLOBAL_DOCK_JOB}" ] && DEP_MSG=" (depends on ${GLOBAL_DOCK_JOB})"
        echo "  ${STEP_NUM}. FlexPepDocking: ${FLEXPEP_JOB} (array 0-${DOCK_ARRAY_MAX})${DEP_MSG}"
        ((STEP_NUM++))
    fi
    if [ -n "${SCORING_JOB}" ]; then
        DEP_MSG=""
        [ -n "${FLEXPEP_JOB}" ] && DEP_MSG=" (depends on ${FLEXPEP_JOB})"
        [ -z "${FLEXPEP_JOB}" ] && [ -n "${GLOBAL_DOCK_JOB}" ] && DEP_MSG=" (depends on ${GLOBAL_DOCK_JOB})"
        echo "  ${STEP_NUM}. Scoring/Select: ${SCORING_JOB}${DEP_MSG}"
    fi
    echo ""
    echo "Monitor progress with:"
    echo "  squeue -u \$USER"
    echo "  watch -n 30 'squeue -u \$USER'"
    echo ""
    echo "Check logs in:"
    echo "  ${LOG_DIR}/"
    echo ""
    echo "Final output will be in:"
    echo "  ${MD_DIR}/"
    echo ""
    
    # Generate scancel command
    ALL_JOBS=$(IFS=,; echo "${SUBMITTED_JOBS[*]}")
    echo "=========================================="
    echo "To CANCEL all submitted jobs, run:"
    echo "=========================================="
    echo ""
    echo "  scancel ${ALL_JOBS}"
    echo ""
fi

echo "You can safely close this terminal - jobs will continue running."
echo ""
