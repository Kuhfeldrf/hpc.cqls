#!/bin/bash
#
# Docking Pipeline
# ================
# Orchestrates the complete peptide docking workflow:
#   1. Predict peptide structures (ESMFold API)
#   2. Global docking (Rosetta docking_protocol) - searches entire receptor surface
#   3. FlexPepDocking refinement (Rosetta) - refines best global poses
#   4. Compile scores and select best structures
#
# Usage:
#   ./run_docking_pipeline.sh <peptide_list.txt> [nstruct_refine] [nstruct_global]
#
# Arguments:
#   peptide_list.txt  - Text file with peptide sequences (one per line)
#   nstruct_refine    - Number of FlexPepDocking refinement structures per peptide (default: 10)
#   nstruct_global    - Number of global docking structures per peptide (default: 50)
#
# Example:
#   ./run_docking_pipeline.sh my_peptides.txt 10 50
#

set -e

# =============================================================================
# CONFIGURATION
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

PEPTIDE_LIST="$1"
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
echo "       DOCKING PIPELINE"
echo "=========================================="
echo "Started: $(date)"
echo ""
echo "Configuration:"
echo "  Peptide list: ${PEPTIDE_LIST}"
echo "  Global docking structures: ${NSTRUCT_GLOBAL}"
echo "  Refinement structures: ${NSTRUCT_REFINE}"
echo "  Work directory: ${WORK_DIR}"
echo "  Global dock output: ${GLOBAL_DOCK_DIR}"
echo "  Final output: ${OUTPUT_DIR}"
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
    echo "  Please keep only one receptor file in the receptor folder."
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

# Export variables for child scripts
export PEPTIDE_LIST_FILE="${PEPTIDE_LIST}"
export NSTRUCT
export RECEPTOR_NAME
export WORK_DIR
export PEPTIDE_DIR
export OUTPUT_DIR

# =============================================================================
# STEP 1: PREDICT PEPTIDE STRUCTURES
# =============================================================================

echo "=========================================="
echo "STEP 1: Predict Peptide Structures"
echo "=========================================="

# Check which peptides need prediction
PEPTIDES_TO_PREDICT=()
for pep in "${PEPTIDES[@]}"; do
    if [ ! -f "${PEPTIDE_DIR}/${pep}.pdb" ]; then
        PEPTIDES_TO_PREDICT+=("${pep}")
    else
        echo "  ✓ ${pep}.pdb already exists"
    fi
done

if [ ${#PEPTIDES_TO_PREDICT[@]} -gt 0 ]; then
    echo "Predicting ${#PEPTIDES_TO_PREDICT[@]} peptide structure(s)..."
    
    # Load Python module
    module load intel-python/24.0.0 2>/dev/null || module load python/3.10 2>/dev/null || true
    
    # Run prediction for each peptide
    for pep in "${PEPTIDES_TO_PREDICT[@]}"; do
        echo "  Predicting: ${pep}"
        python3 "${WORK_DIR}/predict_peptide_structure.py" --sequence "${pep}" --output-dir "${PEPTIDE_DIR}" || {
            echo "  ⚠ Warning: Failed to predict ${pep}"
        }
    done
else
    echo "All peptide structures already exist."
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
    echo "Please provide these structures manually or check prediction errors."
    exit 1
fi

echo "✓ Step 1 complete: All peptide structures available"
echo ""

# =============================================================================
# STEP 2: GLOBAL DOCKING (Rosetta docking_protocol)
# =============================================================================

echo "=========================================="
echo "STEP 2: Global Docking"
echo "=========================================="
echo ""
echo "This step searches the ENTIRE receptor surface for binding poses."
echo "Using Rosetta docking_protocol with randomization."
echo ""

# Calculate array size
ARRAY_MAX=$((NUM_PEPTIDES - 1))

echo "Submitting global docking job array (0-${ARRAY_MAX})..."
echo "  ${NUM_PEPTIDES} peptides × ${NSTRUCT_GLOBAL} global structures each"

# Submit global docking job
GLOBAL_DOCK_JOB=$(sbatch --parsable \
    --array=0-${ARRAY_MAX} \
    --export=ALL,PEPTIDE_LIST_FILE="${PEPTIDE_LIST}",PEPTIDE_DIR="${PEPTIDE_DIR}",OUTPUT_DIR="${GLOBAL_DOCK_DIR}",NSTRUCT_GLOBAL="${NSTRUCT_GLOBAL}" \
    "${WORK_DIR}/global_dock_parallel.sh")

echo "  Global docking job submitted: ${GLOBAL_DOCK_JOB}"
echo "  Waiting for global docking to complete..."
echo "  (This may take 2-4 hours per peptide)"

# Wait for global docking jobs to complete
squeue -j ${GLOBAL_DOCK_JOB} &>/dev/null && {
    while squeue -j ${GLOBAL_DOCK_JOB} &>/dev/null 2>&1; do
        RUNNING=$(squeue -j ${GLOBAL_DOCK_JOB} -h 2>/dev/null | wc -l)
        if [ ${RUNNING} -eq 0 ]; then
            break
        fi
        echo "    Jobs remaining: ${RUNNING}"
        sleep 60
    done
}

echo "✓ Step 2 complete: Global docking finished"
echo ""

# =============================================================================
# STEP 3: FLEXPEP DOCKING REFINEMENT
# =============================================================================

echo "=========================================="
echo "STEP 3: FlexPepDocking Refinement"
echo "=========================================="
echo ""
echo "Refining best poses from global docking with full peptide flexibility."
echo ""

# Calculate total jobs: num_peptides * nstruct_refine
TOTAL_JOBS=$((NUM_PEPTIDES * NSTRUCT_REFINE))
DOCK_ARRAY_MAX=$((TOTAL_JOBS - 1))

echo "Submitting refinement job array (0-${DOCK_ARRAY_MAX})..."
echo "  ${NUM_PEPTIDES} peptides × ${NSTRUCT_REFINE} refinement structures = ${TOTAL_JOBS} jobs"

# Submit FlexPepDocking job - uses global dock output as input
DOCK_JOB=$(sbatch --parsable \
    --array=0-${DOCK_ARRAY_MAX} \
    --export=ALL,PEPTIDE_LIST_FILE="${PEPTIDE_LIST}",NSTRUCT="${NSTRUCT_REFINE}",PEPTIDE_DIR="${PEPTIDE_DIR}",OUTPUT_DIR="${OUTPUT_DIR}",GLOBAL_DOCK_DIR="${GLOBAL_DOCK_DIR}" \
    "${WORK_DIR}/flexpep_dock_parallel.sh")

echo "  Refinement job submitted: ${DOCK_JOB}"
echo "  Waiting for refinement to complete..."

# Wait for refinement jobs to complete
squeue -j ${DOCK_JOB} &>/dev/null && {
    while squeue -j ${DOCK_JOB} &>/dev/null 2>&1; do
        RUNNING=$(squeue -j ${DOCK_JOB} -h 2>/dev/null | wc -l)
        if [ ${RUNNING} -eq 0 ]; then
            break
        fi
        echo "    Jobs remaining: ${RUNNING}"
        sleep 60
    done
}

echo "✓ Step 3 complete: Refinement finished"
echo ""

# =============================================================================
# STEP 4: COMPILE SCORES AND SELECT BEST STRUCTURES
# =============================================================================

echo "=========================================="
echo "STEP 4: Compile Scores & Select Best"
echo "=========================================="

SUMMARY_FILE="${OUTPUT_DIR}/docking_summary_$(date +%Y%m%d_%H%M%S).csv"
BEST_STRUCTURES_DIR="${MD_DIR}"

# Create summary header
echo "peptide,receptor,structure,total_score,reweighted_sc,pep_sc,I_sc,description" > "${SUMMARY_FILE}"

echo "Compiling scores from all peptides..."

# Process each peptide
for pep in "${PEPTIDES[@]}"; do
    COMPLEX_NAME="${pep}_${RECEPTOR_NAME}"
    
    # Find all score files for this peptide
    SCORE_FILES=($(ls "${OUTPUT_DIR}/${COMPLEX_NAME}_scores"*.sc 2>/dev/null || true))
    
    if [ ${#SCORE_FILES[@]} -eq 0 ]; then
        echo "  ⚠ No score files found for ${pep}"
        continue
    fi
    
    echo "  Processing ${pep}..."
    
    # Combine all score files for this peptide
    COMBINED_SCORES="${OUTPUT_DIR}/${COMPLEX_NAME}_all_scores.txt"
    
    # Extract scores from all files (skip header lines)
    > "${COMBINED_SCORES}"
    for sf in "${SCORE_FILES[@]}"; do
        awk 'NR > 2 && NF > 1 {print}' "${sf}" >> "${COMBINED_SCORES}"
    done
    
    # Parse and add to summary
    if [ -s "${COMBINED_SCORES}" ]; then
        while read -r line; do
            # Extract key columns (adjust based on actual score file format)
            # Typical columns: SCORE: total_score reweighted_sc pep_sc I_sc ... description
            total=$(echo "${line}" | awk '{print $2}')
            reweighted=$(echo "${line}" | awk '{print $3}')
            pep_sc=$(echo "${line}" | awk '{print $4}')
            i_sc=$(echo "${line}" | awk '{print $5}')
            desc=$(echo "${line}" | awk '{print $NF}')
            
            echo "${pep},${RECEPTOR_NAME},${desc},${total},${reweighted},${pep_sc},${i_sc},${desc}" >> "${SUMMARY_FILE}"
        done < "${COMBINED_SCORES}"
        
        # Find best structure for this peptide (lowest total_score)
        BEST_LINE=$(sort -k2 -n "${COMBINED_SCORES}" | head -1)
        BEST_SCORE=$(echo "${BEST_LINE}" | awk '{print $2}')
        BEST_DESC=$(echo "${BEST_LINE}" | awk '{print $NF}')
        
        echo "    Best score: ${BEST_SCORE} (${BEST_DESC})"
        
        # Find and copy the best structure file
        BEST_FILE=$(ls "${OUTPUT_DIR}/"*"${BEST_DESC}"*.pdb 2>/dev/null | head -1)
        
        if [ -f "${BEST_FILE}" ]; then
            DEST_FILE="${BEST_STRUCTURES_DIR}/${pep}_${RECEPTOR_NAME}.pdb"
            cp "${BEST_FILE}" "${DEST_FILE}"
            echo "    ✓ Copied to: ${DEST_FILE}"
        else
            echo "    ⚠ Could not find PDB file for ${BEST_DESC}"
        fi
    fi
done

echo ""
echo "Summary file: ${SUMMARY_FILE}"

# Generate ranking across all peptides
echo ""
echo "=========================================="
echo "OVERALL RANKINGS (by total_score)"
echo "=========================================="

# Sort summary by total_score (column 4) and display top entries
echo ""
echo "Top structures across all peptides:"
tail -n +2 "${SUMMARY_FILE}" | sort -t',' -k4 -n | head -20 | \
    awk -F',' '{printf "  %-20s %10s  %s\n", $1, $4, $8}'

echo ""

# =============================================================================
# COMPLETION
# =============================================================================

echo "=========================================="
echo "       PIPELINE COMPLETE"
echo "=========================================="
echo "Completed: $(date)"
echo ""
echo "Results:"
echo "  Summary file: ${SUMMARY_FILE}"
echo "  Best structures copied to: ${BEST_STRUCTURES_DIR}/"
echo ""
echo "Best structures for MD simulation:"
for pep in "${PEPTIDES[@]}"; do
    DEST_FILE="${BEST_STRUCTURES_DIR}/${pep}_${RECEPTOR_NAME}.pdb"
    if [ -f "${DEST_FILE}" ]; then
        echo "  ✓ ${DEST_FILE}"
    else
        echo "  ✗ ${pep}_${RECEPTOR_NAME}.pdb (not found)"
    fi
done
echo ""
echo "To run MD simulations, use:"
echo "  cd ${MD_DIR}"
echo "  ./run_full_pipeline.sh ${PEPTIDES[*]/%/_${RECEPTOR_NAME}}"
echo ""

