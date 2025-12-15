#!/bin/bash
#SBATCH --job-name=dock_scoring
#SBATCH --output=/home/kuhfeldr/hpc.cqls/docking/logs/scoring_%j.out
#SBATCH --error=/home/kuhfeldr/hpc.cqls/docking/logs/scoring_%j.err
#SBATCH --partition=normal
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

# Compile Scores and Select Best Structures
# =========================================
# Final step of the docking pipeline.
# Compiles all FlexPepDocking scores and copies best structures to MD directory.
#
# Environment variables:
#   PEPTIDE_LIST_FILE - Path to peptide list file
#   OUTPUT_DIR        - Directory containing FlexPepDocking output
#   MD_DIR            - Directory to copy best structures to
#   RECEPTOR_NAME     - Name of the receptor (without .pdb extension)

set -e

echo "=========================================="
echo "   COMPILE SCORES & SELECT BEST"
echo "=========================================="
echo "Started: $(date)"
echo ""

# Configuration - use environment variables with absolute path defaults
WORK_DIR="/home/kuhfeldr/hpc.cqls/docking"
OUTPUT_DIR="${OUTPUT_DIR:-/home/kuhfeldr/hpc.cqls/docking/docked_complexes}"
MD_DIR="${MD_DIR:-/home/kuhfeldr/hpc.cqls/md}"
GLOBAL_DOCK_DIR="${GLOBAL_DOCK_DIR:-/home/kuhfeldr/hpc.cqls/docking/global_docked}"

# Get receptor name if not provided
if [ -z "${RECEPTOR_NAME}" ]; then
    RECEPTOR_FILES=($(ls ${WORK_DIR}/receptor/*.pdb 2>/dev/null))
    if [ ${#RECEPTOR_FILES[@]} -gt 0 ]; then
        RECEPTOR_NAME=$(basename "${RECEPTOR_FILES[0]}" .pdb)
    else
        RECEPTOR_NAME="receptor"
    fi
fi

echo "Configuration:"
echo "  Output directory: ${OUTPUT_DIR}"
echo "  MD directory: ${MD_DIR}"
echo "  Receptor name: ${RECEPTOR_NAME}"
echo "  Peptide list: ${PEPTIDE_LIST_FILE}"
echo ""

# Read peptide sequences
PEPTIDES=()
if [ -n "${PEPTIDE_LIST_FILE}" ] && [ -f "${PEPTIDE_LIST_FILE}" ]; then
    while IFS= read -r line || [ -n "$line" ]; do
        line=$(echo "${line}" | tr -d '\r' | xargs)
        if [ -n "${line}" ] && [[ ! "${line}" =~ ^# ]]; then
            PEPTIDES+=("${line}")
        fi
    done < "${PEPTIDE_LIST_FILE}"
else
    echo "ERROR: Peptide list file not found: ${PEPTIDE_LIST_FILE}"
    exit 1
fi

NUM_PEPTIDES=${#PEPTIDES[@]}
echo "Processing ${NUM_PEPTIDES} peptides..."
echo ""

# Create directories
mkdir -p "${MD_DIR}"

# Create summary file
SUMMARY_FILE="${OUTPUT_DIR}/docking_summary_$(date +%Y%m%d_%H%M%S).csv"
echo "peptide,receptor,structure,total_score,reweighted_sc,pep_sc,I_sc,source_file" > "${SUMMARY_FILE}"

echo "=========================================="
echo "Compiling FlexPepDocking Scores"
echo "=========================================="

PEPTIDES_PROCESSED=0
PEPTIDES_FOUND=0

for pep in "${PEPTIDES[@]}"; do
    COMPLEX_NAME="${pep}_${RECEPTOR_NAME}"
    
    echo ""
    echo "Processing: ${pep}"
    
    # Find all score files for this peptide from FlexPepDocking
    SCORE_FILES=($(ls "${OUTPUT_DIR}/${COMPLEX_NAME}_scores"*.sc 2>/dev/null || true))
    
    if [ ${#SCORE_FILES[@]} -eq 0 ]; then
        echo "  ⚠ No FlexPepDocking score files found"
        echo "  Looking for global docking scores as fallback..."
        
        # Fallback to global docking best if no FlexPep results
        if [ -f "${GLOBAL_DOCK_DIR}/${COMPLEX_NAME}_best_global.pdb" ]; then
            echo "  ✓ Found global docking best structure"
            DEST_FILE="${MD_DIR}/${pep}_${RECEPTOR_NAME}.pdb"
            cp "${GLOBAL_DOCK_DIR}/${COMPLEX_NAME}_best_global.pdb" "${DEST_FILE}"
            echo "  ✓ Copied global dock best to: ${DEST_FILE}"
            ((PEPTIDES_FOUND++)) || true
        else
            echo "  ✗ No structures found for ${pep}"
        fi
        continue
    fi
    
    echo "  Found ${#SCORE_FILES[@]} score file(s)"
    
    # Combine all score files for this peptide
    COMBINED_SCORES="${OUTPUT_DIR}/${COMPLEX_NAME}_all_scores.txt"
    > "${COMBINED_SCORES}"
    
    for sf in "${SCORE_FILES[@]}"; do
        awk 'NR > 2 && NF > 1 {print}' "${sf}" >> "${COMBINED_SCORES}"
    done
    
    SCORE_COUNT=$(wc -l < "${COMBINED_SCORES}")
    echo "  Total scores: ${SCORE_COUNT}"
    
    if [ "${SCORE_COUNT}" -eq 0 ]; then
        echo "  ⚠ No scores found in score files"
        continue
    fi
    
    # Parse and add to summary
    while read -r line; do
        total=$(echo "${line}" | awk '{print $2}')
        reweighted=$(echo "${line}" | awk '{print $3}')
        pep_sc=$(echo "${line}" | awk '{print $4}')
        i_sc=$(echo "${line}" | awk '{print $5}')
        desc=$(echo "${line}" | awk '{print $NF}')
        
        echo "${pep},${RECEPTOR_NAME},${desc},${total},${reweighted},${pep_sc},${i_sc},${desc}" >> "${SUMMARY_FILE}"
    done < "${COMBINED_SCORES}"
    
    # Find best structure (lowest total_score)
    BEST_LINE=$(sort -k2 -n "${COMBINED_SCORES}" | head -1)
    BEST_SCORE=$(echo "${BEST_LINE}" | awk '{print $2}')
    BEST_DESC=$(echo "${BEST_LINE}" | awk '{print $NF}')
    
    echo "  Best score: ${BEST_SCORE}"
    echo "  Best structure: ${BEST_DESC}"
    
    # Find and copy the best structure file
    BEST_FILE=""
    
    # Try different naming patterns
    for pattern in "${OUTPUT_DIR}/${BEST_DESC}.pdb" \
                   "${OUTPUT_DIR}/"*"${BEST_DESC}"*.pdb \
                   "${OUTPUT_DIR}/${COMPLEX_NAME}"*"${BEST_DESC}"*.pdb; do
        if ls ${pattern} &>/dev/null; then
            BEST_FILE=$(ls ${pattern} 2>/dev/null | head -1)
            break
        fi
    done
    
    if [ -f "${BEST_FILE}" ]; then
        DEST_FILE="${MD_DIR}/${pep}_${RECEPTOR_NAME}.pdb"
        cp "${BEST_FILE}" "${DEST_FILE}"
        echo "  ✓ Copied to: ${DEST_FILE}"
        ((PEPTIDES_FOUND++)) || true
    else
        echo "  ⚠ Could not find PDB file for ${BEST_DESC}"
        
        # Fallback to global docking best
        if [ -f "${GLOBAL_DOCK_DIR}/${COMPLEX_NAME}_best_global.pdb" ]; then
            echo "  Using global docking best as fallback..."
            DEST_FILE="${MD_DIR}/${pep}_${RECEPTOR_NAME}.pdb"
            cp "${GLOBAL_DOCK_DIR}/${COMPLEX_NAME}_best_global.pdb" "${DEST_FILE}"
            echo "  ✓ Copied global dock best to: ${DEST_FILE}"
            ((PEPTIDES_FOUND++)) || true
        fi
    fi
    
    ((PEPTIDES_PROCESSED++)) || true
done

echo ""
echo "=========================================="
echo "Summary"
echo "=========================================="
echo "  Score file: ${SUMMARY_FILE}"
echo "  Peptides with FlexPep scores: ${PEPTIDES_PROCESSED}"
echo "  Structures copied to MD: ${PEPTIDES_FOUND}"
echo ""

# Generate ranking
echo "=========================================="
echo "OVERALL RANKINGS (by total_score)"
echo "=========================================="
echo ""
echo "Top structures across all peptides:"

tail -n +2 "${SUMMARY_FILE}" | sort -t',' -k4 -n | head -20 | \
    awk -F',' '{printf "  %-20s %12s  %s\n", $1, $4, $3}'

echo ""

# Show what was copied to MD directory
echo "=========================================="
echo "Best Structures for MD Simulation"
echo "=========================================="
echo "Location: ${MD_DIR}/"
echo ""

for pep in "${PEPTIDES[@]}"; do
    DEST_FILE="${MD_DIR}/${pep}_${RECEPTOR_NAME}.pdb"
    if [ -f "${DEST_FILE}" ]; then
        SIZE=$(ls -lh "${DEST_FILE}" | awk '{print $5}')
        echo "  ✓ ${pep}_${RECEPTOR_NAME}.pdb (${SIZE})"
    else
        echo "  ✗ ${pep}_${RECEPTOR_NAME}.pdb (not found)"
    fi
done

echo ""
echo "=========================================="
echo "       PIPELINE COMPLETE"
echo "=========================================="
echo "Completed: $(date)"
echo ""
echo "To run MD simulations:"
echo "  cd ${MD_DIR}"
echo "  # Submit MD jobs for each complex"
echo ""

