#!/bin/bash
#SBATCH --job-name=global_dock
#SBATCH --output=/home/kuhfeldr/hpc.cqls/docking/logs/global_dock_%A_%a.out
#SBATCH --error=/home/kuhfeldr/hpc.cqls/docking/logs/global_dock_%A_%a.err
#SBATCH --partition=normal
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G

# Global Peptide-Protein Docking using Rosetta docking_protocol
# ==============================================================
# Performs TRUE global docking - searches entire receptor surface
# for optimal peptide binding poses.
#
# This should be run BEFORE FlexPepDocking refinement.
#
# Environment variables (set by run_docking_pipeline.sh):
#   PEPTIDE_LIST_FILE - Path to peptide list file (optional)
#   PEPTIDE_DIR       - Directory containing peptide PDB files
#   OUTPUT_DIR        - Output directory for docked structures
#   NSTRUCT_GLOBAL    - Number of global docking structures per peptide (default: 50)

set -e

# Load apptainer module
module load apptainer/1.4.1-gcc-13.4.0

# Paths - use environment variables if set, otherwise use defaults
WORK_DIR="/home/kuhfeldr/hpc.cqls/docking"
RECEPTOR_DIR="${WORK_DIR}/receptor"
LIGAND_DIR="${PEPTIDE_DIR:-${WORK_DIR}/peptides}"
OUTPUT_DIR="${OUTPUT_DIR:-${WORK_DIR}/global_docked}"
ROSETTA_SIF="/home/kuhfeldr/rosetta.sif"

# Configuration - number of structures for global docking
NSTRUCT="${NSTRUCT_GLOBAL:-50}"

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

# Create output and log directories
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

if [ ${#LIGANDS[@]} -eq 0 ]; then
    echo "ERROR: No ligand PDB files found in ${LIGAND_DIR}/"
    exit 1
fi

# Select ligand based on SLURM_ARRAY_TASK_ID
if [ ${SLURM_ARRAY_TASK_ID} -ge ${#LIGANDS[@]} ]; then
    echo "ERROR: Array task ID ${SLURM_ARRAY_TASK_ID} exceeds number of ligands (${#LIGANDS[@]})"
    exit 1
fi

LIGAND="${LIGANDS[${SLURM_ARRAY_TASK_ID}]}"
LIGAND_NAME=$(basename "${LIGAND}" .pdb)

echo "=========================================="
echo "Global Peptide-Protein Docking"
echo "=========================================="
echo "Started: $(date)"
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Container: ${ROSETTA_SIF}"
echo "Receptor: ${RECEPTOR_NAME}"
echo "Peptide: ${LIGAND_NAME}"
echo "Structures to generate: ${NSTRUCT}"
echo ""

# Verify container exists
if [ ! -f "${ROSETTA_SIF}" ]; then
    echo "ERROR: Rosetta container not found at ${ROSETTA_SIF}"
    exit 1
fi

# Function to run Rosetta commands via container
run_rosetta() {
    apptainer exec --bind "${WORK_DIR}:${WORK_DIR}" "${ROSETTA_SIF}" "$@"
}

# Get receptor chain IDs
RECEPTOR_CHAINS=$(grep "^ATOM" "${RECEPTOR}" | awk '{print substr($0,22,1)}' | sort -u | tr -d '\n')
echo "Receptor chains found: ${RECEPTOR_CHAINS}"
echo ""

# Setup file names
COMPLEX_NAME="${LIGAND_NAME}_${RECEPTOR_NAME}"
INPUT_COMPLEX="${OUTPUT_DIR}/${COMPLEX_NAME}_input.pdb"

echo "=========================================="
echo "Step 1: Creating input complex"
echo "=========================================="

# Copy receptor atoms (keep original chains)
grep "^ATOM\|^TER\|^HETATM" "${RECEPTOR}" | grep -v "^END" > "${INPUT_COMPLEX}.tmp"

# Add TER after receptor
echo "TER" >> "${INPUT_COMPLEX}.tmp"

# Add peptide with chain P, renumbered residues starting at 1
awk 'BEGIN {
    res_num = 0
    last_res = ""
    atom_num = 10000
}
/^ATOM|^HETATM/ {
    atom_num++
    curr_res = substr($0, 23, 4)
    gsub(/[ ]+/, "", curr_res)
    
    if (curr_res != last_res) {
        res_num++
        last_res = curr_res
    }
    
    record = substr($0, 1, 6)
    atom_name = substr($0, 13, 4)
    alt_loc = substr($0, 17, 1)
    res_name = substr($0, 18, 3)
    rest = substr($0, 27)
    
    if (record ~ /HETATM/) {
        record = "ATOM  "
    }
    
    printf "%s%5d %s%s%s P%4d%s\n", record, atom_num, atom_name, alt_loc, res_name, res_num, rest
}' "${LIGAND}" >> "${INPUT_COMPLEX}.tmp"

echo "TER" >> "${INPUT_COMPLEX}.tmp"
echo "END" >> "${INPUT_COMPLEX}.tmp"
mv "${INPUT_COMPLEX}.tmp" "${INPUT_COMPLEX}"

echo "Input complex created: ${INPUT_COMPLEX}"

# Verify input complex
RECEPTOR_ATOMS=$(grep "^ATOM" "${INPUT_COMPLEX}" | grep -v " P " | wc -l)
PEPTIDE_ATOMS=$(grep "^ATOM" "${INPUT_COMPLEX}" | grep " P " | wc -l)
echo "  Receptor atoms: ${RECEPTOR_ATOMS}"
echo "  Peptide atoms: ${PEPTIDE_ATOMS}"

if [ "${PEPTIDE_ATOMS}" -eq 0 ]; then
    echo "ERROR: No peptide atoms in input complex!"
    exit 1
fi

echo ""
echo "=========================================="
echo "Step 2: Global Docking (docking_protocol)"
echo "=========================================="
echo ""
echo "This step performs TRUE global docking:"
echo "  - Randomizes peptide position around receptor"
echo "  - Searches entire receptor surface"
echo "  - Generates ${NSTRUCT} candidate poses"
echo ""

# Define partner chains for docking
# Format: receptor_chains_peptide_chain (e.g., "ABCD_P")
PARTNERS="${RECEPTOR_CHAINS}_P"
echo "Docking partners: ${PARTNERS}"
echo ""

# Run global docking using docking_protocol
# Key flags:
#   -partners: Defines which chains dock to which
#   -dock_pert: Translational and rotational perturbation
#   -spin: Random spin of peptide around docking axis
#   -randomize1/2: Randomize starting positions
#   -docking_local_refine: Local refinement after global search

run_rosetta docking_protocol.default.linuxgccrelease \
    -s "${INPUT_COMPLEX}" \
    -partners "${PARTNERS}" \
    -dock_pert 10 20 \
    -spin \
    -randomize1 \
    -randomize2 \
    -nstruct ${NSTRUCT} \
    -out:path:pdb "${OUTPUT_DIR}" \
    -out:file:scorefile "${OUTPUT_DIR}/${COMPLEX_NAME}_global_scores.sc" \
    -out:prefix "${COMPLEX_NAME}_global_" \
    -out:suffix "" \
    -docking_local_refine \
    -ex1 \
    -ex2aro \
    -score:weights ref2015 \
    -overwrite 2>&1 | tee "${OUTPUT_DIR}/${COMPLEX_NAME}_global_dock.log"

EXIT_CODE=$?

echo ""
echo "=========================================="
echo "Step 3: Analyzing global docking results"
echo "=========================================="

if [ ${EXIT_CODE} -ne 0 ]; then
    echo "ERROR: Global docking failed for ${LIGAND_NAME}"
    exit 1
fi

# Find output files
OUTPUT_FILES=($(ls "${OUTPUT_DIR}/${COMPLEX_NAME}_global_"*.pdb 2>/dev/null | grep -v "_input.pdb"))
echo "Generated ${#OUTPUT_FILES[@]} structures"

if [ ${#OUTPUT_FILES[@]} -eq 0 ]; then
    echo "WARNING: No output files generated!"
    exit 1
fi

# Analyze scores and find best structure
if [ -f "${OUTPUT_DIR}/${COMPLEX_NAME}_global_scores.sc" ]; then
    echo ""
    echo "Score analysis:"
    
    # Extract key metrics (I_sc = interface score is most relevant for docking)
    awk 'NR > 2 && NF > 1 {
        total = $2
        I_sc = $6  # Interface score column (may vary)
        count++
        sum_total += total
        if (NR == 3 || total < min_total) { min_total = total; best_struct = $NF }
        if (NR == 3 || total > max_total) max_total = total
    }
    END {
        if (count > 0) {
            print "  Structures analyzed:", count
            print "  total_score: min=" min_total ", max=" max_total ", avg=" sum_total/count
            print "  Best structure:", best_struct
        }
    }' "${OUTPUT_DIR}/${COMPLEX_NAME}_global_scores.sc"
    
    # Get the best structure (lowest total_score)
    BEST_LINE=$(awk 'NR > 2 && NF > 1 {print $0}' "${OUTPUT_DIR}/${COMPLEX_NAME}_global_scores.sc" | \
               sort -k2 -n | head -1)
    
    if [ -n "${BEST_LINE}" ]; then
        BEST_SCORE=$(echo "${BEST_LINE}" | awk '{print $2}')
        BEST_DESC=$(echo "${BEST_LINE}" | awk '{print $NF}')
        
        echo ""
        echo "Best global dock pose: ${BEST_DESC}"
        echo "  total_score: ${BEST_SCORE}"
        
        # Copy best structure for FlexPepDocking refinement
        BEST_FILE=$(ls "${OUTPUT_DIR}/"*"${BEST_DESC}"*.pdb 2>/dev/null | head -1)
        
        if [ -f "${BEST_FILE}" ]; then
            BEST_OUTPUT="${OUTPUT_DIR}/${COMPLEX_NAME}_best_global.pdb"
            cp "${BEST_FILE}" "${BEST_OUTPUT}"
            echo "✓ Best structure saved to: ${BEST_OUTPUT}"
        fi
    fi
    
    # Also save top 5 for diversity
    echo ""
    echo "Saving top 5 poses for potential refinement..."
    TOP_LINES=$(awk 'NR > 2 && NF > 1 {print $0}' "${OUTPUT_DIR}/${COMPLEX_NAME}_global_scores.sc" | \
               sort -k2 -n | head -5)
    
    TOP_NUM=1
    while IFS= read -r line; do
        DESC=$(echo "${line}" | awk '{print $NF}')
        TOP_FILE=$(ls "${OUTPUT_DIR}/"*"${DESC}"*.pdb 2>/dev/null | head -1)
        if [ -f "${TOP_FILE}" ]; then
            TOP_OUTPUT="${OUTPUT_DIR}/${COMPLEX_NAME}_top${TOP_NUM}_global.pdb"
            cp "${TOP_FILE}" "${TOP_OUTPUT}"
            echo "  ✓ Top ${TOP_NUM}: ${TOP_OUTPUT}"
            ((TOP_NUM++))
        fi
    done <<< "${TOP_LINES}"
fi

echo ""
echo "=========================================="
echo "Global docking complete for ${LIGAND_NAME}!"
echo "=========================================="
echo "Completed: $(date)"
echo ""
echo "Output files:"
echo "  Input complex: ${INPUT_COMPLEX}"
echo "  Score file: ${OUTPUT_DIR}/${COMPLEX_NAME}_global_scores.sc"
echo "  Best structure: ${OUTPUT_DIR}/${COMPLEX_NAME}_best_global.pdb"
echo "  Log file: ${OUTPUT_DIR}/${COMPLEX_NAME}_global_dock.log"
echo ""
echo "Next step: Run FlexPepDocking on the best global poses for refinement"
echo ""

