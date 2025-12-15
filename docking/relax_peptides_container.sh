#!/bin/bash
#SBATCH --job-name=relax_pep
#SBATCH --output=/home/kuhfeldr/hpc.cqls/docking/logs/relax_pep_%A_%a.out
#SBATCH --error=/home/kuhfeldr/hpc.cqls/docking/logs/relax_pep_%A_%a.err
#SBATCH --partition=normal
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G

# Peptide Relaxation using Rosetta Container
# ===========================================
# Uses apptainer/singularity with rosetta.sif container
# Only generates relaxed structures - no docking
# Processes peptides from peptides/ directory
#
# Environment variables (set by run_docking_pipeline.sh):
#   PEPTIDE_LIST_FILE - Path to peptide list file (optional)
#   PEPTIDE_DIR       - Directory containing peptide PDB files
#   OUTPUT_DIR        - Output directory for relaxed complexes

set -e

# Load apptainer module
module load apptainer/1.4.1-gcc-13.4.0

# Paths - use environment variables if set, otherwise use defaults
WORK_DIR="/home/kuhfeldr/hpc.cqls/docking"
RECEPTOR_DIR="${WORK_DIR}/receptor"
LIGAND_DIR="${PEPTIDE_DIR:-${WORK_DIR}/peptides}"
OUTPUT_DIR="${OUTPUT_DIR:-${WORK_DIR}/docked_complexes}"
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
echo "Peptide Relaxation (Container Version)"
echo "=========================================="
echo "Started: $(date)"
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Container: ${ROSETTA_SIF}"
echo "Processing: ${LIGAND_NAME}"
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

# Process the selected peptide ligand
COMPLEX_NAME="${LIGAND_NAME}_${RECEPTOR_NAME}"
INPUT_COMPLEX="${OUTPUT_DIR}/${COMPLEX_NAME}_input.pdb"
RELAXED_COMPLEX="${OUTPUT_DIR}/${COMPLEX_NAME}_relaxed.pdb"

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
echo "Step 2: Relaxation (removes clashes)"
echo "=========================================="

echo "Running Rosetta relax to minimize initial clashes..."

run_rosetta relax.default.linuxgccrelease \
    -s "${INPUT_COMPLEX}" \
    -relax:quick \
    -relax:constrain_relax_to_start_coords \
    -relax:coord_constrain_sidechains \
    -relax:ramp_constraints false \
    -out:path:pdb "${OUTPUT_DIR}" \
    -out:prefix "relax_" \
    -out:suffix "" \
    -nstruct 1 \
    -score:weights ref2015 \
    -overwrite 2>&1 | tee "${OUTPUT_DIR}/${COMPLEX_NAME}_relax.log"

EXIT_CODE=$?

# Find relaxed structure
RELAXED_FILE=$(ls "${OUTPUT_DIR}/relax_${COMPLEX_NAME}_input_"*.pdb 2>/dev/null | head -1)

if [ ${EXIT_CODE} -eq 0 ] && [ -f "${RELAXED_FILE}" ]; then
    mv "${RELAXED_FILE}" "${RELAXED_COMPLEX}"
    echo ""
    echo "✓ Relaxation complete: ${RELAXED_COMPLEX}"
    
    # Extract score from log
    SCORE=$(grep "protocols.relax.FastRelax: Final score:" "${OUTPUT_DIR}/${COMPLEX_NAME}_relax.log" 2>/dev/null | tail -1 | awk '{print $NF}')
    if [ -n "${SCORE}" ]; then
        echo "  Final score: ${SCORE}"
    fi
else
    echo ""
    echo "ERROR: Relaxation failed for ${LIGAND_NAME}"
    exit 1
fi

echo ""
echo "=========================================="
echo "Relaxation complete for ${LIGAND_NAME}!"
echo "=========================================="
echo "Completed: $(date)"
echo ""
echo "Output files:"
echo "  Input complex: ${INPUT_COMPLEX}"
echo "  Relaxed structure: ${RELAXED_COMPLEX}"
echo "  Log file: ${OUTPUT_DIR}/${COMPLEX_NAME}_relax.log"
echo ""

