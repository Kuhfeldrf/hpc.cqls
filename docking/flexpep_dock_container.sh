#!/bin/bash
#SBATCH --job-name=flexpep
#SBATCH --output=/home/kuhfeldr/hpc.cqls/docking/logs/flexpep_%A_%a.out
#SBATCH --error=/home/kuhfeldr/hpc.cqls/docking/logs/flexpep_%A_%a.err
#SBATCH --partition=normal
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --array=0-4

# FlexPepDocking using Rosetta Container
# =======================================
# Uses apptainer/singularity with rosetta.sif container
# Processes 5 peptides from alphafold_predictions/

set -e

# Load apptainer module
module load apptainer/1.4.1-gcc-13.4.0

# Paths
WORK_DIR="/home/kuhfeldr/hpc.cqls/docking"
RECEPTOR_DIR="${WORK_DIR}/receptor"
LIGAND_DIR="${WORK_DIR}/peptides"
OUTPUT_DIR="${WORK_DIR}/docked_complexes"
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

# Number of structures to generate
NSTRUCT=10

# Create output and log directories
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${WORK_DIR}/logs"

cd "${WORK_DIR}"

# Get list of ligands (sorted for reproducibility)
LIGANDS=($(ls ${LIGAND_DIR}/*.pdb 2>/dev/null | sort))

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
echo "FlexPepDocking (Container Version)"
echo "=========================================="
echo "Started: $(date)"
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Container: ${ROSETTA_SIF}"
echo "Processing: ${LIGAND_NAME}"
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
echo "Step 2: Pre-relaxation (removes clashes)"
echo "=========================================="

echo "Running Rosetta relax to minimize initial clashes..."

run_rosetta relax.default.linuxgccrelease \
    -s "${INPUT_COMPLEX}" \
    -relax:quick \
    -relax:constrain_relax_to_start_coords \
    -relax:coord_constrain_sidechains \
    -relax:ramp_constraints false \
    -out:path:pdb "${OUTPUT_DIR}" \
    -out:prefix "prerelax_" \
    -out:suffix "" \
    -nstruct 1 \
    -score:weights ref2015 \
    -overwrite 2>&1 | tee "${OUTPUT_DIR}/${COMPLEX_NAME}_relax.log"

# Find relaxed structure
RELAXED_FILE=$(ls "${OUTPUT_DIR}/prerelax_${COMPLEX_NAME}_input_"*.pdb 2>/dev/null | head -1)

if [ -f "${RELAXED_FILE}" ]; then
    mv "${RELAXED_FILE}" "${RELAXED_COMPLEX}"
    echo "✓ Pre-relaxation complete: ${RELAXED_COMPLEX}"
    INPUT_FOR_DOCKING="${RELAXED_COMPLEX}"
else
    echo "⚠ Pre-relaxation failed, using unrelaxed input"
    INPUT_FOR_DOCKING="${INPUT_COMPLEX}"
fi

echo ""
echo "=========================================="
echo "Step 3: FlexPepDocking with full refinement"
echo "=========================================="
echo "  Input: ${INPUT_FOR_DOCKING}"
echo "  Receptor chains: ${RECEPTOR_CHAINS}"
echo "  Peptide chain: P"
echo "  Structures: ${NSTRUCT}"
echo ""

# Run FlexPepDocking with PROPER refinement protocol
run_rosetta FlexPepDocking.default.linuxgccrelease \
    -s "${INPUT_FOR_DOCKING}" \
    -flexPepDocking:receptor_chain "${RECEPTOR_CHAINS}" \
    -flexPepDocking:peptide_chain "P" \
    -nstruct ${NSTRUCT} \
    -out:path:pdb "${OUTPUT_DIR}" \
    -out:file:scorefile "${OUTPUT_DIR}/${COMPLEX_NAME}_scores.sc" \
    -out:prefix "${COMPLEX_NAME}_" \
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
    -overwrite 2>&1 | tee "${OUTPUT_DIR}/${COMPLEX_NAME}_docking.log"

EXIT_CODE=$?

if [ ${EXIT_CODE} -ne 0 ]; then
    echo ""
    echo "⚠ Primary docking failed, trying alternative protocol..."
    
    # Alternative: refine_only mode
    run_rosetta FlexPepDocking.default.linuxgccrelease \
        -s "${INPUT_FOR_DOCKING}" \
        -flexPepDocking:receptor_chain "${RECEPTOR_CHAINS}" \
        -flexPepDocking:peptide_chain "P" \
        -nstruct $((NSTRUCT / 2)) \
        -out:path:pdb "${OUTPUT_DIR}" \
        -out:file:scorefile "${OUTPUT_DIR}/${COMPLEX_NAME}_scores.sc" \
        -out:prefix "${COMPLEX_NAME}_" \
        -out:suffix "_docked" \
        -flexPepDocking:refine_only \
        -ex1 \
        -ex2aro \
        -use_input_sc \
        -score:weights ref2015 \
        -renumber_pdb false \
        -overwrite 2>&1 | tee -a "${OUTPUT_DIR}/${COMPLEX_NAME}_docking.log"
    
    if [ $? -ne 0 ]; then
        echo "ERROR: FlexPepDocking failed for ${LIGAND_NAME}"
        exit 1
    fi
fi

echo ""
echo "=========================================="
echo "Step 4: Analyzing results"
echo "=========================================="

# Find output files
OUTPUT_FILES=($(ls "${OUTPUT_DIR}/${COMPLEX_NAME}_"*_docked*.pdb 2>/dev/null))
echo "Generated ${#OUTPUT_FILES[@]} structures"

if [ ${#OUTPUT_FILES[@]} -eq 0 ]; then
    echo "WARNING: No output files generated!"
else
    # Analyze scores
    if [ -f "${OUTPUT_DIR}/${COMPLEX_NAME}_scores.sc" ]; then
        echo ""
        echo "Score analysis:"
        
        # Extract key metrics
        awk 'NR > 2 && NF > 1 {
            total = $2
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
        }' "${OUTPUT_DIR}/${COMPLEX_NAME}_scores.sc"
        
        # Get best structure and copy to md folder
        BEST_LINE=$(awk 'NR > 2 && NF > 1 {print $0}' "${OUTPUT_DIR}/${COMPLEX_NAME}_scores.sc" | \
                   sort -k2 -n | head -1)
        
        if [ -n "${BEST_LINE}" ]; then
            BEST_SCORE=$(echo "${BEST_LINE}" | awk '{print $2}')
            BEST_DESC=$(echo "${BEST_LINE}" | awk '{print $NF}')
            
            echo ""
            echo "Best structure: ${BEST_DESC}"
            echo "  total_score: ${BEST_SCORE}"
            
            # Find and copy the best file
            BEST_FILE=$(ls "${OUTPUT_DIR}/"*"${BEST_DESC}"*.pdb 2>/dev/null | head -1)
            
            if [ -f "${BEST_FILE}" ]; then
                mkdir -p "/home/kuhfeldr/hpc.cqls/md"
                FINAL_NAME="${LIGAND_NAME}_${RECEPTOR_NAME}_docked_best"
                cp "${BEST_FILE}" "/home/kuhfeldr/hpc.cqls/md/${FINAL_NAME}.pdb"
                echo "✓ Best structure copied to: /home/kuhfeldr/hpc.cqls/md/${FINAL_NAME}.pdb"
            fi
        fi
    fi
fi

echo ""
echo "=========================================="
echo "FlexPepDocking complete for ${LIGAND_NAME}!"
echo "=========================================="
echo "Completed: $(date)"
echo ""
echo "Output files: ${OUTPUT_DIR}/"
echo "Score file: ${OUTPUT_DIR}/${COMPLEX_NAME}_scores.sc"
echo ""

