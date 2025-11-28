#!/bin/bash
#SBATCH --job-name=flexpep_dock_3fxi
#SBATCH --output=/home/pi/kuhfeldr/docking/flexpep_dock_%j.out
#SBATCH --error=/home/pi/kuhfeldr/docking/flexpep_dock_%j.err
#SBATCH --partition=cqls_gpu-1080
#SBATCH --account=cqls
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

# FlexPepDocking script for docking peptides to 3fxi_0001.pdb
# Uses FlexPepDocking (the proper Rosetta tool for peptide-protein docking)
# Outputs are GROMACS-compatible: standard ATOM records and atom naming

set -e

cd /home/pi/kuhfeldr/docking

RECEPTOR="3fxi_0001.pdb"
LIGAND_DIR="alphafold_predictions"
OUTPUT_DIR="docked_complexes"
ROSETTA_BIN="/home/pi/kuhfeldr/rosetta/source/bin/FlexPepDocking.linuxgccrelease"

# Create output directory
mkdir -p "${OUTPUT_DIR}"

echo "=========================================="
echo "FlexPepDocking: Peptides to ${RECEPTOR}"
echo "=========================================="
echo "Started: $(date)"
echo ""

# Check Rosetta binary exists
if [ ! -f "${ROSETTA_BIN}" ]; then
    echo "ERROR: Rosetta FlexPepDocking not found at ${ROSETTA_BIN}"
    exit 1
fi

# Get receptor chain IDs
RECEPTOR_CHAINS=$(grep "^ATOM" "${RECEPTOR}" | awk '{print substr($0,22,1)}' | sort -u | tr -d '\n')
echo "Receptor chains found: ${RECEPTOR_CHAINS}"

# Get list of ligands
LIGANDS=($(ls ${LIGAND_DIR}/*.pdb 2>/dev/null))

if [ ${#LIGANDS[@]} -eq 0 ]; then
    echo "ERROR: No ligand PDB files found in ${LIGAND_DIR}/"
    exit 1
fi

echo "Found ${#LIGANDS[@]} peptide(s) to dock:"
for LIG in "${LIGANDS[@]}"; do
    echo "  - $(basename $LIG)"
done
echo ""

# Process each peptide ligand
for LIGAND in "${LIGANDS[@]}"; do
    LIGAND_NAME=$(basename "${LIGAND}" .pdb)
    COMPLEX_NAME="${LIGAND_NAME}_3fxi"
    INPUT_COMPLEX="${OUTPUT_DIR}/${COMPLEX_NAME}_input.pdb"
    
    echo "=========================================="
    echo "Processing: ${LIGAND_NAME}"
    echo "=========================================="
    
    # Step 1: Create input complex with proper chain assignments
    # Receptor: chains A, B, C, D (as in original)
    # Peptide: chain P (use P for peptide to be explicit)
    echo "Creating input complex with proper chain assignments..."
    
    # Copy receptor atoms (keep original chains)
    grep "^ATOM\|^TER\|^HETATM" "${RECEPTOR}" | grep -v "^END" > "${INPUT_COMPLEX}.tmp"
    
    # Add TER after receptor
    echo "TER" >> "${INPUT_COMPLEX}.tmp"
    
    # Add peptide with chain P, renumbered residues starting at 1
    # Use awk to properly format PDB columns
    awk 'BEGIN {
        res_num = 0
        last_res = ""
        atom_num = 10000  # Start after receptor atoms
    }
    /^ATOM|^HETATM/ {
        atom_num++
        # Get current residue number from input
        curr_res = substr($0, 23, 4)
        gsub(/[ ]+/, "", curr_res)
        
        if (curr_res != last_res) {
            res_num++
            last_res = curr_res
        }
        
        # Format new line with chain P and new residue number
        # PDB format: columns 1-6(record), 7-11(atom#), 13-16(atom name), 17(alt), 18-20(res name), 22(chain), 23-26(res#), rest unchanged
        record = substr($0, 1, 6)
        atom_name = substr($0, 13, 4)
        alt_loc = substr($0, 17, 1)
        res_name = substr($0, 18, 3)
        rest = substr($0, 27)
        
        # Force ATOM record for standard amino acids
        if (record ~ /HETATM/) {
            record = "ATOM  "
        }
        
        printf "%s%5d %s%s%s P%4d%s\n", record, atom_num, atom_name, alt_loc, res_name, res_num, rest
    }' "${LIGAND}" >> "${INPUT_COMPLEX}.tmp"
    
    # Add final TER and END
    echo "TER" >> "${INPUT_COMPLEX}.tmp"
    echo "END" >> "${INPUT_COMPLEX}.tmp"
    
    mv "${INPUT_COMPLEX}.tmp" "${INPUT_COMPLEX}"
    
    echo "Input complex created: ${INPUT_COMPLEX}"
    
    # Verify input complex
    echo "Verifying input complex..."
    RECEPTOR_ATOMS=$(grep "^ATOM" "${INPUT_COMPLEX}" | grep -v " P " | wc -l)
    PEPTIDE_ATOMS=$(grep "^ATOM" "${INPUT_COMPLEX}" | grep " P " | wc -l)
    echo "  Receptor atoms: ${RECEPTOR_ATOMS}"
    echo "  Peptide atoms: ${PEPTIDE_ATOMS}"
    
    if [ "${PEPTIDE_ATOMS}" -eq 0 ]; then
        echo "ERROR: No peptide atoms in input complex!"
        continue
    fi
    
    echo ""
    echo "Running FlexPepDocking..."
    echo "  Receptor chains: ${RECEPTOR_CHAINS}"
    echo "  Peptide chain: P"
    
    # Run FlexPepDocking with GROMACS-compatible output settings
    "${ROSETTA_BIN}" \
        -s "${INPUT_COMPLEX}" \
        -flexPepDocking:receptor_chain "${RECEPTOR_CHAINS}" \
        -flexPepDocking:peptide_chain "P" \
        -nstruct 10 \
        -out:path:pdb "${OUTPUT_DIR}" \
        -out:file:scorefile "${OUTPUT_DIR}/${COMPLEX_NAME}_scores.sc" \
        -out:prefix "${COMPLEX_NAME}_" \
        -out:suffix "_docked" \
        -pep_refine \
        -ex1 -ex2 \
        -use_input_sc \
        -unboundrot "${LIGAND}" \
        -score:weights ref2015 \
        -renumber_pdb false \
        -preserve_header false \
        -packing:repack_only \
        -overwrite
    
    EXIT_CODE=$?
    
    if [ ${EXIT_CODE} -ne 0 ]; then
        echo "ERROR: FlexPepDocking failed for ${LIGAND_NAME} (exit code: ${EXIT_CODE})"
        echo "Trying with lowres_preoptimize for difficult cases..."
        
        # Try with additional optimization for difficult cases
        "${ROSETTA_BIN}" \
            -s "${INPUT_COMPLEX}" \
            -flexPepDocking:receptor_chain "${RECEPTOR_CHAINS}" \
            -flexPepDocking:peptide_chain "P" \
            -nstruct 5 \
            -out:path:pdb "${OUTPUT_DIR}" \
            -out:file:scorefile "${OUTPUT_DIR}/${COMPLEX_NAME}_scores.sc" \
            -out:prefix "${COMPLEX_NAME}_" \
            -out:suffix "_docked" \
            -lowres_preoptimize \
            -pep_refine \
            -ex1 \
            -use_input_sc \
            -score:weights ref2015 \
            -renumber_pdb false \
            -overwrite
        
        if [ $? -ne 0 ]; then
            echo "ERROR: FlexPepDocking failed even with lowres_preoptimize for ${LIGAND_NAME}"
            continue
        fi
    fi
    
    echo "FlexPepDocking completed for ${LIGAND_NAME}"
    
    # Verify outputs are GROMACS-compatible
    echo ""
    echo "Verifying output PDB format..."
    
    # Find output files
    OUTPUT_FILES=($(ls "${OUTPUT_DIR}/${COMPLEX_NAME}_"*_docked*.pdb 2>/dev/null | head -5))
    
    if [ ${#OUTPUT_FILES[@]} -eq 0 ]; then
        echo "WARNING: No output files found for ${LIGAND_NAME}"
        continue
    fi
    
    # Check first output file
    SAMPLE_OUTPUT="${OUTPUT_FILES[0]}"
    echo "Checking: $(basename ${SAMPLE_OUTPUT})"
    
    # Check for HETATM (should be minimal or none for standard amino acids)
    HETATM_COUNT=$(grep "^HETATM" "${SAMPLE_OUTPUT}" 2>/dev/null | wc -l)
    ATOM_COUNT=$(grep "^ATOM" "${SAMPLE_OUTPUT}" 2>/dev/null | wc -l)
    CA_COUNT=$(grep "^ATOM" "${SAMPLE_OUTPUT}" 2>/dev/null | grep " CA " | wc -l)
    
    echo "  ATOM records: ${ATOM_COUNT}"
    echo "  HETATM records: ${HETATM_COUNT}"
    echo "  CA atoms: ${CA_COUNT}"
    
    # Quick GROMACS compatibility check - look for standard atom names
    if [ ${CA_COUNT} -gt 0 ]; then
        echo "  ✓ Contains standard CA atoms (GROMACS-compatible)"
    else
        echo "  ⚠ WARNING: No CA atoms found - may need manual review"
    fi
    
    # Find best structure from scorefile
    if [ -f "${OUTPUT_DIR}/${COMPLEX_NAME}_scores.sc" ]; then
        echo ""
        echo "Processing scores..."
        
        # Get best structure (lowest total_score)
        BEST_LINE=$(awk 'NR > 2 && NF > 1 {print $0}' "${OUTPUT_DIR}/${COMPLEX_NAME}_scores.sc" | \
                   sort -k2 -n | head -1)
        
        if [ -n "${BEST_LINE}" ]; then
            BEST_SCORE=$(echo "${BEST_LINE}" | awk '{print $2}')
            BEST_DESC=$(echo "${BEST_LINE}" | awk '{print $NF}')
            
            echo "Best structure: ${BEST_DESC}"
            echo "Score: ${BEST_SCORE}"
            
            # Find the actual file
            BEST_FILE=$(ls "${OUTPUT_DIR}/"*"${BEST_DESC}"*.pdb 2>/dev/null | head -1)
            
            if [ -f "${BEST_FILE}" ]; then
                # Copy to md/ directory with clean name
                FINAL_NAME="${LIGAND_NAME}_recligand"
                cp "${BEST_FILE}" "/home/pi/kuhfeldr/md/${FINAL_NAME}.pdb"
                echo "✓ Best structure copied to: md/${FINAL_NAME}.pdb"
            fi
        fi
    fi
    
    echo ""
done

echo "=========================================="
echo "FlexPepDocking complete!"
echo "=========================================="
echo "Completed: $(date)"
echo ""
echo "Output files in: ${OUTPUT_DIR}/"
echo "Best structures copied to: /home/pi/kuhfeldr/md/"
echo ""
echo "Next steps:"
echo "  1. Review docking scores in ${OUTPUT_DIR}/*_scores.sc"
echo "  2. Visually inspect docked structures in PyMOL/VMD"
echo "  3. Prepare for MD: cd /home/pi/kuhfeldr/md && ./prepare_all_structures.sh"
echo ""

