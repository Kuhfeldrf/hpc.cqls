#!/bin/bash
# Fix existing receptor-peptide PDB files using GROMACS tools
# This uses gmx editconf and gmx pdb2gmx for reliable format conversion

MD_DIR="/home/pi/kuhfeldr/md"
FILES_TO_FIX=(
    "GLAPYKLRPVAA_recligand_lr1.pdb"
    "LLFKDSAIGF_recligand_lr6.pdb"
    "RPKLPLRYP_recligand_lr2.pdb"
)

echo "=========================================="
echo "Fixing existing PDB files with GROMACS tools"
echo "=========================================="
echo ""

cd "${MD_DIR}"

for PDB_FILE in "${FILES_TO_FIX[@]}"; do
    if [ ! -f "${PDB_FILE}" ]; then
        echo "SKIP: ${PDB_FILE} not found"
        continue
    fi
    
    echo "Processing: ${PDB_FILE}"
    
    # Step 1: Convert HETATM to ATOM for standard amino acids (safe operation)
    TMP_FILE=$(mktemp)
    standard_aa="ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL"
    
    awk -v aa_list="$standard_aa" '
    /^HETATM/ {
        res_name = substr($0, 18, 3)
        gsub(/^[ ]+|[ ]+$/, "", res_name)
        
        # Only convert if it's a standard amino acid
        if (index(aa_list, res_name) > 0) {
            # Replace HETATM with ATOM (maintain exact spacing)
            $0 = "ATOM  " substr($0, 7)
        }
    }
    { print }
    ' "${PDB_FILE}" > "${TMP_FILE}"
    
    # Step 2: Use gmx editconf to validate and convert format
    # This uses GROMACS's robust parser
    TMP_GRO="${TMP_FILE}.gro"
    
    echo "  Validating with gmx editconf..."
    gmx editconf -f "${TMP_FILE}" -o "${TMP_GRO}" -quiet > /dev/null 2>&1
    
    if [ $? -ne 0 ]; then
        echo "  ERROR: gmx editconf failed. Checking with pdb2gmx for specific errors..."
        
        # Try pdb2gmx to see what the actual problem is
        gmx pdb2gmx -f "${TMP_FILE}" -o /dev/null -water none -ignh > "${PDB_FILE}.error.log" 2>&1
        
        if grep -q "atom CA.*not found" "${PDB_FILE}.error.log"; then
            echo "  ERROR: Missing CA atoms detected in ${PDB_FILE}"
            echo "         This file needs manual fixing or regeneration."
            echo "         Error details: ${PDB_FILE}.error.log"
            rm -f "${TMP_FILE}" "${TMP_GRO}"
            continue
        else
            echo "  WARNING: Unknown error. Check ${PDB_FILE}.error.log"
            rm -f "${TMP_FILE}" "${TMP_GRO}"
            continue
        fi
    fi
    
    # Step 3: Convert back to PDB format using editconf
    # This ensures proper formatting
    BACKUP_FILE="${PDB_FILE}.backup"
    FIXED_FILE="${PDB_FILE}"
    
    echo "  Converting back to PDB format..."
    gmx editconf -f "${TMP_GRO}" -o "${FIXED_FILE}" -quiet > /dev/null 2>&1
    
    if [ $? -eq 0 ]; then
        # Backup original
        cp "${PDB_FILE}" "${BACKUP_FILE}"
        
        # Replace with fixed version
        mv "${FIXED_FILE}" "${PDB_FILE}"
        
        echo "  ✓ Fixed: ${PDB_FILE}"
        echo "    Backup saved as: ${BACKUP_FILE}"
        
        # Verify the fix worked
        if grep -q "^ATOM.*CA" "${PDB_FILE}"; then
            echo "    ✓ Contains CA atoms"
        else
            echo "    WARNING: Still no CA atoms found"
        fi
    else
        echo "  ERROR: Failed to convert back to PDB format"
    fi
    
    # Clean up
    rm -f "${TMP_FILE}" "${TMP_GRO}"
    echo ""
done

echo "=========================================="
echo "Processing complete!"
echo ""
echo "Next step: Test with pdb2gmx:"
echo "  cd ${MD_DIR}"
echo "  gmx pdb2gmx -f <fixed_file.pdb> -o test.gro -water none"
echo "=========================================="

