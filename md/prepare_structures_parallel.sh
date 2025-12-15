#!/bin/bash
#SBATCH --job-name=prep_md
#SBATCH --output=logs/prep_md_%A_%a.out
#SBATCH --error=logs/prep_md_%A_%a.err
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --mem=32G
#SBATCH --time=3:00:00
# NOTE: Array size is set dynamically by run_full_pipeline.sh via --array=0-N
#
# GROMACS Parallel Structure Preparation (V2)
# ===========================================
# V2 focus: preserve peptide vs receptor identity as distinct molecule types and
# emit an explicit analysis index (analysis.ndx) that uses moltype selections.
#
# Key best-practice changes vs v1:
# - pdb2gmx uses: -chainsep id_or_ter -merge no
#   This keeps chains as separate molecule types (e.g., Protein_chain_P).
# - After md.tpr is generated, create analysis.ndx:
#     [ Receptor ], [ Peptide ], [ PeptideBackbone ]
#

set -euo pipefail

# Load required modules
module load apptainer/1.4.1-gcc-13.4.0

# Configuration
MD_DIR="/home/kuhfeldr/hpc.cqls/md"
CONTAINER_PATH="/home/kuhfeldr/gromacs_v2025.3.sif"

# Create logs directory
mkdir -p "${MD_DIR}/logs"

# Get structures from environment variable (set by run_full_pipeline.sh) or fallback to file
STRUCTURES=()
if [ -n "${STRUCTURES_LIST}" ]; then
    # Parse comma-separated list from environment variable
    # Use readarray for more reliable parsing
    IFS=',' read -ra STRUCTURES <<< "${STRUCTURES_LIST}"
    # Filter out any empty elements
    TEMP_STRUCTURES=()
    for struct in "${STRUCTURES[@]}"; do
        [ -n "$struct" ] && TEMP_STRUCTURES+=("$struct")
    done
    STRUCTURES=("${TEMP_STRUCTURES[@]}")
else
    # Fallback: read from file (for backward compatibility)
    STRUCTURES_FILE="${MD_DIR}/structures_to_process.txt"
    if [ ! -f "${STRUCTURES_FILE}" ]; then
        echo "ERROR: No structures provided!"
        echo "This script should be run via run_full_pipeline.sh which sets STRUCTURES_LIST"
        echo "  OR provide structures_to_process.txt file"
        exit 1
    fi
    
    # Load structures into array (filter out empty lines)
    while IFS= read -r line; do
        [ -z "$line" ] && continue
        # Remove .pdb extension if present
        struct_name="${line%.pdb}"
        STRUCTURES+=("$struct_name")
    done < "${STRUCTURES_FILE}"
fi

# Get structure for this array task (default to 0 if not run as array job)
SLURM_ARRAY_TASK_ID="${SLURM_ARRAY_TASK_ID:-0}"

# Verify we have structures
if [ ${#STRUCTURES[@]} -eq 0 ]; then
    echo "ERROR: No structures found in STRUCTURES_LIST or file!"
    echo "STRUCTURES_LIST=${STRUCTURES_LIST}"
    exit 1
fi

if [ "${SLURM_ARRAY_TASK_ID}" -ge "${#STRUCTURES[@]}" ]; then
    echo "No work for task ${SLURM_ARRAY_TASK_ID} (array has ${#STRUCTURES[@]} elements: ${STRUCTURES[*]})"
    exit 0
fi

STRUCTURE_NAME="${STRUCTURES[${SLURM_ARRAY_TASK_ID}]}"
# Try with _3fxi suffix first, then without
if [ -f "${MD_DIR}/${STRUCTURE_NAME}_3fxi.pdb" ]; then
    PDB_FILE="${MD_DIR}/${STRUCTURE_NAME}_3fxi.pdb"
    STRUCTURE_NAME="${STRUCTURE_NAME}_3fxi"  # Update to match directory name
    PREP_DIR="${MD_DIR}/${STRUCTURE_NAME}"
elif [ -f "${MD_DIR}/${STRUCTURE_NAME}.pdb" ]; then
    PDB_FILE="${MD_DIR}/${STRUCTURE_NAME}.pdb"
    PREP_DIR="${MD_DIR}/${STRUCTURE_NAME}"
else
    # Try the other way - if structure name already has _3fxi, use it as-is
    if [[ "${STRUCTURE_NAME}" == *_3fxi ]]; then
        PDB_FILE="${MD_DIR}/${STRUCTURE_NAME}.pdb"
        PREP_DIR="${MD_DIR}/${STRUCTURE_NAME}"
    else
        echo "ERROR: PDB file not found for structure: ${STRUCTURES[${SLURM_ARRAY_TASK_ID}]}"
        echo "  Tried: ${MD_DIR}/${STRUCTURE_NAME}_3fxi.pdb"
        echo "  Tried: ${MD_DIR}/${STRUCTURE_NAME}.pdb"
        exit 1
    fi
fi

echo "=========================================="
echo "GROMACS Structure Preparation (Parallel) - V2"
echo "=========================================="
echo "Structure: ${STRUCTURE_NAME}"
echo "PDB File: ${PDB_FILE}"
echo "Prep Dir: ${PREP_DIR}"
echo "Array Task: ${SLURM_ARRAY_TASK_ID}"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Started: $(date)"
echo "Node: $(hostname)"
echo ""

# Check if already prepared
if [ -f "${PREP_DIR}/md.tpr" ]; then
    echo "✓ Structure already prepared (md.tpr exists)"
    echo "  Skipping: ${STRUCTURE_NAME}"
    exit 0
fi

# PDB file existence was already checked above

# Verify container exists
if [ ! -f "${CONTAINER_PATH}" ]; then
    echo "ERROR: GROMACS container not found at: ${CONTAINER_PATH}"
    exit 1
fi

# Create and enter preparation directory
mkdir -p "${PREP_DIR}"
cd "${PREP_DIR}"

# Copy input files
cp "${MD_DIR}"/*.mdp . 2>/dev/null || { echo "ERROR: No .mdp files found"; exit 1; }
cp "${PDB_FILE}" .

# Set up GROMACS container
CPUS=${SLURM_CPUS_PER_TASK:-8}
SINGULARITY="apptainer run --nv -B ${PWD}:/host_pwd --pwd /host_pwd ${CONTAINER_PATH}"

run_gmx() {
    ${SINGULARITY} gmx "$@"
}

echo "GPU: $(nvidia-smi -L 2>/dev/null | head -1 || echo 'N/A')"
echo "CPUs: ${CPUS}"
echo ""

PDB_BASENAME="$(basename "${PDB_FILE}")"

# Step 1: Generate topology
echo "Step 1: Generating topology (pdb2gmx)..."
echo "  V2: Using -chainsep id_or_ter -merge no to preserve chain/molecule identity"
echo -e "1\n1" | run_gmx pdb2gmx \
    -f "${PDB_BASENAME}" \
    -o processed.gro \
    -water tip3p \
    -ignh \
    -chainsep id_or_ter \
    -merge no 2>&1 | tee pdb2gmx.log
[ ${PIPESTATUS[1]} -ne 0 ] && { echo "ERROR: pdb2gmx failed"; exit 1; }
echo "✓ Topology generated"

# Step 2: Create box
echo "Step 2: Creating simulation box..."
run_gmx editconf -f processed.gro -o boxed.gro -c -d 1.2 -bt cubic
echo "✓ Box created"

# Step 3: Solvate
echo "Step 3: Solvating system..."
run_gmx solvate -cp boxed.gro -cs spc216.gro -o solvated.gro -p topol.top
echo "✓ System solvated"

# Step 4: Add ions
echo "Step 4: Adding ions..."
run_gmx grompp -f ions.mdp -c solvated.gro -p topol.top -o ions.tpr -maxwarn 1
echo "SOL" | run_gmx genion -s ions.tpr -o ionized.gro -p topol.top -pname NA -nname CL -neutral -conc 0.15
echo "✓ Ions added"

# Step 5: Energy minimization
echo "Step 5: Energy minimization..."
run_gmx grompp -f minim.mdp -c ionized.gro -p topol.top -o em.tpr -maxwarn 1
run_gmx mdrun -v -deffnm em -ntmpi 1 -ntomp "${CPUS}" -nb gpu
echo "✓ Minimization complete"

# Step 6: NVT equilibration
echo "Step 6: NVT equilibration (250 ps)..."
run_gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr -maxwarn 1
run_gmx mdrun -deffnm nvt -ntmpi 1 -ntomp "${CPUS}" -nb gpu -pme gpu -update gpu -v
echo "✓ NVT complete (250 ps)"

# Step 7: NPT equilibration
echo "Step 7: NPT equilibration (500 ps)..."
run_gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr -maxwarn 1
run_gmx mdrun -deffnm npt -ntmpi 1 -ntomp "${CPUS}" -nb gpu -pme gpu -update gpu -v
echo "✓ NPT complete (500 ps)"

# Step 8: Generate production TPR
echo "Step 8: Generating production MD .tpr..."
run_gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr -maxwarn 1
echo "✓ md.tpr generated"

# Step 9 (V2): Create explicit analysis.ndx using moltype selections
echo "Step 9: Creating analysis.ndx (moltype-based)..."
NDX_OK=false
if [ -f "md.tpr" ] && [ -f "topol_Protein_chain_P.itp" ]; then
    receptor_sel=""
    for itp in topol_Protein_chain_*.itp; do
        [ -f "${itp}" ] || continue
        chain_letter="$(echo "${itp}" | sed -E 's/.*chain_(.)\\.itp/\\1/')"
        moltype="Protein_chain_${chain_letter}"
        if [ "${chain_letter}" = "P" ]; then
            continue
        fi
        if [ -z "${receptor_sel}" ]; then
            receptor_sel="moltype \\\"${moltype}\\\""
        else
            receptor_sel="${receptor_sel} or moltype \\\"${moltype}\\\""
        fi
    done

    if [ -n "${receptor_sel}" ]; then
        cat > analysis_selections.txt << EOF
Receptor = ${receptor_sel}
Peptide = moltype "Protein_chain_P"
PeptideBackbone = moltype "Protein_chain_P" and backbone
EOF
        if run_gmx select -s md.tpr -sf analysis_selections.txt -on analysis.ndx >/dev/null 2>&1; then
            echo "✓ analysis.ndx created (Receptor/Peptide/PeptideBackbone)"
            NDX_OK=true
        fi
    fi
fi

if [ "${NDX_OK}" != true ]; then
    echo "WARNING: Failed to create analysis.ndx via moltype selections"
    echo "         (Downstream analyze_md_v2.sh will fall back to residue-range selection.)"
fi

echo ""
echo "=========================================="
echo "PREPARATION COMPLETE (V2): ${STRUCTURE_NAME}"
echo "=========================================="
echo "Completed: $(date)"
echo "Output: ${PREP_DIR}/md.tpr"
echo ""


