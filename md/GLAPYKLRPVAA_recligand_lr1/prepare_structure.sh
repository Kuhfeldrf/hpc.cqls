#!/bin/bash
#SBATCH --job-name=prep_md
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --partition=cqls_gpu-1080
#SBATCH --account=cqls
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:1
#SBATCH --mem=32G
#SBATCH --time=12:00:00

# GROMACS Structure Preparation Script
# =====================================
# This script prepares a PDB structure for MD simulation:
# 1. Generate topology (pdb2gmx)
# 2. Create simulation box (editconf)
# 3. Solvate system (solvate)
# 4. Add ions (genion)
# 5. Energy minimization (grompp + mdrun)
# 6. NVT equilibration (grompp + mdrun)
# 7. NPT equilibration (grompp + mdrun)
# 8. Generate production MD .tpr file (grompp)
#
# Usage: sbatch prepare_structure.sh <structure_name.pdb>

# Check if filename argument is provided
if [ -z "$1" ]; then
    echo "Error: Please provide a PDB filename as an argument"
    echo "Usage: sbatch prepare_structure.sh <structure_name.pdb>"
    exit 1
fi

# Strip file extension from the input filename
INPUT_NAME=$1
STRUCTURE_NAME="${INPUT_NAME%.*}"
BASE_DIR="/fs1/home/pi/kuhfeldr"
MD_DIR="${BASE_DIR}/md"
# Handle both relative and absolute paths
if [[ "${INPUT_NAME}" == /* ]]; then
    PDB_FILE="${INPUT_NAME}"
else
    PDB_FILE="${MD_DIR}/${INPUT_NAME}"
fi
PREP_DIR="${MD_DIR}/${STRUCTURE_NAME}"

# Create preparation directory first
mkdir -p "${PREP_DIR}"

# Change to preparation directory so SLURM outputs go here
cd "${PREP_DIR}"

# Copy script to preparation directory
SCRIPT_NAME=$(basename "$0")
SCRIPT_PATH="${MD_DIR}/${SCRIPT_NAME}"
if [ -f "${SCRIPT_PATH}" ]; then
    cp "${SCRIPT_PATH}" "${PREP_DIR}/"
fi

# OpenCL environment for GTX 1080 Ti
export GPU_DEVICE_ORDINAL=0
export OCL_ICD_VENDORS=/etc/OpenCL/vendors
export OPENCL_VENDOR_PATH=/etc/OpenCL/vendors

echo "=========================================="
echo "GROMACS STRUCTURE PREPARATION"
echo "Structure: ${STRUCTURE_NAME}"
echo "Job ID: ${SLURM_JOB_ID}"
echo "=========================================="
echo "Job started: $(date)"
echo "Node: $(hostname)"
echo "Working directory: ${PREP_DIR}"
echo ""

# Check if PDB file exists
if [ ! -f "${PDB_FILE}" ]; then
    echo "ERROR: PDB file not found: ${PDB_FILE}"
    exit 1
fi

# Copy all .mdp files to preparation directory
echo "Copying input files..."
cp "${MD_DIR}"/*.mdp . 2>/dev/null || true
PDB_BASENAME=$(basename "${PDB_FILE}")
cp "${PDB_FILE}" .
echo "✓ Input files copied"
echo ""

# Already in preparation directory (changed above)

# Step 1: Generate topology from PDB with AMBER99SB-ILDN forcefield
echo "=========================================="
echo "Step 1: Generating topology (pdb2gmx)"
echo "=========================================="
echo "Forcefield: AMBER99SB-ILDN (option 1)"
echo "Water model: TIP3P (option 1)"
echo ""

# Use AMBER99SB-ILDN (option 1) and TIP3P water (option 1)
PDB_BASENAME=$(basename "${PDB_FILE}")
echo -e "1\n1" | gmx pdb2gmx -f "${PDB_BASENAME}" -o processed.gro -water tip3p -ignh

if [ $? -ne 0 ]; then
    echo "ERROR: pdb2gmx failed"
    exit 1
fi

echo "✓ Topology generated with AMBER99SB-ILDN + TIP3P"
echo ""

# Step 2: Define the simulation box - CUBIC with 1.2 nm buffer
echo "=========================================="
echo "Step 2: Creating simulation box (editconf)"
echo "=========================================="
echo "Box type: Cubic"
echo "Buffer: 1.2 nm"
echo ""

gmx editconf -f processed.gro -o boxed.gro -c -d 1.2 -bt cubic

if [ $? -ne 0 ]; then
    echo "ERROR: editconf failed"
    exit 1
fi

echo "✓ Cubic box created with 1.2 nm buffer"
echo ""

# Step 3: Solvate the system with TIP3P water
echo "=========================================="
echo "Step 3: Solvating system (solvate)"
echo "=========================================="
echo ""

gmx solvate -cp boxed.gro -cs spc216.gro -o solvated.gro -p topol.top

if [ $? -ne 0 ]; then
    echo "ERROR: solvate failed"
    exit 1
fi

NWATERS=$(grep "SOL" topol.top | tail -1 | awk '{print $2}' 2>/dev/null || echo "unknown")
echo "✓ Added $NWATERS water molecules"
echo ""

# Step 4: Add ions - Na:Cl at 0.15 M concentration (1:1 ratio)
echo "=========================================="
echo "Step 4: Adding ions (genion)"
echo "=========================================="
echo "Concentration: 0.15 M NaCl (1:1 ratio)"
echo ""

gmx grompp -f ions.mdp -c solvated.gro -p topol.top -o ions.tpr -maxwarn 1

if [ $? -ne 0 ]; then
    echo "ERROR: grompp for ions failed"
    exit 1
fi

# Add ions: neutralize + 0.15 M concentration
echo "SOL" | gmx genion -s ions.tpr -o ionized.gro -p topol.top \
    -pname NA -nname CL -neutral -conc 0.15

if [ $? -ne 0 ]; then
    echo "ERROR: genion failed"
    exit 1
fi

NNA=$(grep "NA" topol.top | tail -1 | awk '{print $2}' 2>/dev/null || echo "unknown")
NCL=$(grep "CL" topol.top | tail -1 | awk '{print $2}' 2>/dev/null || echo "unknown")
echo "✓ Added $NNA Na+ and $NCL Cl- ions"
echo ""

# Step 5: Energy minimization
echo "=========================================="
echo "Step 5: Energy minimization"
echo "=========================================="
echo ""

gmx grompp -f minim.mdp -c ionized.gro -p topol.top -o em.tpr -maxwarn 1

if [ $? -ne 0 ]; then
    echo "ERROR: grompp for minimization failed"
    exit 1
fi

gmx mdrun -v -deffnm em -ntmpi 1 -ntomp ${SLURM_CPUS_PER_TASK:-4}

if [ $? -ne 0 ]; then
    echo "ERROR: Energy minimization failed"
    exit 1
fi

EPOT=$(grep "Potential Energy" em.log | tail -1 | awk '{print $NF}' 2>/dev/null || echo "N/A")
echo "✓ Energy minimization complete. Final Epot: $EPOT"
echo ""

# Step 6: NVT equilibration at 300K
echo "=========================================="
echo "Step 6: NVT equilibration (300K)"
echo "=========================================="
echo "Duration: 100 ps"
echo ""

gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr -maxwarn 1

if [ $? -ne 0 ]; then
    echo "ERROR: grompp for NVT failed"
    exit 1
fi

gmx mdrun -deffnm nvt -ntmpi 1 -ntomp ${SLURM_CPUS_PER_TASK:-4}

if [ $? -ne 0 ]; then
    echo "ERROR: NVT equilibration failed"
    exit 1
fi

echo "✓ NVT equilibration complete"
echo ""

# Step 7: NPT equilibration at 300K, 1 bar
echo "=========================================="
echo "Step 7: NPT equilibration (300K, 1 bar)"
echo "=========================================="
echo "Duration: 100 ps"
echo ""

gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr -maxwarn 1

if [ $? -ne 0 ]; then
    echo "ERROR: grompp for NPT failed"
    exit 1
fi

gmx mdrun -deffnm npt -ntmpi 1 -ntomp ${SLURM_CPUS_PER_TASK:-4}

if [ $? -ne 0 ]; then
    echo "ERROR: NPT equilibration failed"
    exit 1
fi

echo "✓ NPT equilibration complete"
echo ""

# Step 8: Generate production MD .tpr file
echo "=========================================="
echo "Step 8: Generating production MD .tpr file"
echo "=========================================="
echo "Duration: 100 ns"
echo ""

gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr -maxwarn 1

if [ $? -ne 0 ]; then
    echo "ERROR: grompp for production MD failed"
    exit 1
fi

echo "✓ Production MD .tpr file generated: md.tpr"
echo ""

# Summary
echo "=========================================="
echo "PREPARATION COMPLETE!"
echo "=========================================="
echo "Structure: ${STRUCTURE_NAME}"
echo "Job completed: $(date)"
echo ""
echo "Output directory: ${PREP_DIR}"
echo ""
echo "Generated files:"
echo "  - md.tpr        : Production MD input file (READY FOR MD RUN)"
echo "  - npt.gro       : Final equilibrated structure"
echo "  - npt.cpt       : Checkpoint file"
echo "  - topol.top     : System topology"
echo "  - posre_*.itp   : Position restraint files"
echo ""
echo "Next step: Run production MD with:"
echo "  sbatch run_md_1080ti.sh ${STRUCTURE_NAME}.pdb"
echo "=========================================="

