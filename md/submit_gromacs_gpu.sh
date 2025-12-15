#!/bin/bash
#SBATCH --job-name=md_100ns
#SBATCH --output=md_production_%j.out
#SBATCH --error=md_production_%j.err
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --gres=gpu:l40s:4
#SBATCH --mem=128G
#SBATCH --time=1-00:00:00
#
# =============================================================================
# OPTIMIZED FOR: Production MD on orcaga normal partition (L40S GPUs)
# =============================================================================
# Using: 4x NVIDIA L40S GPUs (48GB each) on orcaga[01-05]
# Expected performance: ~50-80 ns/day (L40S is ~3x faster than A30)
# Time allocation: 1 day (normal partition limit)
# NOTE: 100 ns runs will checkpoint and can be continued in subsequent jobs
#
# GPU OPTIONS available on orcaga cluster:
# -----------------------------------------
# Long partition (7 days) - A30 GPUs on orcaga[11-19]:
#   --gres=gpu:a30:4   (4 GPUs - RECOMMENDED for 100ns runs, no checkpoint needed)
#
# Normal partition (1 day) - L40S GPUs on orcaga[01-05]:
#   --gres=gpu:l40s:4  (4 GPUs - faster but 1-day limit, checkpoint for 100ns)
#
# CPU/Memory: 32 CPUs + 128GB is optimal for 4-GPU MD runs
# (8 threads per GPU, leaves headroom for other processes)

# GPU-Accelerated GROMACS MD Simulation
# =====================================
# Using: GROMACS 2025.3 via local Apptainer container
# Target: 4x NVIDIA L40S GPUs (48GB VRAM each) on orcaga[01-05]
# Container: /home/kuhfeldr/gromacs_v2025.3.sif
#
# Features:
# - Automatic checkpoint detection (skips completed steps)
# - Full GPU acceleration: NB, PME, bonded, update all on GPU
# - Multi-GPU domain decomposition with direct GPU communication
# - Checkpoint every 4 hours for fault tolerance
#
# Simulation Settings:
# - Forcefield: AMBER99SB-ILDN (robust for protein-peptide)
# - Temperature: 300K (NVT → NPT equilibration done)
# - Solvent: TIP3P water in cubic box (1.2 nm buffer)
# - Salt: 0.15 M NaCl (physiological)
# - Duration: 100 ns production
# - Output: Trajectory includes all atoms

# =============================================================================
# SAFETY CHECK - Ensure running via sbatch on a compute node
# =============================================================================
if [ -z "$SLURM_JOB_ID" ]; then
    echo "ERROR: This script must be submitted via sbatch, not run directly!"
    echo ""
    echo "Usage:"
    echo "  cd /path/to/simulation/directory"
    echo "  sbatch submit_gromacs_gpu.sh"
    echo ""
    echo "This will queue the job to run on a GPU compute node."
    exit 1
fi

# =============================================================================
# LOAD REQUIRED MODULES
# =============================================================================
module load apptainer/1.4.1-gcc-13.4.0

# =============================================================================
# CONTAINER SETUP
# =============================================================================
CONTAINER_PATH="/home/kuhfeldr/gromacs_v2025.3.sif"

# Verify container exists
if [ ! -f "$CONTAINER_PATH" ]; then
    echo "ERROR: GROMACS container not found at: $CONTAINER_PATH"
    echo ""
    echo "Please ensure the container file exists. You can download it with:"
    echo "  apptainer pull gromacs_v2025.3.sif docker://nvcr.io/hpc/gromacs:2025.3"
    exit 1
fi

# Set default CPU count (32 CPUs optimal for 4 GPUs = 8 threads/GPU)
CPUS=${SLURM_CPUS_PER_TASK:-32}

# Define apptainer command with GPU support
SINGULARITY="apptainer run --nv -B ${PWD}:/host_pwd --pwd /host_pwd ${CONTAINER_PATH}"

# Wrapper function for GMX commands
run_gmx() {
    ${SINGULARITY} gmx "$@"
}

# =============================================================================
# SYSTEM NAMING HELPERS (for self-describing outputs)
# =============================================================================
# Many downstream tools need to clearly distinguish peptide vs receptor.
# GROMACS default index groups in the .tpr often do NOT preserve chain IDs,
# so we also create an explicit analysis index at the end of the run.
SYS_DIR_BASENAME="$(basename "$(pwd)")"
PEPTIDE_SEQ="$(echo "${SYS_DIR_BASENAME}" | sed -E 's/^([A-Z]+)_.*/\1/')"
if [[ ! "${PEPTIDE_SEQ}" =~ ^[A-Z]+$ ]]; then
    PEPTIDE_SEQ="UNKNOWN"
fi
PEPTIDE_LEN="${#PEPTIDE_SEQ}"
SYS_TAG="${PEPTIDE_SEQ}"

echo "=========================================="
echo "GROMACS 2025.3 - GPU SIMULATION"
echo "=========================================="
echo "Job started at: $(date)"
echo "Running on node: $(hostname)"
echo "SLURM Job ID: $SLURM_JOB_ID"
echo "Working directory: $(pwd)"
echo "Container: $CONTAINER_PATH"
echo "CPUs per task: $SLURM_CPUS_PER_TASK"
echo "System tag: ${SYS_TAG}"
echo "Peptide sequence: ${PEPTIDE_SEQ} (len=${PEPTIDE_LEN})"
echo ""

# Enable direct GPU-to-GPU communication for multi-GPU runs
export GMX_ENABLE_DIRECT_GPU_COMM=1
export GMX_GPU_DD_COMMS=true
export GMX_GPU_PME_PP_COMMS=true
export GMX_FORCE_UPDATE_DEFAULT_GPU=true

# =============================================================================
# GPU DETECTION
# =============================================================================
echo "GPU Information:"
if command -v nvidia-smi &> /dev/null; then
    nvidia-smi -L
    echo ""
fi

if [ -n "$CUDA_VISIBLE_DEVICES" ]; then
    NUM_GPUS=$(echo "$CUDA_VISIBLE_DEVICES" | tr ',' '\n' | wc -l)
    echo "SLURM allocated $NUM_GPUS GPU(s): $CUDA_VISIBLE_DEVICES"
elif command -v nvidia-smi &> /dev/null; then
    NUM_GPUS=$(nvidia-smi -L 2>/dev/null | wc -l)
    echo "Detected $NUM_GPUS GPU(s) via nvidia-smi"
else
    NUM_GPUS=1
    echo "Assuming 1 GPU (could not detect)"
fi

# Build GPU ID string
GPU_IDS=""
for ((i=0; i<NUM_GPUS; i++)); do
    GPU_IDS="${GPU_IDS}${i}"
done
echo "GPU IDs for GROMACS: $GPU_IDS"
echo ""

# =============================================================================
# CHECKPOINT DETECTION - Skip completed steps
# =============================================================================
echo "Checking for existing files (checkpoint detection)..."
echo ""

START_STEP=1

if [ -f "md.gro" ] && [ -f "md.xtc" ]; then
    echo "✓ Production MD already complete (md.gro and md.xtc exist)"
    echo "  Nothing to do. Exiting."
    exit 0
elif [ -f "npt.gro" ] && [ -f "npt.cpt" ] && [ -f "topol.top" ]; then
    echo "✓ Found: npt.gro, npt.cpt, topol.top"
    echo "  → Skipping Steps 1-7 (preparation complete)"
    START_STEP=8
elif [ -f "nvt.gro" ] && [ -f "nvt.cpt" ] && [ -f "topol.top" ]; then
    echo "✓ Found: nvt.gro, nvt.cpt, topol.top"
    echo "  → Starting at Step 7 (NPT equilibration)"
    START_STEP=7
elif [ -f "em.gro" ] && [ -f "topol.top" ]; then
    echo "✓ Found: em.gro, topol.top"
    echo "  → Starting at Step 6 (NVT equilibration)"
    START_STEP=6
elif [ -f "ionized.gro" ] && [ -f "topol.top" ]; then
    echo "✓ Found: ionized.gro, topol.top"
    echo "  → Starting at Step 5 (Energy minimization)"
    START_STEP=5
elif [ -f "solvated.gro" ] && [ -f "topol.top" ]; then
    echo "✓ Found: solvated.gro, topol.top"
    echo "  → Starting at Step 4 (Add ions)"
    START_STEP=4
elif [ -f "boxed.gro" ] && [ -f "topol.top" ]; then
    echo "✓ Found: boxed.gro, topol.top"
    echo "  → Starting at Step 3 (Solvation)"
    START_STEP=3
elif [ -f "processed.gro" ] && [ -f "topol.top" ]; then
    echo "✓ Found: processed.gro, topol.top"
    echo "  → Starting at Step 2 (Box creation)"
    START_STEP=2
else
    echo "  No checkpoint files found. Starting from Step 1."
    START_STEP=1
fi

echo ""
echo "=========================================="
echo "Starting from Step $START_STEP"
echo "=========================================="
echo ""

# =============================================================================
# STEP 1: Generate topology with AMBER99SB-ILDN forcefield
# =============================================================================
if [ $START_STEP -le 1 ]; then
    echo "Step 1: Running pdb2gmx with AMBER99SB-ILDN forcefield..."
    
    # Find PDB file
    PDB_FILE=$(ls *.pdb 2>/dev/null | head -1)
    if [ -z "$PDB_FILE" ]; then
        echo "ERROR: No PDB file found in current directory"
        exit 1
    fi
    echo "Using PDB file: $PDB_FILE"
    
    # Best practice for peptide–protein complexes: keep chains as separate molecule types
    # so downstream selection can use moltype (e.g., Protein_chain_P) instead of residue ranges.
    #
    # -chainsep id_or_ter: split molecules when chain ID changes (and/or TER)
    # -merge no          : do NOT merge chains into one moleculetype
    echo -e "1\n1" | run_gmx pdb2gmx \
        -f "$PDB_FILE" \
        -o processed.gro \
        -water tip3p \
        -ignh \
        -chainsep id_or_ter \
        -merge no

    if [ $? -ne 0 ]; then
        echo "ERROR: pdb2gmx failed"
        exit 1
    fi
    echo "✓ Topology generated"
    echo ""
else
    echo "Step 1: SKIPPED (processed.gro exists)"
fi

# =============================================================================
# STEP 2: Define the simulation box
# =============================================================================
if [ $START_STEP -le 2 ]; then
    echo "Step 2: Creating cubic simulation box (1.2 nm buffer)..."
    
    run_gmx editconf -f processed.gro -o boxed.gro -c -d 1.2 -bt cubic

    if [ $? -ne 0 ]; then
        echo "ERROR: editconf failed"
        exit 1
    fi
    echo "✓ Cubic box created"
    echo ""
else
    echo "Step 2: SKIPPED (boxed.gro exists)"
fi

# =============================================================================
# STEP 3: Solvate the system
# =============================================================================
if [ $START_STEP -le 3 ]; then
    echo "Step 3: Solvating system with TIP3P water..."
    
    run_gmx solvate -cp boxed.gro -cs spc216.gro -o solvated.gro -p topol.top

    if [ $? -ne 0 ]; then
        echo "ERROR: solvate failed"
        exit 1
    fi
    NWATERS=$(grep "SOL" topol.top | tail -1 | awk '{print $2}' 2>/dev/null || echo "unknown")
    echo "✓ Added $NWATERS water molecules"
    echo ""
else
    echo "Step 3: SKIPPED (solvated.gro exists)"
fi

# =============================================================================
# STEP 4: Add ions
# =============================================================================
if [ $START_STEP -le 4 ]; then
    echo "Step 4: Adding ions (0.15 M NaCl)..."
    
    run_gmx grompp -f ions.mdp -c solvated.gro -p topol.top -o ions.tpr -maxwarn 1
    echo "SOL" | run_gmx genion -s ions.tpr -o ionized.gro -p topol.top \
        -pname NA -nname CL -neutral -conc 0.15

    if [ $? -ne 0 ]; then
        echo "ERROR: genion failed"
        exit 1
    fi
    NNA=$(grep "NA" topol.top | tail -1 | awk '{print $2}' 2>/dev/null || echo "unknown")
    NCL=$(grep "CL" topol.top | tail -1 | awk '{print $2}' 2>/dev/null || echo "unknown")
    echo "✓ Added $NNA Na+ and $NCL Cl- ions"
    echo ""
else
    echo "Step 4: SKIPPED (ionized.gro exists)"
fi

# =============================================================================
# STEP 5: Energy minimization
# =============================================================================
if [ $START_STEP -le 5 ]; then
    echo "Step 5: Running energy minimization..."
    
    run_gmx grompp -f minim.mdp -c ionized.gro -p topol.top -o em.tpr -maxwarn 1
    run_gmx mdrun -v -deffnm em -ntmpi 1 -ntomp $CPUS -nb gpu

    if [ $? -ne 0 ]; then
        echo "ERROR: Energy minimization failed"
        exit 1
    fi
    echo "✓ Energy minimization complete"
    echo ""
else
    echo "Step 5: SKIPPED (em.gro exists)"
fi

# =============================================================================
# STEP 6: NVT equilibration
# =============================================================================
if [ $START_STEP -le 6 ]; then
    echo "Step 6: Running NVT equilibration (300K, 250 ps)..."
    
    run_gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr -maxwarn 1
    run_gmx mdrun -deffnm nvt -ntmpi 1 -ntomp $CPUS -nb gpu -pme gpu -update gpu -v

    if [ $? -ne 0 ]; then
        echo "ERROR: NVT equilibration failed"
        exit 1
    fi
    echo "✓ NVT equilibration complete (250 ps)"
    echo ""
else
    echo "Step 6: SKIPPED (nvt.gro exists)"
fi

# =============================================================================
# STEP 7: NPT equilibration
# =============================================================================
if [ $START_STEP -le 7 ]; then
    echo "Step 7: Running NPT equilibration (300K, 1 bar, 500 ps)..."
    
    run_gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr -maxwarn 1
    run_gmx mdrun -deffnm npt -ntmpi 1 -ntomp $CPUS -nb gpu -pme gpu -update gpu -v

    if [ $? -ne 0 ]; then
        echo "ERROR: NPT equilibration failed"
        exit 1
    fi
    echo "✓ NPT equilibration complete (500 ps)"
    echo ""
else
    echo "Step 7: SKIPPED (npt.gro exists)"
fi

# =============================================================================
# STEP 8: Production MD (100 ns)
# =============================================================================
echo ""
echo "=========================================="
echo "Step 8: Running production MD (100 ns)..."
echo "=========================================="
echo "Settings:"
echo "  - Temperature: 300K"
echo "  - Pressure: 1 bar (NPT ensemble)"
echo "  - Duration: 100 ns"
echo "  - GPUs: ${NUM_GPUS}x L40S (multi-GPU domain decomposition)"
echo "  - Checkpoint: Every 60 min (for fault tolerance)"
echo "  - Expected: ~50-80 ns/day with 4x L40S"
echo ""

# Check for checkpoint or generate tpr
if [ -f "md.cpt" ]; then
    echo "Found md.cpt - continuing from checkpoint..."
    CONTINUE_FLAG="-cpi md.cpt"
else
    echo "No checkpoint found - starting fresh production run..."
    run_gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr -maxwarn 1
    if [ $? -ne 0 ]; then
        echo "ERROR: grompp for production MD failed"
        exit 1
    fi
    CONTINUE_FLAG=""
fi

# Calculate threads per rank for multi-GPU
THREADS_PER_RANK=$(($CPUS / $NUM_GPUS))
[ $THREADS_PER_RANK -lt 1 ] && THREADS_PER_RANK=1

echo "Launching MD run..."
echo "  MPI ranks: $NUM_GPUS"
echo "  Threads/rank: $THREADS_PER_RANK"
echo ""

# Run production MD with GPU acceleration
if [ $NUM_GPUS -gt 1 ]; then
    # Multi-GPU run (optimized for 4x L40S)
    # - npme 1: Dedicate 1 GPU to PME for better load balance
    # - nstlist 200: Good balance for GPU neighbor searching
    # - cpt 60: Checkpoint every hour for fault tolerance
    run_gmx mdrun -deffnm md \
        ${CONTINUE_FLAG} \
        -ntmpi ${NUM_GPUS} \
        -ntomp ${THREADS_PER_RANK} \
        -nb gpu \
        -pme gpu \
        -npme 1 \
        -update gpu \
        -bonded gpu \
        -nstlist 200 \
        -pin on \
        -gpu_id ${GPU_IDS} \
        -cpt 60 \
        -v
else
    # Single-GPU run (optimized)
    run_gmx mdrun -deffnm md \
        ${CONTINUE_FLAG} \
        -ntmpi 1 \
        -ntomp $CPUS \
        -nb gpu \
        -pme gpu \
        -update gpu \
        -bonded gpu \
        -nstlist 200 \
        -cpt 60 \
        -v
fi

# Check if MD completed successfully
if [ $? -ne 0 ]; then
    echo ""
    echo "WARNING: GPU run failed. Trying CPU fallback..."
    run_gmx mdrun -deffnm md ${CONTINUE_FLAG} -ntomp $CPUS -cpt 60 -v
fi

# =============================================================================
# SUMMARY
# =============================================================================
echo ""
echo "=========================================="
echo "SIMULATION COMPLETE!"
echo "=========================================="
echo "Job completed at: $(date)"
echo ""
echo "Output files:"
echo "  - md.xtc       : Trajectory (100 ns)"
echo "  - md.edr       : Energy file"
echo "  - md.log       : Detailed log"
echo "  - md.gro       : Final structure"
echo "  - md.cpt       : Checkpoint (for restart)"
echo "  - topol.top    : System topology"
echo ""
echo "System configuration:"
echo "  - Container    : GROMACS 2025.3"
echo "  - Forcefield   : AMBER99SB-ILDN"
echo "  - Water model  : TIP3P"
echo "  - Temperature  : 300 K"
echo "  - Pressure     : 1 bar"
echo "  - Salt         : 0.15 M NaCl"
echo "  - Box          : Cubic (1.2 nm buffer)"
echo "  - Duration     : 100 ns"
echo "  - GPUs         : ${NUM_GPUS}x L40S GPU(s)"
echo "  - Checkpoints  : Every 60 minutes"
echo ""

if [ -f md.log ]; then
    echo "Performance summary:"
    grep -E "Performance|ns/day|hour/ns" md.log | tail -5
fi

echo ""
echo "=========================================="
echo "Post-processing: labeling peptide vs receptor"
echo "=========================================="
#
# Create explicit, named index groups for downstream analysis:
#   [ Receptor ]        residues 1..REC_END
#   [ Peptide ]         residues PEP_START..TOTAL_PROT_RES
#   [ PeptideBackbone ] peptide backbone atoms
#
# Prefer moltype-based selections (robust when pdb2gmx uses -merge no and chain separation),
# and fall back to residue-range selection only if moltypes aren't available.
#
if [ -f "md.tpr" ] && [ -f "md.xtc" ] && [ "${PEPTIDE_SEQ}" != "UNKNOWN" ] && [ "${PEPTIDE_LEN}" -gt 0 ]; then
    NDX_OK=false

    # 1) Moltype-based (preferred)
    if [ -f "topol_Protein_chain_P.itp" ]; then
        receptor_sel=""
        for itp in topol_Protein_chain_*.itp; do
            [ -f "${itp}" ] || continue
            chain_letter="$(echo "${itp}" | sed -E 's/.*chain_(.)\.itp/\1/')"
            moltype="Protein_chain_${chain_letter}"
            if [ "${chain_letter}" = "P" ]; then
                continue
            fi
            if [ -z "${receptor_sel}" ]; then
                receptor_sel="moltype \"${moltype}\""
            else
                receptor_sel="${receptor_sel} or moltype \"${moltype}\""
            fi
        done

        if [ -n "${receptor_sel}" ]; then
            cat > analysis_selections.txt << EOF
Receptor = ${receptor_sel}
Peptide = moltype "Protein_chain_P"
PeptideBackbone = moltype "Protein_chain_P" and backbone
EOF
            if run_gmx select -s md.tpr -sf analysis_selections.txt -on analysis.ndx >/dev/null 2>&1; then
                echo "✓ Wrote analysis.ndx via moltype selections"
                NDX_OK=true
            fi
        fi
    fi

    # 2) Residue-range fallback
    if [ "${NDX_OK}" != true ]; then
        TOTAL_PROT_RES="$(printf "q\n" | run_gmx make_ndx -f md.tpr -o /dev/null 2>/dev/null | awk '/There are:/{if($3=="Protein") {print $2; exit}}')"
        if [[ "${TOTAL_PROT_RES}" =~ ^[0-9]+$ ]] && [ "${TOTAL_PROT_RES}" -gt "${PEPTIDE_LEN}" ]; then
            REC_END="$((TOTAL_PROT_RES - PEPTIDE_LEN))"
            PEP_START="$((REC_END + 1))"
            echo "Protein residues (total): ${TOTAL_PROT_RES}"
            echo "Receptor residues: 1-${REC_END}"
            echo "Peptide residues: ${PEP_START}-${TOTAL_PROT_RES}"

            cat > analysis_selections.txt << EOF
Receptor = resnr 1 to ${REC_END}
Peptide = resnr ${PEP_START} to ${TOTAL_PROT_RES}
PeptideBackbone = (resnr ${PEP_START} to ${TOTAL_PROT_RES}) and backbone
EOF
            if run_gmx select -s md.tpr -sf analysis_selections.txt -on analysis.ndx >/dev/null 2>&1; then
                echo "✓ Wrote analysis.ndx via residue-range fallback"
                NDX_OK=true
            fi
        fi
    fi

    if [ "${NDX_OK}" != true ]; then
        echo "WARNING: Failed to write analysis.ndx (moltype and residue-range methods failed)"
    fi
else
    echo "WARNING: Missing md.tpr/md.xtc or unknown peptide sequence; skipping analysis.ndx generation."
fi

# Tag/copy key outputs so it's obvious which peptide they belong to (non-destructive).
for f in md.xtc md.tpr md.gro md.edr md.log md.cpt; do
    if [ -f "${f}" ]; then
        cp -n "${f}" "${SYS_TAG}_${f}" 2>/dev/null || true
    fi
done
if [ -f "analysis.ndx" ]; then
    cp -n "analysis.ndx" "${SYS_TAG}_analysis.ndx" 2>/dev/null || true
fi

echo ""
echo "=========================================="
echo "Next steps:"
echo "  1. Analyze: ./analyze_md.sh $(pwd)"
echo "  2. Visualize: Download md.xtc for PyMOL/VMD"
echo "=========================================="
