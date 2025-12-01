#!/bin/bash
#SBATCH --job-name=gromacs_gpu
#SBATCH --output=gromacs_gpu_%j.out
#SBATCH --error=gromacs_gpu_%j.err
#SBATCH --partition=cqls_gpu
#SBATCH --account=cqls
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:v100:5
# Note: Script automatically uses ALL allocated GPUs (1, 2, 3, 4, or 5)
#       Change the number above based on availability:
#       --gres=gpu:v100:1  (single GPU, fastest to start)
#       --gres=gpu:v100:2  (2 GPUs, good balance)
#       --gres=gpu:v100:4  (4 GPUs, high performance)
#       --gres=gpu:v100:5  (5 GPUs, maximum performance on cqls-gpu1)
#SBATCH --mem=64G
#SBATCH --time=24:00:00

# GPU-Accelerated GROMACS MD Simulation using NVIDIA NGC Container
# ==================================================================
# Using: NVIDIA NGC GROMACS 2023.2 via Singularity
# Target: V100 GPUs (5x) with optimized multi-GPU settings
# Reference: https://catalog.ngc.nvidia.com/orgs/hpc/containers/gromacs
#
# Features:
# - Automatic checkpoint detection (skips completed steps)
# - Multi-GPU support with direct GPU communication
#
# Settings:
# - Forcefield: AMBER99SB-ILDN
# - Temperature: 300K
# - Solvent: TIP3P water in cubic box (1.2 nm buffer)
# - Salt: 0.15 M NaCl (1:1 ratio)
# - Duration: 100 ns
# - Waters included in trajectory

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

echo "=========================================="
echo "NGC GROMACS 2023.2 - MULTI-GPU SIMULATION"
echo "=========================================="
echo "Job started at: $(date)"
echo "Running on node: $(hostname)"
echo "SLURM Job ID: $SLURM_JOB_ID"
echo "GPUs allocated: $CUDA_VISIBLE_DEVICES"
echo "CPUs per task: $SLURM_CPUS_PER_TASK"
echo ""

# =============================================================================
# SINGULARITY/NGC CONTAINER SETUP
# =============================================================================
GROMACS_TAG="2023.2"

# Set default CPU count if not set
CPUS=${SLURM_CPUS_PER_TASK:-16}

# Detect architecture - NGC container only supports x86_64
ARCH=$(uname -m)
if [ "$ARCH" = "ppc64le" ]; then
    echo "ERROR: NGC GROMACS container requires x86_64 architecture!"
    echo ""
    echo "Detected: Power9 (ppc64le) - NGC container not available for this architecture"
    echo ""
    echo "Solution: Submit to x86_64 partition (cqls_gpu) instead of cqls_ppc64le"
    echo ""
    exit 1
fi

# Enable direct GPU-to-GPU communication for multi-GPU runs
export GMX_ENABLE_DIRECT_GPU_COMM=1

# Force GPU-aware MPI (if available)
export GMX_GPU_DD_COMMS=true
export GMX_GPU_PME_PP_COMMS=true
export GMX_FORCE_UPDATE_DEFAULT_GPU=true

# Define Singularity command with GPU support
SINGULARITY="singularity run --nv -B ${PWD}:/host_pwd --pwd /host_pwd docker://nvcr.io/hpc/gromacs:${GROMACS_TAG}"

# Wrapper function for GMX commands via Singularity
run_gmx() {
    ${SINGULARITY} gmx "$@"
}

echo "Container: nvcr.io/hpc/gromacs:${GROMACS_TAG}"
echo "Singularity command: $SINGULARITY"
echo ""

# Check GPU availability
echo "GPU Information:"
if command -v nvidia-smi &> /dev/null; then
    nvidia-smi
else
    echo "WARNING: nvidia-smi not found!"
fi
echo ""

# =============================================================================
# GPU ALLOCATION DETECTION
# =============================================================================
# Use CUDA_VISIBLE_DEVICES (set by SLURM) to determine allocated GPUs
# This is more reliable than counting all physical GPUs on the node
#
# Example: If SLURM allocates GPUs 0,1,2,4 on a 5-GPU node:
#   CUDA_VISIBLE_DEVICES="0,1,2,4"
#   NUM_GPUS=4
#   GPU_IDS="0123" (CUDA remaps to 0-based sequential IDs)

if [ -n "$CUDA_VISIBLE_DEVICES" ]; then
    # Count GPUs from CUDA_VISIBLE_DEVICES (comma-separated list)
    NUM_GPUS=$(echo "$CUDA_VISIBLE_DEVICES" | tr ',' '\n' | wc -l)
    echo "SLURM allocated $NUM_GPUS GPU(s) via CUDA_VISIBLE_DEVICES: $CUDA_VISIBLE_DEVICES"
elif command -v nvidia-smi &> /dev/null; then
    # Fallback: count all visible GPUs
    NUM_GPUS=$(nvidia-smi -L 2>/dev/null | wc -l)
    echo "Detected $NUM_GPUS GPU(s) via nvidia-smi (CUDA_VISIBLE_DEVICES not set)"
else
    NUM_GPUS=0
fi

# Ensure at least 1 GPU for calculations (avoid division by zero)
if [ "$NUM_GPUS" -lt 1 ]; then
    echo "ERROR: No GPUs detected on this node!"
    echo "Make sure you requested GPUs with --gres=gpu:..."
    exit 1
fi

# Build GPU ID string using SEQUENTIAL IDs (0, 1, 2, ...)
# CUDA_VISIBLE_DEVICES remaps physical GPUs to sequential virtual IDs
# So if SLURM gives us GPUs 0,1,2,4, GROMACS sees them as 0,1,2,3
GPU_IDS=""
for ((i=0; i<NUM_GPUS; i++)); do
    GPU_IDS="${GPU_IDS}${i}"
done
echo "GPU IDs for GROMACS: $GPU_IDS (sequential, remapped by CUDA)"
echo ""

# =============================================================================
# CHECKPOINT DETECTION - Skip completed steps
# =============================================================================
echo "Checking for existing files (checkpoint detection)..."
echo ""

# Determine starting step based on existing files
START_STEP=1

# Check from the end backwards to find the latest completed step
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
    echo "  → Skipping Steps 1-6, starting at Step 7 (NPT equilibration)"
    START_STEP=7
elif [ -f "em.gro" ] && [ -f "topol.top" ]; then
    echo "✓ Found: em.gro, topol.top"
    echo "  → Skipping Steps 1-5, starting at Step 6 (NVT equilibration)"
    START_STEP=6
elif [ -f "ionized.gro" ] && [ -f "topol.top" ]; then
    echo "✓ Found: ionized.gro, topol.top"
    echo "  → Skipping Steps 1-4, starting at Step 5 (Energy minimization)"
    START_STEP=5
elif [ -f "solvated.gro" ] && [ -f "topol.top" ]; then
    echo "✓ Found: solvated.gro, topol.top"
    echo "  → Skipping Steps 1-3, starting at Step 4 (Add ions)"
    START_STEP=4
elif [ -f "boxed.gro" ] && [ -f "topol.top" ]; then
    echo "✓ Found: boxed.gro, topol.top"
    echo "  → Skipping Steps 1-2, starting at Step 3 (Solvation)"
    START_STEP=3
elif [ -f "processed.gro" ] && [ -f "topol.top" ]; then
    echo "✓ Found: processed.gro, topol.top"
    echo "  → Skipping Step 1, starting at Step 2 (Box creation)"
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
    echo "Forcefield: AMBER99SB-ILDN (option 1)"
    echo "Water model: TIP3P (option 1)"
    echo ""

    # Use AMBER99SB-ILDN (option 1) and TIP3P water (option 1)
    echo -e "1\n1" | run_gmx pdb2gmx -f 3fxi.pdb -o processed.gro -water tip3p -ignh

    if [ $? -ne 0 ]; then
        echo "Error in pdb2gmx step"
        exit 1
    fi

    echo "✓ Topology generated with AMBER99SB-ILDN + TIP3P"
    echo ""
else
    echo "Step 1: SKIPPED (processed.gro exists)"
fi

# =============================================================================
# STEP 2: Define the simulation box - CUBIC with 1.2 nm buffer
# =============================================================================
if [ $START_STEP -le 2 ]; then
    echo "Step 2: Creating cubic simulation box (1.2 nm buffer)..."
    echo ""

    run_gmx editconf -f processed.gro -o boxed.gro -c -d 1.2 -bt cubic

    if [ $? -ne 0 ]; then
        echo "Error in editconf step"
        exit 1
    fi

    echo "✓ Cubic box created - protein fully encapsulated"
    echo ""
else
    echo "Step 2: SKIPPED (boxed.gro exists)"
fi

# =============================================================================
# STEP 3: Solvate the system with TIP3P water
# =============================================================================
if [ $START_STEP -le 3 ]; then
    echo "Step 3: Solvating system with TIP3P water..."
    echo ""

    run_gmx solvate -cp boxed.gro -cs spc216.gro -o solvated.gro -p topol.top

    if [ $? -ne 0 ]; then
        echo "Error in solvate step"
        exit 1
    fi

    NWATERS=$(grep "SOL" topol.top | tail -1 | awk '{print $2}')
    echo "✓ Added $NWATERS water molecules"
    echo ""
else
    echo "Step 3: SKIPPED (solvated.gro exists)"
fi

# =============================================================================
# STEP 4: Add ions - Na:Cl at 0.15 M concentration (1:1 ratio)
# =============================================================================
if [ $START_STEP -le 4 ]; then
    echo "Step 4: Adding ions (0.15 M NaCl, 1:1 ratio)..."
    echo ""

    run_gmx grompp -f ions.mdp -c solvated.gro -p topol.top -o ions.tpr -maxwarn 1

    echo "SOL" | run_gmx genion -s ions.tpr -o ionized.gro -p topol.top \
        -pname NA -nname CL -neutral -conc 0.15

    if [ $? -ne 0 ]; then
        echo "Error in genion step"
        exit 1
    fi

    NNA=$(grep "NA" topol.top | tail -1 | awk '{print $2}')
    NCL=$(grep "CL" topol.top | tail -1 | awk '{print $2}')
    echo "✓ Added $NNA Na+ and $NCL Cl- ions"
    echo ""
else
    echo "Step 4: SKIPPED (ionized.gro exists)"
fi

# =============================================================================
# STEP 5: Energy minimization (single GPU - fast)
# =============================================================================
if [ $START_STEP -le 5 ]; then
    echo "Step 5: Running energy minimization..."
    echo ""

    run_gmx grompp -f minim.mdp -c ionized.gro -p topol.top -o em.tpr -maxwarn 1
    run_gmx mdrun -v -deffnm em -ntomp $CPUS -nb gpu

    if [ $? -ne 0 ]; then
        echo "Error in energy minimization"
        exit 1
    fi

    echo "✓ Energy minimization complete"
    echo ""
else
    echo "Step 5: SKIPPED (em.gro exists)"
fi

# =============================================================================
# STEP 6: NVT equilibration at 300K (single GPU)
# =============================================================================
if [ $START_STEP -le 6 ]; then
    echo "Step 6: Running NVT equilibration (300K)..."
    echo ""

    run_gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr -maxwarn 1
    run_gmx mdrun -deffnm nvt -ntomp $CPUS -nb gpu -pme gpu -update gpu -v

    if [ $? -ne 0 ]; then
        echo "Error in NVT equilibration"
        exit 1
    fi

    echo "✓ NVT equilibration complete"
    echo ""
else
    echo "Step 6: SKIPPED (nvt.gro exists)"
fi

# =============================================================================
# STEP 7: NPT equilibration at 300K, 1 bar (single GPU)
# =============================================================================
if [ $START_STEP -le 7 ]; then
    echo "Step 7: Running NPT equilibration (300K, 1 bar)..."
    echo ""

    run_gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr -maxwarn 1
    run_gmx mdrun -deffnm npt -ntomp $CPUS -nb gpu -pme gpu -update gpu -v

    if [ $? -ne 0 ]; then
        echo "Error in NPT equilibration"
        exit 1
    fi

    echo "✓ NPT equilibration complete"
    echo ""
else
    echo "Step 7: SKIPPED (npt.gro exists)"
fi

# =============================================================================
# STEP 8: Production MD (100 ns) - MULTI-GPU ACCELERATED with NGC optimizations
# =============================================================================
echo ""
echo "=========================================="
echo "Step 8: Running production MD (100 ns) - MULTI-GPU ACCELERATED..."
echo "=========================================="
echo "Settings:"
echo "  - Temperature: 300K"
echo "  - Duration: 100 ns"
echo "  - Waters: Included in trajectory"
echo "  - GPUs: ${NUM_GPUS}x V100 with direct GPU communication"
echo "  - Container: NGC GROMACS ${GROMACS_TAG}"
echo ""

# Check if we need to generate md.tpr or can continue from checkpoint
if [ -f "md.cpt" ]; then
    echo "Found md.cpt - continuing from checkpoint..."
    CONTINUE_FLAG="-cpi md.cpt"
else
    echo "No checkpoint found - starting fresh production run..."
    run_gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr -maxwarn 1
    CONTINUE_FLAG=""
fi

# Multi-GPU mdrun with NGC optimizations
# =======================================
# -ntmpi: Number of MPI ranks (match GPU count for multi-GPU)
# -ntomp: OpenMP threads per rank
# -nb gpu: Nonbonded on GPU
# -pme gpu: PME on GPU
# -npme 1: Dedicate 1 rank to PME (for multi-GPU)
# -update gpu: Coordinate updates on GPU
# -bonded gpu: Bonded interactions on GPU
# -nstlist 300: Neighbor list update frequency (optimized for GPU)
# -pin on: Pin threads to cores
# -gpu_id: Specify GPU IDs

echo "Launching multi-GPU MD run with NGC GROMACS..."

# Calculate threads per MPI rank
THREADS_PER_RANK=$(($CPUS / $NUM_GPUS))
if [ $THREADS_PER_RANK -lt 1 ]; then
    THREADS_PER_RANK=1
fi

echo "Command: gmx mdrun -ntmpi ${NUM_GPUS} -ntomp ${THREADS_PER_RANK} ..."
echo ""

${SINGULARITY} gmx mdrun -deffnm md \
    ${CONTINUE_FLAG} \
    -ntmpi ${NUM_GPUS} \
    -ntomp ${THREADS_PER_RANK} \
    -nb gpu \
    -pme gpu \
    -npme 1 \
    -update gpu \
    -bonded gpu \
    -nstlist 300 \
    -pin on \
    -gpu_id ${GPU_IDS} \
    -v

# Check exit status and retry with fallback if needed
if [ $? -ne 0 ]; then
    echo ""
    echo "WARNING: Multi-GPU run failed. Trying single-GPU fallback..."
    
    ${SINGULARITY} gmx mdrun -deffnm md \
        ${CONTINUE_FLAG} \
        -ntmpi 1 \
        -ntomp $CPUS \
        -nb gpu \
        -pme gpu \
        -update gpu \
        -bonded gpu \
        -nstlist 200 \
        -v
    
    if [ $? -ne 0 ]; then
        echo ""
        echo "ERROR: GPU run failed. Falling back to CPU-only..."
        ${SINGULARITY} gmx mdrun -deffnm md ${CONTINUE_FLAG} -ntomp $CPUS -v
    fi
fi

echo ""
echo "=========================================="
echo "SIMULATION COMPLETE!"
echo "=========================================="
echo "Job completed at: $(date)"
echo ""

echo "Output files:"
echo "  - md.xtc       : Trajectory (100 ns, includes waters)"
echo "  - md.edr       : Energy file"
echo "  - md.log       : Detailed log"
echo "  - md.gro       : Final structure"
echo "  - md.cpt       : Checkpoint (for restart)"
echo "  - topol.top    : System topology"
echo ""

echo "System configuration:"
echo "  - Container    : NGC GROMACS ${GROMACS_TAG}"
echo "  - Forcefield   : AMBER99SB-ILDN"
echo "  - Water model  : TIP3P"
echo "  - Temperature  : 300 K (26.85°C)"
echo "  - Pressure     : 1 bar"
echo "  - Salt conc.   : 0.15 M NaCl (1:1 ratio)"
echo "  - Box shape    : Cubic (1.2 nm buffer)"
echo "  - Duration     : 100 ns"
echo "  - GPUs         : ${NUM_GPUS}x V100"
echo "  - Direct GPU   : Enabled"
echo "  - Started from : Step ${START_STEP}"
echo ""

echo "Performance summary:"
if [ -f md.log ]; then
    grep "Performance" md.log | tail -5
    echo ""
    grep -E "ns/day|hour/ns" md.log | tail -2
fi

echo ""
echo "GPU utilization during run:"
if command -v nvidia-smi &> /dev/null; then
    nvidia-smi --query-gpu=utilization.gpu,memory.used,memory.total --format=csv
fi

echo ""
echo "=========================================="
echo "Next steps:"
echo "  1. Check simulation quality: bash quick_analysis.sh"
echo "  2. Visualize: Download md.xtc and view in PyMOL/VMD"
echo "  3. Analyze binding: See INTERPRETATION_GUIDE.txt"
echo "=========================================="
