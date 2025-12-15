#!/bin/bash
# =============================================================================
# FULL MD PIPELINE V25 WITH JOB DEPENDENCIES
# =============================================================================
# This script orchestrates the complete MD workflow:
#   1. Prepare all structures in parallel (array job)
#   2. Run production MD for each structure (waits for prep)
#   3. Analyze all trajectories (waits for all MD jobs)
#   4. Generate plots (waits for analysis)
#
# Usage: ./run_full_pipeline.sh STRUCT1 STRUCT2 [STRUCT3 ...]  # Specify structures explicitly
#        ./run_full_pipeline.sh --dry-run STRUCT1 STRUCT2      # Show what would be submitted
#        ./run_full_pipeline.sh --skip-prep STRUCT1 STRUCT2   # Skip prep if already done
#
# =============================================================================

set -e

# Configuration
MD_DIR="/home/kuhfeldr/hpc.cqls/md"
cd "${MD_DIR}"

# Parse arguments
DRY_RUN=false
SKIP_PREP=false
SKIP_MD=false
STRUCTURES=()

for arg in "$@"; do
    case $arg in
        --dry-run)
            DRY_RUN=true
            echo "DRY RUN MODE - No jobs will be submitted"
            ;;
        --skip-prep)
            SKIP_PREP=true
            echo "Skipping preparation step"
            ;;
        --skip-md)
            SKIP_MD=true
            echo "Skipping MD step"
            ;;
        --*)
            echo "Unknown option: $arg"
            exit 1
            ;;
        *)
            # Treat as structure name (strip .pdb if provided)
            STRUCTURES+=("${arg%.pdb}")
            ;;
    esac
done

# Verify we have structures to process
if [ ${#STRUCTURES[@]} -eq 0 ]; then
    echo "ERROR: No structures provided!"
    echo ""
    echo "Usage: ./run_full_pipeline.sh STRUCT1 STRUCT2 [STRUCT3 ...]"
    echo ""
    echo "Examples:"
    echo "  ./run_full_pipeline.sh AVAVVKKGGSFQL_3fxi GLAPYKLRPVAA_3fxi"
    echo "  ./run_full_pipeline.sh --dry-run AVAVVKKGGSFQL_3fxi GLAPYKLRPVAA_3fxi"
    echo "  ./run_full_pipeline.sh --skip-prep AVAVVKKGGSFQL_3fxi GLAPYKLRPVAA_3fxi"
    exit 1
fi

# Normalize structure names (handle _3fxi suffix)
for i in "${!STRUCTURES[@]}"; do
    struct_name="${STRUCTURES[$i]}"
    # If structure name doesn't have _3fxi suffix, check if PDB file with _3fxi exists
    if [[ ! "$struct_name" == *_3fxi ]]; then
        if [ -f "${MD_DIR}/${struct_name}_3fxi.pdb" ]; then
            struct_name="${struct_name}_3fxi"
        fi
    fi
    STRUCTURES[$i]="$struct_name"
done

# Create comma-separated structure list for environment variable
STRUCTURES_LIST=$(IFS=','; echo "${STRUCTURES[*]}")
export STRUCTURES_LIST
echo "Structure list: ${STRUCTURES_LIST}"

NUM_STRUCTURES=${#STRUCTURES[@]}

echo ""
echo "============================================================"
echo "         FULL MD SIMULATION PIPELINE V25"
echo "============================================================"
echo "Started: $(date)"
echo "Working directory: ${MD_DIR}"
echo "Structures to process: ${NUM_STRUCTURES}"
echo ""
for i in "${!STRUCTURES[@]}"; do
    echo "  [$i] ${STRUCTURES[$i]}"
done
echo "============================================================"
echo ""

# Create logs directory
mkdir -p "${MD_DIR}/logs"

# Helper function to submit jobs
submit_job() {
    local cmd="$1"
    if [ "$DRY_RUN" = true ]; then
        echo "[DRY RUN] Would run: $cmd"
        echo "DRY_RUN_JOB_$$"
    else
        eval "$cmd"
    fi
}

# =============================================================================
# STEP 1: PREPARATION (Parallel Array Job)
# =============================================================================
echo "============================================================"
echo "STEP 1: Structure Preparation"
echo "============================================================"

PREP_JOB_ID=""

if [ "$SKIP_PREP" = true ]; then
    echo "Skipping preparation (--skip-prep flag set)"
    echo ""
    # Check if structures are actually prepared
    ALL_PREPARED=true
    for STRUCT in "${STRUCTURES[@]}"; do
        if [ ! -f "${MD_DIR}/${STRUCT}/md.tpr" ]; then
            echo "WARNING: ${STRUCT}/md.tpr not found!"
            ALL_PREPARED=false
        fi
    done
    if [ "$ALL_PREPARED" = false ]; then
        echo ""
        echo "ERROR: Some structures are not prepared. Remove --skip-prep flag."
        exit 1
    fi
else
    echo "Submitting parallel preparation job..."
    echo "  Script: prepare_structures_parallel.sh"
    echo "  Array: 0-$((NUM_STRUCTURES-1))"
    echo ""
    
    # Override the array size in the script dynamically and pass structures list
    # Use ALL to pass all environment variables (including STRUCTURES_LIST)
    PREP_JOB_ID=$(submit_job "sbatch --parsable --array=0-$((NUM_STRUCTURES-1)) --export=ALL prepare_structures_parallel.sh")
    echo "✓ Preparation job submitted: ${PREP_JOB_ID}"
    echo ""
fi

# =============================================================================
# STEP 2: PRODUCTION MD (One job per structure, depends on prep)
# =============================================================================
echo "============================================================"
echo "STEP 2: Production MD Simulations"
echo "============================================================"

MD_JOB_IDS=()

if [ "$SKIP_MD" = true ]; then
    echo "Skipping MD (--skip-md flag set)"
    echo ""
else
    for i in "${!STRUCTURES[@]}"; do
        STRUCT="${STRUCTURES[$i]}"
        STRUCT_DIR="${MD_DIR}/${STRUCT}"
        
        echo "Submitting MD job for: ${STRUCT}"
        
        # Build dependency string
        DEP_STRING=""
        if [ -n "$PREP_JOB_ID" ] && [ "$PREP_JOB_ID" != "DRY_RUN_JOB_$$" ]; then
            # Wait for corresponding prep array task
            DEP_STRING="--dependency=aftercorr:${PREP_JOB_ID}"
        fi
        
        # Submit MD job
        MD_CMD="sbatch --parsable ${DEP_STRING} --chdir=${STRUCT_DIR} ${MD_DIR}/submit_gromacs_gpu.sh"
        MD_JOB_ID=$(submit_job "$MD_CMD")
        MD_JOB_IDS+=("$MD_JOB_ID")
        
        echo "  ✓ Job ${MD_JOB_ID} (depends on prep task $i)"
    done
    echo ""
fi

# =============================================================================
# STEP 3: ANALYSIS (One job per structure, depends on its MD job)
# =============================================================================
echo "============================================================"
echo "STEP 3: Trajectory Analysis"
echo "============================================================"

# Create analysis submission script (runs analyze_md.sh via SLURM)
cat > "${MD_DIR}/submit_analyze_single.sh" << 'ANALYZE_SCRIPT'
#!/bin/bash
#SBATCH --job-name=analyze_md
#SBATCH --output=logs/analyze_%j.out
#SBATCH --error=logs/analyze_%j.err
#SBATCH --partition=short
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --mem=32G

# Analyze a single simulation directory
# Usage: sbatch submit_analyze_single.sh /path/to/simulation_dir

module load apptainer/1.4.1-gcc-13.4.0

SIM_DIR="$1"
if [ -z "$SIM_DIR" ]; then
    echo "ERROR: No simulation directory provided"
    exit 1
fi

cd /home/kuhfeldr/hpc.cqls/md
./analyze_md.sh "$SIM_DIR"
ANALYZE_SCRIPT

chmod +x "${MD_DIR}/submit_analyze_single.sh"

ANALYSIS_JOB_IDS=()

for i in "${!STRUCTURES[@]}"; do
    STRUCT="${STRUCTURES[$i]}"
    STRUCT_DIR="${MD_DIR}/${STRUCT}"
    
    echo "Submitting analysis job for: ${STRUCT}"
    
    # Build dependency string (wait for corresponding MD job)
    DEP_STRING=""
    if [ ${#MD_JOB_IDS[@]} -gt 0 ] && [ "${MD_JOB_IDS[$i]}" != "DRY_RUN_JOB_$$" ]; then
        DEP_STRING="--dependency=afterok:${MD_JOB_IDS[$i]}"
    fi
    
    # Submit analysis job
    ANALYSIS_CMD="sbatch --parsable ${DEP_STRING} ${MD_DIR}/submit_analyze_single.sh ${STRUCT_DIR}"
    ANALYSIS_JOB_ID=$(submit_job "$ANALYSIS_CMD")
    ANALYSIS_JOB_IDS+=("$ANALYSIS_JOB_ID")
    
    echo "  ✓ Job ${ANALYSIS_JOB_ID} (depends on MD job ${MD_JOB_IDS[$i]:-N/A})"
done
echo ""

# =============================================================================
# STEP 4: PLOTTING (Single job, depends on ALL analysis jobs)
# =============================================================================
echo "============================================================"
echo "STEP 4: Generate Plots"
echo "============================================================"

# Use existing submit_plot_all.sh script
if [ ! -f "${MD_DIR}/submit_plot_all.sh" ]; then
    echo "ERROR: ${MD_DIR}/submit_plot_all.sh not found!"
    exit 1
fi

# Build dependency string for all analysis jobs
PLOT_DEP=""
if [ ${#ANALYSIS_JOB_IDS[@]} -gt 0 ]; then
    # Filter out dry run placeholders
    REAL_ANALYSIS_IDS=()
    for id in "${ANALYSIS_JOB_IDS[@]}"; do
        if [[ "$id" != "DRY_RUN_JOB_"* ]]; then
            REAL_ANALYSIS_IDS+=("$id")
        fi
    done
    
    if [ ${#REAL_ANALYSIS_IDS[@]} -gt 0 ]; then
        PLOT_DEP="--dependency=afterok:$(IFS=:; echo "${REAL_ANALYSIS_IDS[*]}")"
    fi
fi

echo "Submitting plotting job..."
# Pass structures as environment variable to plotting script (use ALL for consistency)
PLOT_CMD="sbatch --parsable ${PLOT_DEP} --export=ALL ${MD_DIR}/submit_plot_all.sh"
PLOT_JOB_ID=$(submit_job "$PLOT_CMD")
echo "✓ Plot job submitted: ${PLOT_JOB_ID}"
echo ""

# =============================================================================
# SUMMARY
# =============================================================================
echo "============================================================"
echo "                    PIPELINE SUMMARY"
echo "============================================================"
echo ""
echo "Jobs submitted:"
echo ""

if [ -n "$PREP_JOB_ID" ]; then
    echo "  PREPARATION (array job):"
    echo "    Job ID: ${PREP_JOB_ID}"
    echo "    Tasks: 0-$((NUM_STRUCTURES-1))"
    echo ""
fi

if [ ${#MD_JOB_IDS[@]} -gt 0 ]; then
    echo "  PRODUCTION MD:"
    for i in "${!STRUCTURES[@]}"; do
        echo "    ${STRUCTURES[$i]}: ${MD_JOB_IDS[$i]:-skipped}"
    done
    echo ""
fi

if [ ${#ANALYSIS_JOB_IDS[@]} -gt 0 ]; then
    echo "  ANALYSIS:"
    for i in "${!STRUCTURES[@]}"; do
        echo "    ${STRUCTURES[$i]}: ${ANALYSIS_JOB_IDS[$i]:-skipped}"
    done
    echo ""
fi

echo "  PLOTTING:"
echo "    Job ID: ${PLOT_JOB_ID}"
echo ""

echo "============================================================"
echo "                  DEPENDENCY CHAIN"
echo "============================================================"
echo ""
echo "  PREP (parallel) → MD (parallel per structure) → ANALYSIS → PLOT"
echo ""
echo "  Each MD job waits for its corresponding prep task"
echo "  Each analysis job waits for its MD job"
echo "  Plot job waits for ALL analysis jobs"
echo ""
echo "============================================================"
echo ""
echo "Monitor progress:"
echo "  squeue -u $USER"
echo ""
echo "View job details:"
echo "  scontrol show job <JOB_ID>"
echo ""
echo "Cancel all jobs:"
if [ -n "$PREP_JOB_ID" ]; then
    echo "  scancel ${PREP_JOB_ID} ${MD_JOB_IDS[*]} ${ANALYSIS_JOB_IDS[*]} ${PLOT_JOB_ID}"
fi
echo ""
echo "Logs will be in: ${MD_DIR}/logs/"
echo ""
echo "============================================================"
echo "Pipeline submitted at: $(date)"
echo "============================================================"

