#!/bin/bash
# Script to prepare all recligand structures for MD simulation
# This runs the preprocessing workflow for each structure

cd /home/pi/kuhfeldr/md

# Create log file with timestamp
LOG_FILE="preparation_log_$(date +%Y%m%d_%H%M%S).txt"

# Function to print and log
log_and_echo() {
    echo "$1" | tee -a "$LOG_FILE"
}

log_and_echo "=========================================="
log_and_echo "Preparing all structures for MD simulation"
log_and_echo "=========================================="
log_and_echo "Started: $(date)"
log_and_echo ""

# List of structures to prepare
STRUCTURES=(
    "LLFKDSAIGF_3fxi_LLFKDSAIGF_3fxi_input_docked_0003.pdb"
    "GLAPYKLRPVAA_3fxi_GLAPYKLRPVAA_3fxi_input_docked_0005.pdb"
    "RPKLPLRYP_3fxi_RPKLPLRYP_3fxi_input_docked_0008.pdb"
)

# Submit preparation jobs
JOB_IDS=()
for STRUCTURE in "${STRUCTURES[@]}"; do
    if [ ! -f "$STRUCTURE" ]; then
        log_and_echo "WARNING: $STRUCTURE not found, skipping..."
        continue
    fi
    
    STRUCTURE_NAME="${STRUCTURE%.*}"
    
    # Check if already prepared (md.tpr exists)
    if [ -f "${STRUCTURE_NAME}/md.tpr" ]; then
        log_and_echo "SKIP: ${STRUCTURE_NAME} already prepared (md.tpr exists)"
        log_and_echo ""
        continue
    fi
    
    log_and_echo "Submitting preparation job for ${STRUCTURE}..."
    JOB_ID=$(sbatch prepare_structure_wrapper.sh "$STRUCTURE" | awk '{print $4}')
    JOB_IDS+=($JOB_ID)
    log_and_echo "  Job ID: $JOB_ID"
    log_and_echo ""
done

log_and_echo "=========================================="
log_and_echo "All preparation jobs submitted!"
log_and_echo "Job IDs: ${JOB_IDS[*]}"
log_and_echo ""
log_and_echo "Check status with: squeue -u $USER"
log_and_echo "Monitor progress: tail -f <structure_name>/<structure_name>_prep_<job_id>.out"
log_and_echo ""
log_and_echo "After preparation completes, run MD simulations with:"
log_and_echo "  ./submit_all_recligand.sh"
log_and_echo ""
log_and_echo "Log saved to: $LOG_FILE"
log_and_echo "=========================================="

