#!/bin/bash
# Script to submit all three recligand MD simulations in parallel

cd /home/pi/kuhfeldr/md

# Create log file with timestamp
LOG_FILE="job_submissions_$(date +%Y%m%d_%H%M%S).txt"

# Function to print and log
log_and_echo() {
    echo "$1" | tee -a "$LOG_FILE"
}

log_and_echo "Submitting 3 recligand MD simulations..."
log_and_echo ""

# Submit each job
log_and_echo "Submitting GLAPYKLRPVAA_recligand_lr1.pdb..."
JOB1=$(sbatch run_md_1080ti.sh GLAPYKLRPVAA_recligand_lr1.pdb | awk '{print $4}')
log_and_echo "  Job ID: $JOB1"
log_and_echo ""

log_and_echo "Submitting LLFKDSAIGF_recligand_lr6.pdb..."
JOB2=$(sbatch run_md_1080ti.sh LLFKDSAIGF_recligand_lr6.pdb | awk '{print $4}')
log_and_echo "  Job ID: $JOB2"
log_and_echo ""

log_and_echo "Submitting RPKLPLRYP_recligand_lr2.pdb..."
JOB3=$(sbatch run_md_1080ti.sh RPKLPLRYP_recligand_lr2.pdb | awk '{print $4}')
log_and_echo "  Job ID: $JOB3"
log_and_echo ""

log_and_echo "=========================================="
log_and_echo "All jobs submitted successfully!"
log_and_echo "Job IDs: $JOB1, $JOB2, $JOB3"
log_and_echo ""
log_and_echo "Check status with: squeue -u $USER"
log_and_echo "Check specific job with: squeue -j <JOB_ID>"
log_and_echo "Log saved to: $LOG_FILE"
log_and_echo "=========================================="

