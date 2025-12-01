#!/bin/bash
# Script to submit 100ns MD simulations for top 3 peptide-receptor complexes
# Based on 10ns analysis results

cd /home/pi/kuhfeldr/md

# Create log file with timestamp
LOG_FILE="job_submissions_$(date +%Y%m%d_%H%M%S).txt"

# Function to print and log
log_and_echo() {
    echo "$1" | tee -a "$LOG_FILE"
}

log_and_echo "=========================================="
log_and_echo "Submitting 100ns MD simulations"
log_and_echo "Based on 10ns binding analysis results"
log_and_echo "=========================================="
log_and_echo ""

# Structure 1: LLFKDSAIGF_0003 - Best stability (lowest RMSD, closest contact)
log_and_echo "1. LLFKDSAIGF_0003 - Best binding stability"
log_and_echo "   Submitting LLFKDSAIGF_3fxi_LLFKDSAIGF_3fxi_input_docked_0003..."
JOB1=$(sbatch run_md_1080ti.sh LLFKDSAIGF_3fxi_LLFKDSAIGF_3fxi_input_docked_0003.pdb | awk '{print $4}')
log_and_echo "   Job ID: $JOB1"
log_and_echo ""

# Structure 2: GLAPYKLRPVAA_0005 - Strongest H-bond network
log_and_echo "2. GLAPYKLRPVAA_0005 - Strongest H-bond network (avg 10.3)"
log_and_echo "   Submitting GLAPYKLRPVAA_3fxi_GLAPYKLRPVAA_3fxi_input_docked_0005..."
JOB2=$(sbatch run_md_1080ti.sh GLAPYKLRPVAA_3fxi_GLAPYKLRPVAA_3fxi_input_docked_0005.pdb | awk '{print $4}')
log_and_echo "   Job ID: $JOB2"
log_and_echo ""

# Structure 3: RPKLPLRYP_0008 - Best docking score (monitoring stability)
log_and_echo "3. RPKLPLRYP_0008 - Best docking score (-772)"
log_and_echo "   Submitting RPKLPLRYP_3fxi_RPKLPLRYP_3fxi_input_docked_0008..."
JOB3=$(sbatch run_md_1080ti.sh RPKLPLRYP_3fxi_RPKLPLRYP_3fxi_input_docked_0008.pdb | awk '{print $4}')
log_and_echo "   Job ID: $JOB3"
log_and_echo ""

log_and_echo "=========================================="
log_and_echo "All 3 jobs submitted successfully!"
log_and_echo "Job IDs: $JOB1, $JOB2, $JOB3"
log_and_echo ""
log_and_echo "Expected runtime: ~11 days per simulation"
log_and_echo "Check status with: squeue -u $USER"
log_and_echo "Log saved to: $LOG_FILE"
log_and_echo "=========================================="
