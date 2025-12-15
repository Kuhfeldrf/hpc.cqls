#!/bin/bash
#SBATCH --job-name=plot_md
#SBATCH --output=logs/plot_%j.out
#SBATCH --error=logs/plot_%j.err
#SBATCH --partition=short
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

# Generate plots for all simulations
# Usage: sbatch --export=STRUCTURES_LIST=STRUCT1,STRUCT2 submit_plot_all.sh
#        OR: sbatch submit_plot_all.sh STRUCT1 STRUCT2 (command line args)

MD_DIR="/home/kuhfeldr/hpc.cqls/md"
cd "${MD_DIR}"

# Load Python environment if needed
module load intel-python/24.0.0

echo "============================================================"
echo "Generating plots for all simulations"
echo "============================================================"
echo ""

# Get structures from environment variable (set by run_full_pipeline.sh) or command line arguments
if [ -n "${STRUCTURES_LIST}" ]; then
    # Parse comma-separated list from environment variable
    IFS=',' read -ra STRUCTURES <<< "${STRUCTURES_LIST}"
else
    # Get from command line arguments
    STRUCTURES=("$@")
fi

if [ ${#STRUCTURES[@]} -eq 0 ]; then
    echo "ERROR: No structures provided!"
    echo "Usage: sbatch --export=STRUCTURES_LIST=STRUCT1,STRUCT2 submit_plot_all.sh"
    echo "   OR: sbatch submit_plot_all.sh STRUCT1 STRUCT2"
    exit 1
fi

for STRUCT in "${STRUCTURES[@]}"; do
    [ -z "$STRUCT" ] && continue
    STRUCT_DIR="${MD_DIR}/${STRUCT}"
    
    if [ -d "${STRUCT_DIR}/analysis_output" ]; then
        echo "Generating plots for: ${STRUCT}"
        
        python3 "${MD_DIR}/plot_simulations.py" "${STRUCT_DIR}" --peptide "${STRUCT%%_3fxi*}"
        echo ""
    else
        echo "WARNING: No analysis output found for ${STRUCT}"
    fi
done

echo "============================================================"
echo "Plot generation complete!"
echo "============================================================"
