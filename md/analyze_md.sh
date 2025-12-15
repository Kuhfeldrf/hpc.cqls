#!/bin/bash

# GROMACS MD Analysis Script with Automated PBC Checking
# ==========================================
# Performs analysis of GROMACS MD trajectories:
# - Automatic PBC artifact detection and clustering
# - Backbone RMSD (receptor stability)
# - Peptide RMSD (binding stability & conformational flexibility)
# - RMSF (flexibility per residue)
# - Radius of Gyration (compactness)
# - Hydrogen Bonds (peptide-receptor)
#
# RMSD Analysis Strategy for Binding Studies:
# 1. Receptor RMSD (fit receptor, calc receptor) - Is the protein stable?
# 2. Peptide Internal RMSD (fit peptide, calc peptide) - Does peptide change conformation?
# 3. Peptide Binding RMSD (fit receptor, calc peptide) - Does peptide stay in binding site?
#
# For extended analyses (DSSP, MinDist), use analyze_md_extended.sh
# Submit via: sbatch submit_analyze_md.sh <directory>

set -o pipefail

# =============================================================================
# CONFIGURATION
# =============================================================================
ORIGINAL_DIR=$(pwd)
MAIN_DIR="analysis_output"
CONTAINER_PATH="/home/kuhfeldr/gromacs_v2025.3.sif"

# PBC checking thresholds
PBC_RMSD_THRESHOLD=0.4  # nm - spike in RMSD that indicates PBC artifacts
PBC_CHECK_ENABLED=true   # Set to false to disable automatic PBC checking

# Detect if running interactively (not via SLURM)
INTERACTIVE=false
if [ -z "$SLURM_JOB_ID" ] && [ -t 1 ]; then
    INTERACTIVE=true
fi

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

log() {
    local msg="$1"
    echo "$msg"
    [ -n "$LOG_FILE" ] && echo "$msg" >> "$LOG_FILE"
}

log_error() {
    local msg="ERROR: $1"
    local err_file="$2"
    echo "$msg" >&2
    [ -n "$LOG_FILE" ] && echo "$msg" >> "$LOG_FILE"
    
    if [ "$INTERACTIVE" = true ] && [ -n "$err_file" ] && [ -f "$err_file" ]; then
        echo "--- Error details ---" >&2
        cat "$err_file" >&2
        echo "---------------------" >&2
    fi
}

file_ok() {
    [ -f "$1" ] && [ -s "$1" ]
}

get_index_group() {
    local group_name="$1"
    local index_file="$2"
    if [ -f "$index_file" ]; then
        awk -v target="$group_name" '
            /^\[.*\]/ {
                gsub(/[\[\] ]/, "", $0)
                if ($0 == target) { print group_num; exit }
                group_num++
            }
        ' "$index_file"
    fi
}

xvg_stats() {
    local file="$1"
    local col="${2:-2}"
    grep -v "^[#@]" "$file" | awk -v c="$col" '
        BEGIN { min=999999; max=-999999; sum=0; n=0 }
        { 
            if($c < min) min=$c
            if($c > max) max=$c
            sum+=$c
            n++
            if(n==1) first=$c
            last=$c
        }
        END { 
            if(n>0) printf "%.3f %.3f %.3f %.3f %.3f", first, last, sum/n, min, max
            else print "N/A N/A N/A N/A N/A"
        }'
}

get_peptide_len() {
    local dirname="$1"
    local peptide_seq=$(echo "$dirname" | sed -E 's/^([A-Z]+)_.*/\1/')
    if [[ "$peptide_seq" =~ ^[A-Z]+$ ]]; then
        echo "${#peptide_seq}"
    else
        echo 0
    fi
}

create_index_file() {
    local tpr_file="$1"
    local index_file="$2"
    local existing_file="${3:-}"  # Optional: existing index file to read from
    
    echo ""
    if [ -n "$existing_file" ] && [ -f "$existing_file" ]; then
        echo "--- Adding Receptor and Peptide groups to existing index file ---"
    else
        echo "--- Creating new index file with Receptor and Peptide groups ---"
    fi
    
    local MAKE_NDX_ERR=".make_ndx_error.log"
    
    # Use splitch to split Protein group by chain - creates Protein_chain1, chain2, etc.
    # Based on topology: chains 1-4 are receptor (A,B,C,D), chain 5 is peptide (P)
    # Then combine the receptor chains and name them
    local make_ndx_input="splitch 1
17 | 18 | 19 | 20
name 21 Receptor
22
name 23 Peptide
q
"
    
    # Build make_ndx command - only use -n if existing file is provided and exists
    local make_ndx_cmd="make_ndx -f $tpr_file"
    if [ -n "$existing_file" ] && [ -f "$existing_file" ]; then
        make_ndx_cmd="$make_ndx_cmd -n $existing_file"
        log "Reading existing index file: $existing_file"
    else
        log "Creating new index file"
    fi
    make_ndx_cmd="$make_ndx_cmd -o $index_file"
    
    log "Running gmx make_ndx to add Receptor and Peptide groups..."
    if run_gmx_with_input "$make_ndx_input" $make_ndx_cmd 2>"$MAKE_NDX_ERR"; then
        if [ -f "$index_file" ] && [ -s "$index_file" ]; then
            # Verify groups exist
            local peptide_grp=$(get_index_group "Peptide" "$index_file")
            local receptor_grp=$(get_index_group "Receptor" "$index_file")
            
            if [ -n "$peptide_grp" ] && [ -n "$receptor_grp" ]; then
                rm -f "$MAKE_NDX_ERR" 2>/dev/null
                log "Successfully added groups to index file - Peptide: $peptide_grp | Receptor: $receptor_grp"
                return 0
            else
                log "Warning: Groups not found after creation"
                return 1
            fi
        else
            log_error "Failed to update index file" "$MAKE_NDX_ERR"
            return 1
        fi
    else
        log_error "Failed to add groups to index file" "$MAKE_NDX_ERR"
        return 1
    fi
}

# =============================================================================
# GROMACS CONTAINER SETUP
# =============================================================================
setup_gromacs() {
    USE_CONTAINER=false
    
    if command -v gmx &> /dev/null; then
        echo "Using direct GROMACS installation"
        return 0
    fi
    
    if [ ! -f "$CONTAINER_PATH" ]; then
        echo "WARNING: Neither gmx command nor container found at $CONTAINER_PATH"
        return 1
    fi
    
    echo "GROMACS container found at: $CONTAINER_PATH"
    
    # Try to load apptainer module
    if command -v module &> /dev/null; then
        echo "Attempting to load apptainer module..."
        module load apptainer/1.4.1-gcc-13.4.0 2>&1 || module load apptainer 2>&1 || true
    fi
    
    # Set up container command
    if command -v apptainer &> /dev/null; then
        USE_CONTAINER=true
        CONTAINER_CMD="apptainer"
        echo "Using GROMACS container via apptainer"
    elif command -v singularity &> /dev/null; then
        USE_CONTAINER=true
        CONTAINER_CMD="singularity"
        echo "Using GROMACS container via singularity"
    else
        echo "ERROR: Container found but neither apptainer nor singularity is available"
        echo "Submit via sbatch: sbatch submit_analyze_md.sh <directory>"
        exit 1
    fi
}

run_gmx() {
    if [ "$USE_CONTAINER" = true ]; then
        local CURRENT_DIR=$(pwd)
        
        # Use 'exec' with full path to gmx binary instead of 'run'
        # This bypasses the container's runscript which checks for GPUs
        # The CPU-only gmx binary is at /usr/local/gromacs/avx2_256/bin/gmx
        ${CONTAINER_CMD} exec -B "${CURRENT_DIR}:/host_pwd" --pwd /host_pwd "${CONTAINER_PATH}" \
            /usr/local/gromacs/avx2_256/bin/gmx "$@"
    else
        gmx "$@"
    fi
}

# Run GMX command with input - uses temporary input file
# This works around SLURM batch job stdin issues with containers
run_gmx_with_input() {
    local input="$1"
    shift
    
    # Write input to temporary file (works for both container and direct)
    local INPUT_FILE=".gmx_input_$$.txt"
    printf "%s\n" "$input" > "$INPUT_FILE"
    
    if [ "$USE_CONTAINER" = true ]; then
        local CURRENT_DIR=$(pwd)
        
        # Run with input redirection inside container
        # Since we're using --pwd /host_pwd, we can use relative path
        ${CONTAINER_CMD} exec -B "${CURRENT_DIR}:/host_pwd" --pwd /host_pwd "${CONTAINER_PATH}" \
            /usr/local/gromacs/avx2_256/bin/gmx "$@" < "$INPUT_FILE"
        local rc=$?
    else
        # Direct execution with input file
        gmx "$@" < "$INPUT_FILE"
        local rc=$?
    fi
    
    rm -f "$INPUT_FILE"
    return $rc
}

# Run GMX command and pipe output to another GMX command
# This handles the complexity of piping between container commands
run_gmx_pipe() {
    local input1="$1"
    shift
    local cmd1=("$@")
    
    if [ "$USE_CONTAINER" = true ]; then
        local CURRENT_DIR=$(pwd)
        
        # Create a wrapper script that handles the entire pipeline
        cat > .gmx_pipe_wrapper.sh << 'PIPE_EOF'
#!/bin/bash
# This script runs inside the container and handles the full pipeline
PIPE_EOF
        
        # Add the actual pipeline command
        echo "echo \"$input1\" | /usr/local/gromacs/avx2_256/bin/gmx ${cmd1[*]}" >> .gmx_pipe_wrapper.sh
        chmod +x .gmx_pipe_wrapper.sh
        
        # Execute inside container
        ${CONTAINER_CMD} exec -B "${CURRENT_DIR}:/host_pwd" --pwd /host_pwd "${CONTAINER_PATH}" \
            /host_pwd/.gmx_pipe_wrapper.sh
        local rc=$?
        rm -f .gmx_pipe_wrapper.sh
        return $rc
    else
        # Direct execution
        echo "$input1" | gmx "${cmd1[@]}"
    fi
}

# =============================================================================
# PBC ARTIFACT DETECTION
# =============================================================================

check_pbc_artifacts() {
    local tpr_file="$1"
    local xtc_file="$2"
    local threshold="${3:-$PBC_RMSD_THRESHOLD}"
    
    echo ""
    echo "--- Checking for PBC artifacts ---"
    
    # Create temporary check file
    local CHECK_XVG=".pbc_check_rmsd.xvg"
    local CHECK_ERR=".pbc_check_error.log"
    
    # Run quick RMSD to detect jumps
    echo "Running preliminary RMSD check..."
    if run_gmx_with_input "4
4" rms -s "$tpr_file" -f "$xtc_file" -o "$CHECK_XVG" 2>"$CHECK_ERR"; then
        rm -f "$CHECK_ERR"
    else
        log "Warning: PBC check RMSD calculation failed, assuming clustering needed"
        rm -f "$CHECK_XVG" "$CHECK_ERR"
        return 0  # Return 0 (needs clustering) if check fails
    fi
    
    if [ ! -f "$CHECK_XVG" ] || [ ! -s "$CHECK_XVG" ]; then
        log "Warning: PBC check produced no output, assuming clustering needed"
        return 0
    fi
    
    # Analyze RMSD for spikes using awk
    local result=$(awk -v threshold="$threshold" '
        BEGIN { 
            max_jump = 0
            n_spikes = 0
            prev_rmsd = -1
        }
        /^[^#@]/ && NF >= 2 {
            rmsd = $2
            if (prev_rmsd >= 0) {
                jump = (rmsd > prev_rmsd) ? (rmsd - prev_rmsd) : (prev_rmsd - rmsd)
                if (jump > max_jump) max_jump = jump
                if (jump > threshold) n_spikes++
            }
            prev_rmsd = rmsd
        }
        END {
            needs_clustering = (max_jump > threshold) ? 1 : 0
            printf "%d %.4f %d", needs_clustering, max_jump, n_spikes
        }
    ' "$CHECK_XVG")
    
    read -r needs_clustering max_jump n_spikes <<< "$result"
    
    # Clean up check files
    rm -f "$CHECK_XVG" "$CHECK_ERR"
    
    # Report results
    if [ "$needs_clustering" -eq 1 ]; then
        log "PBC artifacts detected:"
        log "  - Max RMSD jump: ${max_jump} nm (threshold: ${threshold} nm)"
        log "  - Number of spikes: ${n_spikes}"
        log "  → Applying PBC clustering"
        return 0  # Needs clustering
    else
        log "No PBC artifacts detected (max jump: ${max_jump} nm)"
        log "  → Using original trajectory"
        return 1  # No clustering needed
    fi
}

apply_pbc_clustering() {
    local tpr_file="$1"
    local input_xtc="$2"
    local output_xtc="$3"
    local index_file="$4"
    
    echo ""
    echo "--- Applying PBC clustering ---"
    
    local CLUSTER_ERR=".pbc_cluster_error.log"
    
    # Determine clustering and centering groups
    # For multi-chain systems (TLR4-MD2-peptide), we want to cluster on the entire protein complex
    local cluster_group="1"  # Protein (includes all chains)
    local center_group="1"   # Protein (includes all chains)
    local output_group="0"   # System (everything including water/ions)
    local index_flag=""
    
    if [ -n "$index_file" ] && [ -f "$index_file" ]; then
        index_flag="-n $index_file"
        
        # Check for Complex group (manually created group with all interacting molecules)
        local complex_group=$(get_index_group "Complex" "$index_file")
        if [ -n "$complex_group" ]; then
            cluster_group="$complex_group"
            center_group="$complex_group"
            log "Using 'Complex' group ($complex_group) for PBC clustering and centering"
        else
            # Use Protein group which should include all protein chains (receptor + peptide)
            log "Using 'Protein' group (1) for PBC clustering and centering"
            log "  This includes all protein chains: receptor, MD2, and peptide"
        fi
    else
        log "Using default 'Protein' group (1) for PBC clustering and centering"
    fi
    
    # Apply clustering with THREE selections:
    # 1. Group for clustering (which molecules to keep together)
    # 2. Group for centering (what to center in box)
    # 3. Group for output (what to write - always System to keep everything)
    
    # Create input file with three group selections
    local INPUT_FILE=".pbc_input.txt"
    printf "%s\n%s\n%s\n" "$cluster_group" "$center_group" "$output_group" > "$INPUT_FILE"
    
    if [ "$USE_CONTAINER" = true ]; then
        local CURRENT_DIR=$(pwd)
        # Since we're using --pwd /host_pwd, we can use relative path
        # Suppress verbose clustering output by redirecting both stdout and stderr
        ${CONTAINER_CMD} exec -B "${CURRENT_DIR}:/host_pwd" --pwd /host_pwd "${CONTAINER_PATH}" \
            /usr/local/gromacs/avx2_256/bin/gmx trjconv -s "$tpr_file" -f "$input_xtc" \
            -o "$output_xtc" $index_flag -pbc cluster -center < "$INPUT_FILE" >"$CLUSTER_ERR" 2>&1
        local rc=$?
    else
        # Suppress verbose clustering output by redirecting both stdout and stderr
        gmx trjconv -s "$tpr_file" -f "$input_xtc" -o "$output_xtc" $index_flag \
            -pbc cluster -center < "$INPUT_FILE" >"$CLUSTER_ERR" 2>&1
        local rc=$?
    fi
    
    rm -f "$INPUT_FILE"
    
    if [ $rc -eq 0 ] && file_ok "$output_xtc"; then
        rm -f "$CLUSTER_ERR"
        log "PBC clustering successful"
        log "  Clustering group: $cluster_group"
        log "  Centering group: $center_group"
        log "  Output: $output_xtc"
        return 0
    else
        log_error "PBC clustering failed" "$CLUSTER_ERR"
        return 1
    fi
}

# =============================================================================
# CREATE PEPTIDE-FITTED REFERENCE STRUCTURE
# =============================================================================
# For peptide internal RMSD, we need a reference structure where the peptide
# is used as the fitting group. This extracts the first frame with peptide-only fitting.

create_peptide_reference() {
    local tpr_file="$1"
    local xtc_file="$2"
    local index_file="$3"
    local peptide_group="$4"
    local output_gro="$5"
    
    local REF_ERR=".peptide_ref_error.log"
    
    log "Creating peptide reference structure for internal RMSD..."
    
    # Extract first frame, outputting only peptide atoms
    # This will be used as the reference for peptide internal RMSD
    if run_gmx_with_input "$peptide_group" trjconv -s "$tpr_file" -f "$xtc_file" \
        -n "$index_file" -o "$output_gro" -dump 0 2>"$REF_ERR"; then
        if file_ok "$output_gro"; then
            rm -f "$REF_ERR"
            log "Peptide reference structure created: $output_gro"
            return 0
        fi
    fi
    
    log_error "Failed to create peptide reference structure" "$REF_ERR"
    return 1
}

# =============================================================================
# MAIN ANALYSIS FUNCTION
# =============================================================================
analyze_directory() {
    local dir="$1"
    dir="${dir%/}"
    
    if [ ! -d "$dir" ]; then
        echo "Error: Directory '$dir' does not exist"
        return 1
    fi
    
    dir=$(cd "$dir" 2>/dev/null && pwd) || { echo "Error: Cannot access directory '$1'"; return 1; }
    local name=$(basename "$dir")
    
    echo ""
    echo "=========================================="
    echo "Analyzing: $name"
    echo "=========================================="
    
    # Set output directories (absolute paths for shell operations)
    local BASE_OUTPUT_DIR="${dir}/${MAIN_DIR}"
    local OUTPUT_DIR="${BASE_OUTPUT_DIR}/output_files"
    local LOG_DIR="${BASE_OUTPUT_DIR}/processing_logs"
    local ANALYSIS_DIR="${BASE_OUTPUT_DIR}/analysis_data"
    
    local SIM_ANALYSIS_DIR="${ANALYSIS_DIR}/${name}"
    local SIM_LOG_DIR="${LOG_DIR}/${name}"
    mkdir -p "$SIM_ANALYSIS_DIR" "$SIM_LOG_DIR" "$OUTPUT_DIR"
    
    # Relative paths from simulation dir (for GROMACS container commands)
    local REL_ANALYSIS_DIR="${MAIN_DIR}/analysis_data/${name}"
    
    # Initialize log file
    LOG_FILE="${SIM_LOG_DIR}/analysis_$(date +%Y%m%d_%H%M%S).log"
    log "Analysis started: $(date)"
    log "Directory: $dir"
    
    # =========================================================================
    # CHANGE TO SIMULATION DIRECTORY
    # =========================================================================
    cd "$dir" || { echo "Error: Cannot cd to $dir"; return 1; }
    
    # =========================================================================
    # FILE VERIFICATION
    # =========================================================================
    if [ ! -f "md.tpr" ]; then
        log_error "md.tpr not found in $dir"
        cd "$ORIGINAL_DIR"
        return 1
    fi
    
    # Find trajectory file
    local XTC_FILE=""
    if [ -f "md_nopbc.xtc" ]; then
        XTC_FILE="md_nopbc.xtc"
    elif [ -f "md.xtc" ]; then
        XTC_FILE="md.xtc"
    elif [ -f "traj.xtc" ]; then
        XTC_FILE="traj.xtc"
    else
        log_error "No trajectory file found (tried md_nopbc.xtc, md.xtc, traj.xtc)"
        cd "$ORIGINAL_DIR"
        return 1
    fi
    log "Using trajectory: $XTC_FILE"
    
    # =========================================================================
    # INDEX FILE DETECTION AND CREATION
    # =========================================================================
    INDEX_FILE=""
    PEPTIDE_GROUP=""
    RECEPTOR_GROUP=""
    
    if [ -f "index.ndx" ]; then
        INDEX_FILE="index.ndx"
        log "Found custom index file: $INDEX_FILE"
        
        PEPTIDE_GROUP=$(get_index_group "Peptide" "$INDEX_FILE")
        RECEPTOR_GROUP=$(get_index_group "Receptor" "$INDEX_FILE")
        
        if [ -n "$PEPTIDE_GROUP" ] && [ -n "$RECEPTOR_GROUP" ]; then
            log "Custom groups found - Peptide: $PEPTIDE_GROUP | Receptor: $RECEPTOR_GROUP"
        else
            log "Warning: Custom index file present but missing Peptide/Receptor groups"
            log "Attempting to add required groups to existing index file..."
            if create_index_file "md.tpr" "index.ndx" "index.ndx"; then
                INDEX_FILE="index.ndx"
                PEPTIDE_GROUP=$(get_index_group "Peptide" "$INDEX_FILE")
                RECEPTOR_GROUP=$(get_index_group "Receptor" "$INDEX_FILE")
            fi
        fi
    else
        log "No custom index file found, attempting to create one..."
        if create_index_file "md.tpr" "index.ndx" ""; then
            INDEX_FILE="index.ndx"
            PEPTIDE_GROUP=$(get_index_group "Peptide" "$INDEX_FILE")
            RECEPTOR_GROUP=$(get_index_group "Receptor" "$INDEX_FILE")
            log "Index file created - Peptide: $PEPTIDE_GROUP | Receptor: $RECEPTOR_GROUP"
        else
            log "Failed to create index file, using default groups"
            log "Note: Peptide-specific analyses will be skipped"
        fi
    fi
    
    # =========================================================================
    # PBC ARTIFACT DETECTION AND CORRECTION
    # =========================================================================
    ANALYSIS_XTC="$XTC_FILE"  # Default to original trajectory
    CLUSTERED_XTC="${MAIN_DIR}/trajectory_clustered.xtc"
    
    if [ "$PBC_CHECK_ENABLED" = true ]; then
        if check_pbc_artifacts "md.tpr" "$XTC_FILE" "$PBC_RMSD_THRESHOLD"; then
            # PBC artifacts detected - check if clustered trajectory already exists
            if file_ok "$CLUSTERED_XTC"; then
                log "Found existing PBC-corrected trajectory: $CLUSTERED_XTC"
                ANALYSIS_XTC="$CLUSTERED_XTC"
                log "Using existing PBC-corrected trajectory for analysis"
            else
                # No existing clustered trajectory, create one
                log "No existing clustered trajectory found, creating new one..."
                if apply_pbc_clustering "md.tpr" "$XTC_FILE" "$CLUSTERED_XTC" "$INDEX_FILE"; then
                    ANALYSIS_XTC="$CLUSTERED_XTC"
                    log "Using newly created PBC-corrected trajectory for analysis"
                else
                    log "Warning: PBC clustering failed, using original trajectory"
                    ANALYSIS_XTC="$XTC_FILE"
                fi
            fi
        else
            # No artifacts detected, use original
            ANALYSIS_XTC="$XTC_FILE"
        fi
    else
        log "PBC checking disabled, using original trajectory"
    fi
    
    # =========================================================================
    # ATOM COUNT VERIFICATION (for custom index files)
    # =========================================================================
    get_traj_for_custom_index() {
        # If using custom index with different atom counts, we need the right trajectory
        if [ -n "$INDEX_FILE" ] && [ "$ANALYSIS_XTC" != "$XTC_FILE" ]; then
            # Check if clustered trajectory has correct atom count
            # Use awk instead of grep -P to avoid lookbehind issues
            # gmx check outputs: "Step X, Y atoms" - extract the number before "atoms"
            local tpr_atoms=$(run_gmx check -s md.tpr 2>&1 | awk '/atoms/ {for(i=1;i<=NF;i++) if($i=="atoms" && $(i-1)~/^[0-9]+$/) {print $(i-1); exit}}' | head -1)
            local xtc_atoms=$(run_gmx check -f "$ANALYSIS_XTC" 2>&1 | awk '/atoms/ {for(i=1;i<=NF;i++) if($i=="atoms" && $(i-1)~/^[0-9]+$/) {print $(i-1); exit}}' | head -1)
            
            if [ -n "$tpr_atoms" ] && [ -n "$xtc_atoms" ] && [ "$tpr_atoms" != "$xtc_atoms" ]; then
                log "Warning: Atom count mismatch (TPR: $tpr_atoms, XTC: $xtc_atoms)"
                log "Using original trajectory for custom index analysis"
                echo "$XTC_FILE"
            else
                echo "$ANALYSIS_XTC"
            fi
        else
            echo "$ANALYSIS_XTC"
        fi
    }
    
    # =========================================================================
    # ANALYSIS SECTION 1: RMSD
    # =========================================================================
    echo ""
    echo "--- RMSD Analysis ---"
    log ""
    log "RMSD Analysis Strategy:"
    log "  1. Receptor RMSD: Fit receptor backbone → Calculate receptor backbone RMSD"
    log "     → Measures: Receptor structural stability"
    log "  2. Peptide Internal RMSD: Fit peptide → Calculate peptide RMSD"
    log "     → Measures: Peptide conformational changes (independent of binding pose)"
    log "  3. Peptide Binding RMSD: Fit receptor backbone → Calculate peptide RMSD"
    log "     → Measures: Peptide movement relative to receptor (binding stability)"
    log ""
    
    # 1a. Receptor Backbone RMSD
    # Fit: Receptor backbone | Calculate: Receptor backbone
    # Purpose: Monitor receptor structural stability throughout simulation
    local RMSD_RECEPTOR_XVG="${REL_ANALYSIS_DIR}/rmsd_receptor_backbone.xvg"
    local RMSD_RECEPTOR_ERR="${SIM_LOG_DIR}/rmsd_receptor_error.log"
    
    if file_ok "$RMSD_RECEPTOR_XVG"; then
        echo "Using existing receptor RMSD file"
    elif [ -n "$INDEX_FILE" ] && file_ok "$INDEX_FILE" && [ -n "$RECEPTOR_GROUP" ]; then
        echo "Calculating receptor backbone RMSD (fit: receptor, calc: receptor)..."
        local TRAJ_FOR_ANALYSIS=$(get_traj_for_custom_index)
        # First group = fit group, Second group = calculation group
        # Both are receptor for receptor stability measurement
        if run_gmx_with_input "${RECEPTOR_GROUP}
${RECEPTOR_GROUP}" rms -s md.tpr -f "$TRAJ_FOR_ANALYSIS" -n "$INDEX_FILE" -o "$RMSD_RECEPTOR_XVG" 2>"$RMSD_RECEPTOR_ERR" && file_ok "$RMSD_RECEPTOR_XVG"; then
            rm -f "$RMSD_RECEPTOR_ERR" 2>/dev/null
        else
            log_error "Failed to calculate receptor RMSD" "$RMSD_RECEPTOR_ERR"
        fi
    else
        echo "Calculating protein backbone RMSD (Backbone group)..."
        if run_gmx_with_input "4
4" rms -s md.tpr -f "$ANALYSIS_XTC" -o "$RMSD_RECEPTOR_XVG" 2>"$RMSD_RECEPTOR_ERR" && file_ok "$RMSD_RECEPTOR_XVG"; then
            rm -f "$RMSD_RECEPTOR_ERR" 2>/dev/null
        else
            log_error "Failed to calculate receptor RMSD" "$RMSD_RECEPTOR_ERR"
        fi
    fi
    if file_ok "$RMSD_RECEPTOR_XVG"; then
        read -r INIT FINAL AVG MIN MAX <<< $(xvg_stats "$RMSD_RECEPTOR_XVG")
        log "Receptor RMSD (fit: receptor, calc: receptor)"
        log "  Initial: $INIT nm | Final: $FINAL nm | Avg: $AVG nm | Range: [$MIN, $MAX] nm"
    fi
    
    # 1b. Peptide Internal RMSD (peptide vs peptide)
    # Fit: Peptide | Calculate: Peptide
    # Purpose: Measure peptide conformational changes independent of its position in binding site
    # This shows if the peptide changes shape during simulation
    local RMSD_PEPTIDE_INT_XVG="${REL_ANALYSIS_DIR}/rmsd_peptide_internal.xvg"
    local RMSD_PEPTIDE_INT_ERR="${SIM_LOG_DIR}/rmsd_peptide_internal_error.log"
    
    if file_ok "$RMSD_PEPTIDE_INT_XVG"; then
        echo "Using existing peptide internal RMSD file"
    elif [ -n "$INDEX_FILE" ] && file_ok "$INDEX_FILE" && [ -n "$PEPTIDE_GROUP" ]; then
        echo "Calculating peptide internal RMSD (fit: peptide, calc: peptide)..."
        local TRAJ_FOR_ANALYSIS=$(get_traj_for_custom_index)
        
        # For internal RMSD, we fit peptide to peptide and calculate peptide RMSD
        # This removes translational/rotational motion and shows only conformational changes
        # First group = fit group (Peptide), Second group = calculation group (Peptide)
        if run_gmx_with_input "${PEPTIDE_GROUP}
${PEPTIDE_GROUP}" rms -s md.tpr -f "$TRAJ_FOR_ANALYSIS" -n "$INDEX_FILE" -o "$RMSD_PEPTIDE_INT_XVG" 2>"$RMSD_PEPTIDE_INT_ERR" && file_ok "$RMSD_PEPTIDE_INT_XVG"; then
            rm -f "$RMSD_PEPTIDE_INT_ERR" 2>/dev/null
        else
            log_error "Failed to calculate peptide internal RMSD" "$RMSD_PEPTIDE_INT_ERR"
        fi
    else
        log "Note: Peptide internal RMSD skipped (no custom index file)"
    fi
    if file_ok "$RMSD_PEPTIDE_INT_XVG"; then
        read -r INIT FINAL AVG MIN MAX <<< $(xvg_stats "$RMSD_PEPTIDE_INT_XVG")
        log "Peptide Internal RMSD (fit: peptide, calc: peptide)"
        log "  Initial: $INIT nm | Final: $FINAL nm | Avg: $AVG nm | Range: [$MIN, $MAX] nm"
        log "  Interpretation: Low values = peptide maintains conformation; High values = conformational flexibility"
    fi
    
    # 1c. Peptide Binding RMSD (receptor vs peptide) - MOST IMPORTANT FOR BINDING STUDIES
    # Fit: Receptor backbone | Calculate: Peptide
    # Purpose: Measure if peptide stays in the binding site
    # This is the KEY metric for binding stability - shows peptide movement relative to receptor
    local RMSD_PEPTIDE_XVG="${REL_ANALYSIS_DIR}/rmsd_peptide_binding.xvg"
    local RMSD_PEPTIDE_ERR="${SIM_LOG_DIR}/rmsd_peptide_binding_error.log"
    
    if file_ok "$RMSD_PEPTIDE_XVG"; then
        echo "Using existing peptide binding RMSD file"
    elif [ -n "$INDEX_FILE" ] && file_ok "$INDEX_FILE" && [ -n "$PEPTIDE_GROUP" ] && [ -n "$RECEPTOR_GROUP" ]; then
        echo "Calculating peptide binding RMSD (fit: receptor, calc: peptide)... [MOST IMPORTANT]"
        local TRAJ_FOR_ANALYSIS=$(get_traj_for_custom_index)
        
        # First group = fit group (Receptor - aligns the receptor in each frame to reference)
        # Second group = calculation group (Peptide - calculates how far peptide moved)
        # This shows peptide displacement from initial binding pose
        if run_gmx_with_input "${RECEPTOR_GROUP}
${PEPTIDE_GROUP}" rms -s md.tpr -f "$TRAJ_FOR_ANALYSIS" -n "$INDEX_FILE" -o "$RMSD_PEPTIDE_XVG" 2>"$RMSD_PEPTIDE_ERR" && file_ok "$RMSD_PEPTIDE_XVG"; then
            rm -f "$RMSD_PEPTIDE_ERR" 2>/dev/null
        else
            log_error "Failed to calculate peptide binding RMSD" "$RMSD_PEPTIDE_ERR"
        fi
    else
        log "Note: Peptide binding RMSD skipped (no custom index file or missing groups)"
    fi
    if file_ok "$RMSD_PEPTIDE_XVG"; then
        read -r INIT FINAL AVG MIN MAX <<< $(xvg_stats "$RMSD_PEPTIDE_XVG")
        log "Peptide Binding RMSD (fit: receptor, calc: peptide) *** KEY METRIC ***"
        log "  Initial: $INIT nm | Final: $FINAL nm | Avg: $AVG nm | Range: [$MIN, $MAX] nm"
        log "  Interpretation: Low values = stable binding; High values = peptide dissociating/relocating"
    fi
    
    # =========================================================================
    # ANALYSIS SECTION 2: RMSF
    # =========================================================================
    echo ""
    echo "--- RMSF Analysis ---"
    
    # Receptor RMSF
    local RMSF_RECEPTOR_XVG="${REL_ANALYSIS_DIR}/rmsf_receptor.xvg"
    local RMSF_RECEPTOR_ERR="${SIM_LOG_DIR}/rmsf_receptor_error.log"
    
    if file_ok "$RMSF_RECEPTOR_XVG"; then
        echo "Using existing receptor RMSF file"
    elif [ -n "$INDEX_FILE" ] && file_ok "$INDEX_FILE" && [ -n "$RECEPTOR_GROUP" ]; then
        echo "Calculating receptor RMSF..."
        local TRAJ_FOR_ANALYSIS=$(get_traj_for_custom_index)
        if run_gmx_with_input "$RECEPTOR_GROUP" rmsf -s md.tpr -f "$TRAJ_FOR_ANALYSIS" -n "$INDEX_FILE" -o "$RMSF_RECEPTOR_XVG" -res 2>"$RMSF_RECEPTOR_ERR" && file_ok "$RMSF_RECEPTOR_XVG"; then
            rm -f "$RMSF_RECEPTOR_ERR" 2>/dev/null
        else
            log_error "Failed to calculate receptor RMSF" "$RMSF_RECEPTOR_ERR"
        fi
    else
        echo "Calculating protein RMSF (receptor+peptide combined)..."
        if run_gmx_with_input "1" rmsf -s md.tpr -f "$ANALYSIS_XTC" -o "$RMSF_RECEPTOR_XVG" -res 2>"$RMSF_RECEPTOR_ERR" && file_ok "$RMSF_RECEPTOR_XVG"; then
            rm -f "$RMSF_RECEPTOR_ERR" 2>/dev/null
        else
            log_error "Failed to calculate protein RMSF" "$RMSF_RECEPTOR_ERR"
        fi
    fi
    if file_ok "$RMSF_RECEPTOR_XVG"; then
        read -r _ _ AVG _ MAX <<< $(xvg_stats "$RMSF_RECEPTOR_XVG")
        log "Receptor RMSF - Avg: $AVG nm | Max: $MAX nm"
    fi
    
    # Peptide RMSF
    local RMSF_PEPTIDE_XVG="${REL_ANALYSIS_DIR}/rmsf_peptide.xvg"
    local RMSF_PEPTIDE_ERR="${SIM_LOG_DIR}/rmsf_peptide_error.log"
    
    # Function to filter RMSF to only include valid peptide residues
    filter_peptide_rmsf() {
        local rmsf_file="$1"
        if [ ! -f "$rmsf_file" ]; then
            return 1
        fi
        
        local PEPLEN=$(get_peptide_len "$name")
        if [ "$PEPLEN" -le 0 ]; then
            return 1
        fi
        
        local max_reasonable=$((PEPLEN * 3))
        local filtered_file="${rmsf_file}.filtered"
        local header_file="${rmsf_file}.header"
        
        grep "^[#@]" "$rmsf_file" > "$header_file" 2>/dev/null
        
        awk -v max_res="$max_reasonable" '
            NF >= 2 && ($1+0) > 0 && ($1+0) <= max_res {
                print $0
            }
        ' "$rmsf_file" > "$filtered_file" 2>/dev/null
        
        local orig_count=$(grep -v "^[#@]" "$rmsf_file" | grep -v "^$" | wc -l)
        local filt_count=$(grep -v "^[#@]" "$filtered_file" | grep -v "^$" | wc -l)
        
        if [ "$filt_count" -lt "$orig_count" ] && [ "$filt_count" -gt 0 ]; then
            cat "$header_file" "$filtered_file" > "${filtered_file}.final" 2>/dev/null
            mv "${filtered_file}.final" "$rmsf_file"
            rm -f "$filtered_file" "$header_file" 2>/dev/null
            log "Filtered peptide RMSF: removed $((orig_count - filt_count)) outlier residue(s) (peptide length: $PEPLEN)"
            return 0
        fi
        rm -f "$filtered_file" "$header_file" "${filtered_file}.final" 2>/dev/null
        return 1
    }
    
    if file_ok "$RMSF_PEPTIDE_XVG"; then
        echo "Using existing peptide RMSF file"
        filter_peptide_rmsf "$RMSF_PEPTIDE_XVG"
    elif [ -n "$INDEX_FILE" ] && file_ok "$INDEX_FILE" && [ -n "$PEPTIDE_GROUP" ]; then
        echo "Calculating peptide RMSF..."
        local TRAJ_FOR_ANALYSIS=$(get_traj_for_custom_index)
        if run_gmx_with_input "$PEPTIDE_GROUP" rmsf -s md.tpr -f "$TRAJ_FOR_ANALYSIS" -n "$INDEX_FILE" -o "$RMSF_PEPTIDE_XVG" -res 2>"$RMSF_PEPTIDE_ERR" && file_ok "$RMSF_PEPTIDE_XVG"; then
            filter_peptide_rmsf "$RMSF_PEPTIDE_XVG"
            rm -f "$RMSF_PEPTIDE_ERR" 2>/dev/null
        else
            log "Warning: Peptide RMSF failed (peptide group may be too small for RMSF calculation)"
            log_error "Failed to calculate peptide RMSF" "$RMSF_PEPTIDE_ERR"
        fi
    else
        log "Note: Peptide RMSF skipped (no custom index file)"
    fi
    if file_ok "$RMSF_PEPTIDE_XVG"; then
        read -r _ _ AVG _ MAX <<< $(xvg_stats "$RMSF_PEPTIDE_XVG")
        log "Peptide RMSF - Avg: $AVG nm | Max: $MAX nm"
    fi
    
    # =========================================================================
    # ANALYSIS SECTION 3: Radius of Gyration
    # =========================================================================
    echo ""
    echo "--- Radius of Gyration ---"
    
    # Receptor RoG
    local ROG_REC_XVG="${REL_ANALYSIS_DIR}/rog_receptor.xvg"
    local ROG_ERR="${SIM_LOG_DIR}/rog_error.log"
    if file_ok "$ROG_REC_XVG"; then
        echo "Using existing receptor RoG file"
    elif [ -n "$INDEX_FILE" ] && file_ok "$INDEX_FILE" && [ -n "$RECEPTOR_GROUP" ]; then
        echo "Calculating receptor radius of gyration..."
        local TRAJ_FOR_ANALYSIS=$(get_traj_for_custom_index)
        if run_gmx_with_input "$RECEPTOR_GROUP" gyrate -s md.tpr -f "$TRAJ_FOR_ANALYSIS" -n "$INDEX_FILE" -o "$ROG_REC_XVG" 2>"$ROG_ERR" && file_ok "$ROG_REC_XVG"; then
            rm -f "$ROG_ERR" 2>/dev/null
        else
            log_error "Failed to calculate receptor RoG" "$ROG_ERR"
        fi
    else
        echo "Calculating protein RoG (receptor+peptide)..."
        if run_gmx_with_input "1" gyrate -s md.tpr -f "$ANALYSIS_XTC" -o "$ROG_REC_XVG" 2>"$ROG_ERR" && file_ok "$ROG_REC_XVG"; then
            rm -f "$ROG_ERR" 2>/dev/null
        else
            log_error "Failed to calculate protein RoG" "$ROG_ERR"
        fi
    fi
    if file_ok "$ROG_REC_XVG"; then
        read -r INIT FINAL AVG _ _ <<< $(xvg_stats "$ROG_REC_XVG")
        log "Receptor RoG - Initial: $INIT nm | Final: $FINAL nm | Avg: $AVG nm"
    fi
    
    # Peptide RoG
    local ROG_PEP_XVG="${REL_ANALYSIS_DIR}/rog_peptide.xvg"
    if file_ok "$ROG_PEP_XVG"; then
        echo "Using existing peptide RoG file"
    elif [ -n "$INDEX_FILE" ] && file_ok "$INDEX_FILE" && [ -n "$PEPTIDE_GROUP" ]; then
        echo "Calculating peptide radius of gyration..."
        local TRAJ_FOR_ANALYSIS=$(get_traj_for_custom_index)
        if run_gmx_with_input "$PEPTIDE_GROUP" gyrate -s md.tpr -f "$TRAJ_FOR_ANALYSIS" -n "$INDEX_FILE" -o "$ROG_PEP_XVG" 2>"$ROG_ERR" && file_ok "$ROG_PEP_XVG"; then
            rm -f "$ROG_ERR" 2>/dev/null
        else
            log_error "Failed to calculate peptide RoG" "$ROG_ERR"
        fi
    else
        log "Note: Peptide RoG skipped (no custom index file)"
    fi
    if file_ok "$ROG_PEP_XVG"; then
        read -r INIT FINAL AVG _ _ <<< $(xvg_stats "$ROG_PEP_XVG")
        log "Peptide RoG - Initial: $INIT nm | Final: $FINAL nm | Avg: $AVG nm"
    fi
    
    # =========================================================================
    # ANALYSIS SECTION 4: Hydrogen Bonds
    # =========================================================================
    local HBONDS_XVG="${REL_ANALYSIS_DIR}/hbonds.xvg"
    local HBOND_ERR="${SIM_LOG_DIR}/hbond_error.log"
    echo ""
    echo "--- Hydrogen Bonds (Peptide-Receptor) ---"
    if file_ok "$HBONDS_XVG"; then
        echo "Using existing hydrogen bonds file"
    elif [ -n "$INDEX_FILE" ] && file_ok "$INDEX_FILE" && [ -n "$PEPTIDE_GROUP" ] && [ -n "$RECEPTOR_GROUP" ]; then
        echo "Calculating intermolecular hydrogen bonds..."
        local TRAJ_FOR_ANALYSIS=$(get_traj_for_custom_index)
        if run_gmx_with_input "${PEPTIDE_GROUP}
${RECEPTOR_GROUP}" hbond -s md.tpr -f "$TRAJ_FOR_ANALYSIS" -n "$INDEX_FILE" -num "$HBONDS_XVG" 2>"$HBOND_ERR" && file_ok "$HBONDS_XVG"; then
            rm -f "$HBOND_ERR" 2>/dev/null
        else
            log_error "Failed to calculate hydrogen bonds" "$HBOND_ERR"
        fi
    else
        log "Note: Hydrogen bond analysis skipped (no custom index file)"
    fi
    if file_ok "$HBONDS_XVG"; then
        read -r INIT FINAL AVG MIN MAX <<< $(xvg_stats "$HBONDS_XVG")
        log "Hydrogen Bonds - Initial: ${INIT%.*} | Final: ${FINAL%.*} | Avg: $AVG | Max: ${MAX%.*}"
    fi
    
    # =========================================================================
    # FINALIZE
    # =========================================================================
    log ""
    log "=========================================="
    log "Analysis completed: $(date)"
    log ""
    log "Output files saved to: ${SIM_ANALYSIS_DIR}/"
    log "  - rmsd_receptor_backbone.xvg  : Receptor stability"
    log "  - rmsd_peptide_internal.xvg   : Peptide conformational changes"
    log "  - rmsd_peptide_binding.xvg    : Peptide binding stability (KEY METRIC)"
    log "  - rmsf_receptor.xvg           : Receptor flexibility per residue"
    log "  - rmsf_peptide.xvg            : Peptide flexibility per residue"
    log "  - rog_receptor.xvg            : Receptor compactness"
    log "  - rog_peptide.xvg             : Peptide compactness"
    log "  - hbonds.xvg                  : Peptide-receptor hydrogen bonds"
    
    # Copy results to output directory
    for xvg_file in "${SIM_ANALYSIS_DIR}"/*.xvg; do
        [ -f "$xvg_file" ] && cp "$xvg_file" "${OUTPUT_DIR}/${name}_$(basename "$xvg_file")" 2>/dev/null
    done
    
    # Cleanup temp files
    rm -f .gmx_input*.txt .gmx_wrapper.sh .pbc_*.xvg .pbc_*.log .pbc_input.txt 2>/dev/null
    
    cd "$ORIGINAL_DIR"
}

# =============================================================================
# MAIN EXECUTION
# =============================================================================
setup_gromacs

if [ $# -gt 0 ]; then
    for dir_arg in "$@"; do
        analyze_directory "$dir_arg"
    done
else
    for dir in *_docked_*/; do
        [ -d "$dir" ] && analyze_directory "$dir"
    done
fi

echo ""
echo "=========================================="
echo "Analysis complete!"
echo "Results in ${MAIN_DIR}/ within each analyzed directory"
echo ""
echo "RMSD Interpretation Guide:"
echo "  - Receptor RMSD: <0.3 nm stable, 0.3-0.5 nm moderate flexibility, >0.5 nm significant motion"
echo "  - Peptide Internal RMSD: Low = maintains conformation, High = flexible/changes shape"
echo "  - Peptide Binding RMSD: <0.3 nm stable binding, 0.3-0.5 nm minor shifts, >0.5 nm significant movement"
echo ""
echo "For extended analyses (DSSP, MinDist), run:"
echo "  ./analyze_md_extended.sh <directory>"
echo "=========================================="
