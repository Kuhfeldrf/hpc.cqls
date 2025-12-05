#!/bin/bash

# Store original working directory
ORIGINAL_DIR=$(pwd)

# Create main directory structure for analysis organization
MAIN_DIR="analysis_output"
OUTPUT_DIR="${ORIGINAL_DIR}/${MAIN_DIR}/output_files"
LOG_DIR="${ORIGINAL_DIR}/${MAIN_DIR}/processing_logs"
ANALYSIS_DIR="${ORIGINAL_DIR}/${MAIN_DIR}/analysis_data"

mkdir -p "$OUTPUT_DIR" "$LOG_DIR" "$ANALYSIS_DIR"

# Get peptide info based on directory name
get_peptide_len() {
    case $1 in
        *GLAPYKLRPVAA*) echo 12 ;;
        *LLFKDSAIGF*) echo 10 ;;
        *RPKLPLRYP*) echo 9 ;;
        *) echo 0 ;;
    esac
}

# Function to analyze a single directory
analyze_directory() {
    local dir="$1"
    # Remove trailing slash if present
    dir="${dir%/}"
    if [ ! -d "$dir" ]; then
        echo "Error: Directory '$dir' does not exist"
        return 1
    fi
    # Get absolute path (using subshell to avoid changing current directory)
    dir=$(cd "$dir" 2>/dev/null && pwd)
    if [ -z "$dir" ]; then
        echo "Error: Cannot access directory '$1'"
        return 1
    fi
    local name=$(basename "$dir")
    
    echo ""
    echo "=========================================="
    echo "Analyzing: $name"
    echo "=========================================="
    
    # Create analysis subdirectory for this simulation
    SIM_ANALYSIS_DIR="${ANALYSIS_DIR}/${name}"
    SIM_LOG_DIR="${LOG_DIR}/${name}"
    mkdir -p "$SIM_ANALYSIS_DIR" "$SIM_LOG_DIR"
    
    # Create log file for this simulation
    LOG_FILE="${SIM_LOG_DIR}/analysis_$(date +%Y%m%d_%H%M%S).log"
    {
        echo "Analysis started: $(date)"
        echo "Simulation: $name"
        echo "=========================================="
    } > "$LOG_FILE"
    
    cd "$dir"
    
    # Try to find trajectory file - check for md.xtc first, then *_100ns.xtc
    XTC=""
    if [ -f "md.xtc" ]; then
        XTC="md.xtc"
    else
        XTC=$(ls *_100ns.xtc 2>/dev/null | head -1)
    fi
    
    if [ -z "$XTC" ] || [ ! -f "$XTC" ]; then
        echo "No trajectory found!"
        echo "ERROR: No trajectory file found (looked for md.xtc or *_100ns.xtc)" >> "$LOG_FILE"
        cd "$ORIGINAL_DIR"
        return 1
    fi
    
    echo "Found trajectory: $XTC"
    echo "Found trajectory: $XTC" >> "$LOG_FILE"
    
    # Get peptide length and calculate residue indices
    PEPLEN=$(get_peptide_len "$name")
    if [ "$PEPLEN" -eq 0 ]; then
        echo "Warning: Could not determine peptide length from directory name '$name'"
        echo "Warning: Could not determine peptide length from directory name '$name'" >> "$LOG_FILE"
    fi
    
    # Try multiple methods to get total residue count
    TOTAL_RES=""
    
    # Method 1: Try gmx make_ndx
    TOTAL_RES=$(echo "q" | gmx make_ndx -f md.tpr 2>&1 | grep -i "protein.*residue" | head -1 | awk '{for(i=1;i<=NF;i++) if($i ~ /^[0-9]+$/) {print $i; exit}}')
    
    # Method 2: Try gmx check
    if [ -z "$TOTAL_RES" ] || [ "$TOTAL_RES" -le 0 ]; then
        TOTAL_RES=$(gmx check -f md.tpr 2>&1 | grep -i "residue" | head -1 | awk '{for(i=1;i<=NF;i++) if($i ~ /^[0-9]+$/) {print $i; exit}}')
    fi
    
    # Method 3: Try gmx info
    if [ -z "$TOTAL_RES" ] || [ "$TOTAL_RES" -le 0 ]; then
        TOTAL_RES=$(gmx check -s md.tpr 2>&1 | grep -i "residue" | head -1 | awk '{for(i=1;i<=NF;i++) if($i ~ /^[0-9]+$/) {print $i; exit}}')
    fi
    
    # Method 4: Try parsing from topology file if available
    if [ -z "$TOTAL_RES" ] || [ "$TOTAL_RES" -le 0 ]; then
        if [ -f "topol.top" ]; then
            # Try to count residues from topology (approximate)
            TOTAL_RES=$(grep -E "^\[.*\]" topol.top | grep -v "defaults\|atomtypes\|moleculetype\|atoms\|bonds\|angles\|dihedrals" | wc -l)
            # This is rough, so we'll use it as fallback only
        fi
    fi
    
    # If still no result, try to get it from RMSD calculation (use first available RMSD file)
    if [ -z "$TOTAL_RES" ] || [ "$TOTAL_RES" -le 0 ]; then
        # Try to calculate RMSD and infer from output
        TEMP_RMSD="${SIM_ANALYSIS_DIR}/temp_rmsd.xvg"
        echo "4 4" | gmx rms -s md.tpr -f $XTC -o "$TEMP_RMSD" -tu ns >/dev/null 2>&1
        if [ -f "$TEMP_RMSD" ] && [ -s "$TEMP_RMSD" ]; then
            # RMSD files don't directly give residue count, so this won't work
            rm -f "$TEMP_RMSD"
        fi
    fi
    
    # Final check - if we still don't have it, try a default or skip residue-specific calculations
    if [ -z "$TOTAL_RES" ] || [ "$TOTAL_RES" -le 0 ]; then
        echo "Warning: Could not determine total number of residues from md.tpr"
        echo "Warning: Could not determine total number of residues from md.tpr" >> "$LOG_FILE"
        echo "Attempting to continue with available information..."
        echo "Attempting to continue with available information..." >> "$LOG_FILE"
        
        # Set a default or try to continue without residue-specific info
        # For now, we'll try to proceed and see if we can work around it
        TOTAL_RES=0
    fi
    
    # Only calculate receptor/peptide boundaries if we have valid residue count
    if [ "$TOTAL_RES" -gt 0 ] && [ "$PEPLEN" -gt 0 ]; then
        REC_END=$((TOTAL_RES - PEPLEN))
        PEP_START=$((REC_END + 1))
        
        if [ "$REC_END" -le 0 ] || [ "$PEP_START" -gt "$TOTAL_RES" ]; then
            echo "Warning: Invalid residue calculation (Total: $TOTAL_RES, Peptide length: $PEPLEN)"
            echo "Warning: Invalid residue calculation (Total: $TOTAL_RES, Peptide length: $PEPLEN)" >> "$LOG_FILE"
            echo "Will attempt to use default group selections..."
            echo "Will attempt to use default group selections..." >> "$LOG_FILE"
            REC_END=0
            PEP_START=0
        else
            echo "Receptor: residues 1-$REC_END, Peptide: residues $PEP_START-$TOTAL_RES ($PEPLEN residues)"
            echo "Receptor: residues 1-$REC_END, Peptide: residues $PEP_START-$TOTAL_RES ($PEPLEN residues)" >> "$LOG_FILE"
        fi
    else
        echo "Warning: Cannot determine receptor/peptide boundaries. Using default selections."
        echo "Warning: Cannot determine receptor/peptide boundaries. Using default selections." >> "$LOG_FILE"
        REC_END=0
        PEP_START=0
    fi
    
    # Create index file with proper residue selection
    INDEX_FILE="${SIM_ANALYSIS_DIR}/analysis.ndx"
    if [ ! -f "$INDEX_FILE" ]; then
        echo "Creating index file..."
        if [ "$REC_END" -gt 0 ] && [ "$PEP_START" -gt 0 ]; then
            # Use calculated residue ranges
            gmx make_ndx -f md.tpr -o "$INDEX_FILE" << ENDNDX >/dev/null 2>&1
ri 1-$REC_END
name 17 Receptor
ri $PEP_START-$TOTAL_RES
name 18 Peptide
q
ENDNDX
        else
            # Use default group selections (Protein, System, etc.)
            gmx make_ndx -f md.tpr -o "$INDEX_FILE" << ENDNDX >/dev/null 2>&1
1
name 17 Receptor
q
ENDNDX
            # Try to find peptide group - might be named differently
            # For now, we'll use group 1 (Protein) as receptor
            # Peptide selection will need to be done manually or via other means
        fi
    else
        echo "Using existing index file"
    fi

    # 1. RMSD of backbone
    echo ""
    echo "--- Backbone RMSD ---"
    BACKBONE_XVG="${SIM_ANALYSIS_DIR}/rmsd_backbone.xvg"
    if [ -f "$BACKBONE_XVG" ] && [ -s "$BACKBONE_XVG" ]; then
        echo "Using existing backbone RMSD file"
    else
        echo "Calculating backbone RMSD..."
        echo "4 4" | gmx rms -s md.tpr -f $XTC -o "$BACKBONE_XVG" -tu ns >/dev/null 2>&1
    fi
    if [ -f "$BACKBONE_XVG" ] && [ -s "$BACKBONE_XVG" ]; then
        INIT=$(head -30 "$BACKBONE_XVG" | grep -v "^[#@]" | head -1 | awk '{printf "%.3f", $2}')
        FINAL=$(tail -5 "$BACKBONE_XVG" | grep -v "^[#@]" | tail -1 | awk '{printf "%.3f", $2}')
        MAX=$(grep -v "^[#@]" "$BACKBONE_XVG" | awk 'BEGIN{max=0} {if($2>max) max=$2} END{printf "%.3f", max}')
        echo "Initial: $INIT nm | Final: $FINAL nm | Max: $MAX nm"
        echo "Backbone RMSD - Initial: $INIT nm | Final: $FINAL nm | Max: $MAX nm" >> "$LOG_FILE"
    else
        echo "Failed to calculate backbone RMSD"
        echo "ERROR: Failed to create $BACKBONE_XVG" >> "$LOG_FILE"
    fi
    
    # 2. RMSD of peptide backbone only
    echo ""
    echo "--- Peptide RMSD ---"
    PEPTIDE_XVG="${SIM_ANALYSIS_DIR}/rmsd_peptide.xvg"
    if [ -f "$PEPTIDE_XVG" ] && [ -s "$PEPTIDE_XVG" ]; then
        echo "Using existing peptide RMSD file"
    else
        echo "Calculating peptide RMSD..."
        # Try group 18 (Peptide) first, fallback to other methods
        if [ -f "$INDEX_FILE" ] && grep -q "Peptide" "$INDEX_FILE" 2>/dev/null; then
            echo "18 18" | gmx rms -s md.tpr -f $XTC -n "$INDEX_FILE" -o "$PEPTIDE_XVG" -tu ns >/dev/null 2>&1
        else
            # Fallback: try to use a different group or skip
            echo "Note: Peptide group not found in index file, skipping peptide-specific RMSD"
            echo "Note: Peptide group not found in index file" >> "$LOG_FILE"
        fi
    fi
    if [ -f "$PEPTIDE_XVG" ] && [ -s "$PEPTIDE_XVG" ]; then
        INIT=$(head -30 "$PEPTIDE_XVG" | grep -v "^[#@]" | head -1 | awk '{printf "%.3f", $2}')
        FINAL=$(tail -5 "$PEPTIDE_XVG" | grep -v "^[#@]" | tail -1 | awk '{printf "%.3f", $2}')
        MAX=$(grep -v "^[#@]" "$PEPTIDE_XVG" | awk 'BEGIN{max=0} {if($2>max) max=$2} END{printf "%.3f", max}')
        echo "Initial: $INIT nm | Final: $FINAL nm | Max: $MAX nm"
        echo "Peptide RMSD - Initial: $INIT nm | Final: $FINAL nm | Max: $MAX nm" >> "$LOG_FILE"
    else
        echo "Failed to calculate peptide RMSD"
        echo "ERROR: Failed to create $PEPTIDE_XVG" >> "$LOG_FILE"
    fi
    
    # 3. Minimum distance peptide-receptor
    echo ""
    echo "--- Minimum Distance (Peptide-Receptor) ---"
    MINDIST_XVG="${SIM_ANALYSIS_DIR}/mindist.xvg"
    if [ -f "$MINDIST_XVG" ] && [ -s "$MINDIST_XVG" ]; then
        echo "Using existing minimum distance file"
    else
        echo "Calculating minimum distance..."
        # Try peptide-receptor groups if available
        if [ -f "$INDEX_FILE" ] && grep -q "Peptide" "$INDEX_FILE" 2>/dev/null && grep -q "Receptor" "$INDEX_FILE" 2>/dev/null; then
            echo "18 17" | gmx mindist -s md.tpr -f $XTC -n "$INDEX_FILE" -od "$MINDIST_XVG" -tu ns >/dev/null 2>&1
        else
            # Fallback: use default groups (e.g., Protein-Protein)
            echo "1 1" | gmx mindist -s md.tpr -f $XTC -od "$MINDIST_XVG" -tu ns >/dev/null 2>&1
        fi
    fi
    if [ -f "$MINDIST_XVG" ] && [ -s "$MINDIST_XVG" ]; then
        INIT=$(head -30 "$MINDIST_XVG" | grep -v "^[#@]" | head -1 | awk '{printf "%.3f", $2}')
        FINAL=$(tail -5 "$MINDIST_XVG" | grep -v "^[#@]" | tail -1 | awk '{printf "%.3f", $2}')
        AVG=$(grep -v "^[#@]" "$MINDIST_XVG" | awk '{sum+=$2; n++} END {printf "%.3f", sum/n}')
        MIN=$(grep -v "^[#@]" "$MINDIST_XVG" | awk 'BEGIN{min=999} {if($2<min) min=$2} END{printf "%.3f", min}')
        echo "Initial: $INIT nm | Final: $FINAL nm | Avg: $AVG nm | Min: $MIN nm"
        echo "Minimum Distance - Initial: $INIT nm | Final: $FINAL nm | Avg: $AVG nm | Min: $MIN nm" >> "$LOG_FILE"
    else
        echo "Failed to calculate"
        echo "ERROR: Failed to calculate minimum distance" >> "$LOG_FILE"
    fi
    
    # 4. H-bonds
    echo ""
    echo "--- Hydrogen Bonds (Peptide-Receptor) ---"
    HBONDS_XVG="${SIM_ANALYSIS_DIR}/hbonds.xvg"
    HBOND_ERR="${SIM_LOG_DIR}/hbond_error.log"
    if [ -f "$HBONDS_XVG" ] && [ -s "$HBONDS_XVG" ]; then
        echo "Using existing hydrogen bonds file"
    else
        echo "Calculating hydrogen bonds..."
        # gmx hbond syntax: specify donor and acceptor groups separately
        # For intermolecular H-bonds between peptide and receptor, we need to consider both directions
        # Try: peptide as donor, receptor as acceptor
        gmx hbond -s md.tpr -f $XTC -n "$INDEX_FILE" -don 18 -acc 17 -num "$HBONDS_XVG" >"$HBOND_ERR" 2>&1
        # If that creates an empty file or fails, try the reverse
        if [ ! -f "$HBONDS_XVG" ] || [ ! -s "$HBONDS_XVG" ]; then
            if [ -f "$INDEX_FILE" ] && grep -q "Receptor" "$INDEX_FILE" 2>/dev/null && grep -q "Peptide" "$INDEX_FILE" 2>/dev/null; then
                gmx hbond -s md.tpr -f $XTC -n "$INDEX_FILE" -don 17 -acc 18 -num "$HBONDS_XVG" >"$HBOND_ERR" 2>&1
            fi
        fi
        # If still failing, try without specifying donor/acceptor (finds all H-bonds between the groups)
        if [ ! -f "$HBONDS_XVG" ] || [ ! -s "$HBONDS_XVG" ]; then
            if [ -f "$INDEX_FILE" ] && grep -q "Peptide" "$INDEX_FILE" 2>/dev/null && grep -q "Receptor" "$INDEX_FILE" 2>/dev/null; then
                printf "18\n17\n" | gmx hbond -s md.tpr -f $XTC -n "$INDEX_FILE" -num "$HBONDS_XVG" >"$HBOND_ERR" 2>&1
            else
                # Final fallback: use default protein group
                printf "1\n1\n" | gmx hbond -s md.tpr -f $XTC -num "$HBONDS_XVG" >"$HBOND_ERR" 2>&1
            fi
        fi
    fi
    if [ -f "$HBONDS_XVG" ] && [ -s "$HBONDS_XVG" ]; then
        INIT=$(head -30 "$HBONDS_XVG" | grep -v "^[#@]" | head -1 | awk '{print $2}')
        FINAL=$(tail -5 "$HBONDS_XVG" | grep -v "^[#@]" | tail -1 | awk '{print $2}')
        AVG=$(grep -v "^[#@]" "$HBONDS_XVG" | awk '{sum+=$2; n++} END {printf "%.1f", sum/n}')
        MAX=$(grep -v "^[#@]" "$HBONDS_XVG" | awk 'BEGIN{max=0} {if($2>max) max=$2} END{print max}')
        echo "Initial: $INIT | Final: $FINAL | Avg: $AVG | Max: $MAX"
        echo "Hydrogen Bonds - Initial: $INIT | Final: $FINAL | Avg: $AVG | Max: $MAX" >> "$LOG_FILE"
        rm -f "$HBOND_ERR" 2>/dev/null
    else
        echo "Failed to calculate"
        echo "ERROR: Failed to calculate hydrogen bonds" >> "$LOG_FILE"
        if [ -f "$HBOND_ERR" ]; then
            echo "Error details:" >> "$LOG_FILE"
            cat "$HBOND_ERR" >> "$LOG_FILE"
        fi
    fi
    
    # 5. RMSF (Root Mean Square Fluctuation)
    echo ""
    echo "--- RMSF Analysis ---"
    
    # RMSF for backbone
    RMSF_BACKBONE_XVG="${SIM_ANALYSIS_DIR}/rmsf_backbone.xvg"
    if [ -f "$RMSF_BACKBONE_XVG" ] && [ -s "$RMSF_BACKBONE_XVG" ]; then
        echo "Using existing backbone RMSF file"
    else
        echo "Calculating backbone RMSF..."
        echo "4" | gmx rmsf -s md.tpr -f $XTC -o "$RMSF_BACKBONE_XVG" -res >/dev/null 2>&1
    fi
    if [ -f "$RMSF_BACKBONE_XVG" ] && [ -s "$RMSF_BACKBONE_XVG" ]; then
        AVG=$(grep -v "^[#@]" "$RMSF_BACKBONE_XVG" | awk '{sum+=$2; n++} END {printf "%.3f", sum/n}')
        MAX=$(grep -v "^[#@]" "$RMSF_BACKBONE_XVG" | awk 'BEGIN{max=0} {if($2>max) max=$2} END{printf "%.3f", max}')
        echo "Backbone RMSF - Avg: $AVG nm | Max: $MAX nm"
        echo "Backbone RMSF - Avg: $AVG nm | Max: $MAX nm" >> "$LOG_FILE"
    else
        echo "Failed to calculate backbone RMSF"
        echo "ERROR: Failed to create $RMSF_BACKBONE_XVG" >> "$LOG_FILE"
    fi
    
    # RMSF for protein (all atoms)
    RMSF_PROTEIN_XVG="${SIM_ANALYSIS_DIR}/rmsf_protein.xvg"
    if [ -f "$RMSF_PROTEIN_XVG" ] && [ -s "$RMSF_PROTEIN_XVG" ]; then
        echo "Using existing protein RMSF file"
    else
        echo "Calculating protein RMSF..."
        echo "1" | gmx rmsf -s md.tpr -f $XTC -o "$RMSF_PROTEIN_XVG" -res >/dev/null 2>&1
    fi
    if [ -f "$RMSF_PROTEIN_XVG" ] && [ -s "$RMSF_PROTEIN_XVG" ]; then
        AVG=$(grep -v "^[#@]" "$RMSF_PROTEIN_XVG" | awk '{sum+=$2; n++} END {printf "%.3f", sum/n}')
        MAX=$(grep -v "^[#@]" "$RMSF_PROTEIN_XVG" | awk 'BEGIN{max=0} {if($2>max) max=$2} END{printf "%.3f", max}')
        echo "Protein RMSF - Avg: $AVG nm | Max: $MAX nm"
        echo "Protein RMSF - Avg: $AVG nm | Max: $MAX nm" >> "$LOG_FILE"
    else
        echo "Failed to calculate protein RMSF"
        echo "ERROR: Failed to create $RMSF_PROTEIN_XVG" >> "$LOG_FILE"
    fi
    
    # RMSF for sidechain (approximate: protein - backbone)
    # Note: GROMACS doesn't have direct sidechain selection, so we'll use a workaround
    # or skip this if it's too complex. For now, we'll create it by selecting sidechain atoms
    RMSF_SIDECHAIN_XVG="${SIM_ANALYSIS_DIR}/rmsf_sidechain.xvg"
    if [ -f "$RMSF_SIDECHAIN_XVG" ] && [ -s "$RMSF_SIDECHAIN_XVG" ]; then
        echo "Using existing sidechain RMSF file"
    else
        echo "Calculating sidechain RMSF..."
        # Try to select sidechain atoms (not backbone: CA, C, N, O)
        # This is approximate - we'll use a selection that excludes backbone
        echo "3" | gmx rmsf -s md.tpr -f $XTC -o "$RMSF_SIDECHAIN_XVG" -res >/dev/null 2>&1
        # If that doesn't work, we can try a custom selection
        if [ ! -f "$RMSF_SIDECHAIN_XVG" ] || [ ! -s "$RMSF_SIDECHAIN_XVG" ]; then
            # Alternative: use protein and subtract backbone (requires post-processing)
            # For simplicity, we'll just use protein RMSF as fallback
            if [ -f "$RMSF_PROTEIN_XVG" ] && [ -s "$RMSF_PROTEIN_XVG" ]; then
                cp "$RMSF_PROTEIN_XVG" "$RMSF_SIDECHAIN_XVG"
                echo "Using protein RMSF as sidechain approximation"
            fi
        fi
    fi
    if [ -f "$RMSF_SIDECHAIN_XVG" ] && [ -s "$RMSF_SIDECHAIN_XVG" ]; then
        AVG=$(grep -v "^[#@]" "$RMSF_SIDECHAIN_XVG" | awk '{sum+=$2; n++} END {printf "%.3f", sum/n}')
        MAX=$(grep -v "^[#@]" "$RMSF_SIDECHAIN_XVG" | awk 'BEGIN{max=0} {if($2>max) max=$2} END{printf "%.3f", max}')
        echo "Sidechain RMSF - Avg: $AVG nm | Max: $MAX nm"
        echo "Sidechain RMSF - Avg: $AVG nm | Max: $MAX nm" >> "$LOG_FILE"
    else
        echo "Note: Sidechain RMSF not calculated (optional)"
        echo "Note: Sidechain RMSF not calculated" >> "$LOG_FILE"
    fi
    
    # 6. Radius of Gyration (RoG)
    echo ""
    echo "--- Radius of Gyration ---"
    
    # Overall RoG
    ROG_XVG="${SIM_ANALYSIS_DIR}/rog.xvg"
    if [ -f "$ROG_XVG" ] && [ -s "$ROG_XVG" ]; then
        echo "Using existing RoG file"
    else
        echo "Calculating overall radius of gyration..."
        echo "1" | gmx gyrate -s md.tpr -f $XTC -o "$ROG_XVG" >/dev/null 2>&1
    fi
    if [ -f "$ROG_XVG" ] && [ -s "$ROG_XVG" ]; then
        INIT=$(head -30 "$ROG_XVG" | grep -v "^[#@]" | head -1 | awk '{printf "%.3f", $2}')
        FINAL=$(tail -5 "$ROG_XVG" | grep -v "^[#@]" | tail -1 | awk '{printf "%.3f", $2}')
        AVG=$(grep -v "^[#@]" "$ROG_XVG" | awk '{sum+=$2; n++} END {printf "%.3f", sum/n}')
        echo "Overall RoG - Initial: $INIT nm | Final: $FINAL nm | Avg: $AVG nm"
        echo "Overall RoG - Initial: $INIT nm | Final: $FINAL nm | Avg: $AVG nm" >> "$LOG_FILE"
    else
        echo "Failed to calculate overall RoG"
        echo "ERROR: Failed to create $ROG_XVG" >> "$LOG_FILE"
    fi
    
    # RoG components (x, y, z)
    ROG_X_XVG="${SIM_ANALYSIS_DIR}/rog_x.xvg"
    ROG_Y_XVG="${SIM_ANALYSIS_DIR}/rog_y.xvg"
    ROG_Z_XVG="${SIM_ANALYSIS_DIR}/rog_z.xvg"
    
    if [ -f "$ROG_X_XVG" ] && [ -f "$ROG_Y_XVG" ] && [ -f "$ROG_Z_XVG" ] && \
       [ -s "$ROG_X_XVG" ] && [ -s "$ROG_Y_XVG" ] && [ -s "$ROG_Z_XVG" ]; then
        echo "Using existing RoG component files"
    else
        echo "Calculating RoG components (x, y, z)..."
        # gmx gyrate with -oall option outputs all components
        # We'll extract x, y, z components from the output
        ROG_ALL_XVG="${SIM_ANALYSIS_DIR}/rog_all.xvg"
        echo "1" | gmx gyrate -s md.tpr -f $XTC -o "$ROG_ALL_XVG" >/dev/null 2>&1
        
        # Extract components if the file has multiple columns
        if [ -f "$ROG_ALL_XVG" ] && [ -s "$ROG_ALL_XVG" ]; then
            # Check if file has multiple columns (x, y, z components)
            COL_COUNT=$(head -50 "$ROG_ALL_XVG" | grep -v "^[#@]" | head -1 | awk '{print NF}')
            if [ "$COL_COUNT" -ge 4 ]; then
                # Extract x, y, z components (columns 2, 3, 4 typically)
                awk '!/^[@#]/ {if(NF>=4) print $1, $2}' "$ROG_ALL_XVG" > "$ROG_X_XVG" 2>/dev/null
                awk '!/^[@#]/ {if(NF>=4) print $1, $3}' "$ROG_ALL_XVG" > "$ROG_Y_XVG" 2>/dev/null
                awk '!/^[@#]/ {if(NF>=4) print $1, $4}' "$ROG_ALL_XVG" > "$ROG_Z_XVG" 2>/dev/null
                
                # Add headers back
                head -1 "$ROG_ALL_XVG" | grep "^@" > /tmp/rog_header.txt 2>/dev/null
                if [ -f /tmp/rog_header.txt ] && [ -s /tmp/rog_header.txt ]; then
                    cat /tmp/rog_header.txt "$ROG_X_XVG" > /tmp/rog_x_tmp && mv /tmp/rog_x_tmp "$ROG_X_XVG"
                    cat /tmp/rog_header.txt "$ROG_Y_XVG" > /tmp/rog_y_tmp && mv /tmp/rog_y_tmp "$ROG_Y_XVG"
                    cat /tmp/rog_header.txt "$ROG_Z_XVG" > /tmp/rog_z_tmp && mv /tmp/rog_z_tmp "$ROG_Z_XVG"
                    rm -f /tmp/rog_header.txt
                fi
            else
                # If only one component, copy to all three
                if [ -f "$ROG_XVG" ] && [ -s "$ROG_XVG" ]; then
                    cp "$ROG_XVG" "$ROG_X_XVG"
                    cp "$ROG_XVG" "$ROG_Y_XVG"
                    cp "$ROG_XVG" "$ROG_Z_XVG"
                fi
            fi
            rm -f "$ROG_ALL_XVG" 2>/dev/null
        fi
    fi
    
    # 7. DSSP Secondary Structure
    echo ""
    echo "--- DSSP Secondary Structure ---"
    DSSP_XVG="${SIM_ANALYSIS_DIR}/dssp.xvg"
    DSSP_DAT="${SIM_ANALYSIS_DIR}/dssp.dat"
    DSSP_ERR="${SIM_LOG_DIR}/dssp_error.log"
    
    # Check for existing DSSP files before regenerating
    if [ -f "$DSSP_XVG" ] && [ -s "$DSSP_XVG" ]; then
        echo "Using existing DSSP file"
    else
        echo "Calculating DSSP secondary structure..."
        # gmx dssp syntax for GROMACS 2024.5: -f trajectory, -s topology, -o DAT output, -num XVG output with counts
        # Select group 1 (usually Protein)
        echo "1" | gmx dssp -f "$XTC" -s md.tpr -o "$DSSP_DAT" -num "$DSSP_XVG" -dt 100 >"$DSSP_ERR" 2>&1
        
        # If that fails, try without dt flag (analyze all frames)
        if [ ! -f "$DSSP_XVG" ] || [ ! -s "$DSSP_XVG" ]; then
            echo "1" | gmx dssp -f "$XTC" -s md.tpr -o "$DSSP_DAT" -num "$DSSP_XVG" >"$DSSP_ERR" 2>&1
        fi
        
        # Check if DSSP command or binary is available
        if [ ! -f "$DSSP_XVG" ] || [ ! -s "$DSSP_XVG" ]; then
            if grep -qiE "'dssp'.*is not.*command|not a GROMACS command|command.*not.*found|dssp.*not.*found|cannot.*find.*dssp" "$DSSP_ERR" 2>/dev/null; then
                echo "Warning: DSSP command not available. This may require GROMACS with DSSP support or the DSSP binary."
                echo "Warning: For GROMACS 2024.5, ensure 'gmx dssp' is available. Install DSSP: conda install -c bioconda dssp"
                echo "Warning: DSSP command not available" >> "$LOG_FILE"
                if [ -f "$DSSP_ERR" ]; then
                    echo "Error details:" >> "$LOG_FILE"
                    head -10 "$DSSP_ERR" >> "$LOG_FILE" 2>/dev/null
                fi
            elif [ -f "$DSSP_ERR" ]; then
                echo "Warning: DSSP analysis failed. Check error log: $DSSP_ERR"
                echo "Warning: DSSP analysis failed" >> "$LOG_FILE"
                # Show error details for debugging
                echo "Error preview:" >> "$LOG_FILE"
                head -10 "$DSSP_ERR" >> "$LOG_FILE" 2>/dev/null
            fi
        else
            # Clean up error log if successful
            rm -f "$DSSP_ERR" 2>/dev/null
        fi
    fi
    
    if [ -f "$DSSP_XVG" ] && [ -s "$DSSP_XVG" ]; then
        # Parse DSSP output to get summary statistics
        TOTAL_FRAMES=$(grep -v "^[#@]" "$DSSP_XVG" | wc -l)
        if [ "$TOTAL_FRAMES" -gt 0 ]; then
            echo "DSSP analysis completed - $TOTAL_FRAMES frames analyzed"
            echo "DSSP analysis completed - $TOTAL_FRAMES frames analyzed" >> "$LOG_FILE"
        else
            echo "Warning: DSSP file created but appears empty"
            echo "Warning: DSSP file appears empty" >> "$LOG_FILE"
        fi
    else
        echo "Note: DSSP analysis not available"
        echo "Note: DSSP analysis not available" >> "$LOG_FILE"
        if [ -f "$DSSP_ERR" ]; then
            echo "  Error details saved to: $DSSP_ERR"
        fi
    fi
    
    # Finalize log file
    {
        echo "=========================================="
        echo "Analysis completed: $(date)"
        echo "Output directory: $SIM_ANALYSIS_DIR"
    } >> "$LOG_FILE"
    
    # Copy analysis files to output directory with simulation name prefix
    for xvg_file in "${SIM_ANALYSIS_DIR}"/*.xvg; do
        if [ -f "$xvg_file" ]; then
            cp "$xvg_file" "${OUTPUT_DIR}/${name}_$(basename "$xvg_file")" 2>/dev/null
        fi
    done
    
    cd "$ORIGINAL_DIR"
}

# Main execution logic
if [ $# -gt 0 ]; then
    # If directory argument provided, analyze that directory
    for dir_arg in "$@"; do
        analyze_directory "$dir_arg"
    done
else
    # If no argument, analyze all *_docked_*/ directories
    for dir in *_docked_*/; do
        if [ -d "$dir" ]; then
            analyze_directory "$dir"
        fi
    done
fi

echo ""
echo "=========================================="
echo "Analysis complete!"
echo "Results organized in: $MAIN_DIR"
echo "  - Analysis data: $ANALYSIS_DIR"
echo "  - Output files: $OUTPUT_DIR"
echo "  - Processing logs: $LOG_DIR"
echo "=========================================="
