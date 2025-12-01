#!/bin/bash

# Get peptide info based on directory name
get_peptide_len() {
    case $1 in
        *GLAPYKLRPVAA*) echo 12 ;;
        *LLFKDSAIGF*) echo 10 ;;
        *RPKLPLRYP*) echo 9 ;;
        *) echo 0 ;;
    esac
}

for dir in *_docked_*/; do
    echo ""
    echo "=========================================="
    name=$(basename "$dir")
    echo "Analyzing: $name"
    echo "=========================================="
    cd "$dir"
    
    XTC=$(ls *_100ns.xtc 2>/dev/null)
    if [ -z "$XTC" ]; then
        echo "No trajectory found!"
        cd ..
        continue
    fi
    
    # Get peptide length and calculate residue indices
    PEPLEN=$(get_peptide_len "$name")
    TOTAL_RES=$(echo "q" | gmx make_ndx -f md.tpr 2>&1 | grep "Protein residues" | awk '{print $3}')
    REC_END=$((TOTAL_RES - PEPLEN))
    PEP_START=$((REC_END + 1))
    
    echo "Receptor: residues 1-$REC_END, Peptide: residues $PEP_START-$TOTAL_RES ($PEPLEN residues)"
    
    # Create index file with proper residue selection
    gmx make_ndx -f md.tpr -o analysis.ndx << ENDNDX 2>/dev/null
ri 1-$REC_END
name 17 Receptor
ri $PEP_START-$TOTAL_RES
name 18 Peptide
q
ENDNDX

    # 1. RMSD of backbone
    echo ""
    echo "--- Backbone RMSD ---"
    echo "4 4" | gmx rms -s md.tpr -f $XTC -o rmsd_backbone.xvg -tu ns 2>/dev/null
    INIT=$(head -30 rmsd_backbone.xvg | grep -v "^[#@]" | head -1 | awk '{printf "%.3f", $2}')
    FINAL=$(tail -5 rmsd_backbone.xvg | grep -v "^[#@]" | tail -1 | awk '{printf "%.3f", $2}')
    MAX=$(grep -v "^[#@]" rmsd_backbone.xvg | awk 'BEGIN{max=0} {if($2>max) max=$2} END{printf "%.3f", max}')
    echo "Initial: $INIT nm | Final: $FINAL nm | Max: $MAX nm"
    
    # 2. RMSD of peptide backbone only
    echo ""
    echo "--- Peptide RMSD ---"
    echo "18 18" | gmx rms -s md.tpr -f $XTC -n analysis.ndx -o rmsd_peptide.xvg -tu ns 2>/dev/null
    INIT=$(head -30 rmsd_peptide.xvg | grep -v "^[#@]" | head -1 | awk '{printf "%.3f", $2}')
    FINAL=$(tail -5 rmsd_peptide.xvg | grep -v "^[#@]" | tail -1 | awk '{printf "%.3f", $2}')
    MAX=$(grep -v "^[#@]" rmsd_peptide.xvg | awk 'BEGIN{max=0} {if($2>max) max=$2} END{printf "%.3f", max}')
    echo "Initial: $INIT nm | Final: $FINAL nm | Max: $MAX nm"
    
    # 3. Minimum distance peptide-receptor
    echo ""
    echo "--- Minimum Distance (Peptide-Receptor) ---"
    echo "18 17" | gmx mindist -s md.tpr -f $XTC -n analysis.ndx -od mindist.xvg -tu ns 2>/dev/null
    if [ -f mindist.xvg ]; then
        INIT=$(head -30 mindist.xvg | grep -v "^[#@]" | head -1 | awk '{printf "%.3f", $2}')
        FINAL=$(tail -5 mindist.xvg | grep -v "^[#@]" | tail -1 | awk '{printf "%.3f", $2}')
        AVG=$(grep -v "^[#@]" mindist.xvg | awk '{sum+=$2; n++} END {printf "%.3f", sum/n}')
        MIN=$(grep -v "^[#@]" mindist.xvg | awk 'BEGIN{min=999} {if($2<min) min=$2} END{printf "%.3f", min}')
        echo "Initial: $INIT nm | Final: $FINAL nm | Avg: $AVG nm | Min: $MIN nm"
    else
        echo "Failed to calculate"
    fi
    
    # 4. H-bonds
    echo ""
    echo "--- Hydrogen Bonds (Peptide-Receptor) ---"
    echo "18 17" | gmx hbond -s md.tpr -f $XTC -n analysis.ndx -num hbonds.xvg 2>/dev/null
    if [ -f hbonds.xvg ]; then
        INIT=$(head -30 hbonds.xvg | grep -v "^[#@]" | head -1 | awk '{print $2}')
        FINAL=$(tail -5 hbonds.xvg | grep -v "^[#@]" | tail -1 | awk '{print $2}')
        AVG=$(grep -v "^[#@]" hbonds.xvg | awk '{sum+=$2; n++} END {printf "%.1f", sum/n}')
        MAX=$(grep -v "^[#@]" hbonds.xvg | awk 'BEGIN{max=0} {if($2>max) max=$2} END{print max}')
        echo "Initial: $INIT | Final: $FINAL | Avg: $AVG | Max: $MAX"
    else
        echo "Failed to calculate"
    fi
    
    cd ..
done
