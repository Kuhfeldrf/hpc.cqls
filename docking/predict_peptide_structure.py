#!/usr/bin/env python3
"""
Peptide Structure Prediction using ESMFold/AlphaFold API
========================================================
Predicts 3D structure from amino acid sequence using ESMFold API
(fast alternative to AlphaFold) or local ColabFold if available.

Usage:
    python predict_peptide_structure.py peptide_prediction_list.txt
    python predict_peptide_structure.py --sequence YPFPGP
    python predict_peptide_structure.py peptide_prediction_list.txt --output-dir peptides

Input file format (one sequence per line):
    GLAPYKLRPVAA
    LLFKDSAIGF
    YPFPGP
"""

import sys
import os
import argparse
import time
from pathlib import Path

# Check for required packages
try:
    import requests
except ImportError:
    print("Error: 'requests' package not found. Install with: pip install requests")
    sys.exit(1)


def validate_sequence(sequence):
    """
    Validate amino acid sequence.
    
    Returns:
    - (bool, str): (is_valid, error_message)
    """
    sequence = sequence.upper().strip()
    valid_aa = set("ACDEFGHIKLMNPQRSTVWY")
    
    if not sequence:
        return False, "Empty sequence"
    
    invalid = [aa for aa in sequence if aa not in valid_aa]
    if invalid:
        return False, f"Invalid amino acids: {set(invalid)}"
    
    if len(sequence) < 5:
        return False, "Sequence too short (minimum 5 amino acids)"
    
    if len(sequence) > 2000:
        return False, "Sequence too long (maximum 2000 amino acids)"
    
    return True, ""


def predict_esmfold_api(sequence, output_pdb):
    """
    Predict structure using ESMFold API (Meta AI).
    Free, fast alternative to AlphaFold.
    
    Parameters:
    - sequence: Amino acid sequence
    - output_pdb: Path to save output PDB file
    
    Returns:
    - str: Path to predicted PDB file, or None if failed
    """
    api_url = "https://api.esmatlas.com/foldSequence/v1/pdb/"
    
    print(f"   Submitting to ESMFold API...")
    
    try:
        response = requests.post(api_url, data=sequence, timeout=300)
        
        if response.status_code == 200:
            with open(output_pdb, 'w') as f:
                f.write(response.text)
            return str(output_pdb)
        else:
            print(f"   API request failed (status {response.status_code})")
            if response.text:
                print(f"   Response: {response.text[:200]}")
            return None
            
    except requests.exceptions.Timeout:
        print("   API request timed out (>5 minutes)")
        return None
    except requests.exceptions.RequestException as e:
        print(f"   API request error: {e}")
        return None


def predict_colabfold_local(sequence, output_pdb):
    """
    Predict using local ColabFold installation (if available).
    
    Parameters:
    - sequence: Amino acid sequence
    - output_pdb: Path to save output PDB file
    
    Returns:
    - str: Path to predicted PDB file, or None if failed
    """
    import subprocess
    import tempfile
    import shutil
    
    # Check if colabfold_batch is available
    result = subprocess.run(["which", "colabfold_batch"], capture_output=True)
    
    if result.returncode != 0:
        print("   ColabFold not found locally")
        return None
    
    print("   Using local ColabFold...")
    
    # Create temporary FASTA file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as tmp_fasta:
        tmp_fasta.write(f">peptide\n{sequence}\n")
        tmp_fasta_path = tmp_fasta.name
    
    try:
        output_dir = Path(output_pdb).parent
        cmd = ["colabfold_batch", tmp_fasta_path, str(output_dir)]
        
        print(f"   Running: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
        
        if result.returncode == 0:
            # Find output PDB file
            predicted_files = list(output_dir.glob(f"*{sequence}*.pdb"))
            if predicted_files:
                if predicted_files[0] != Path(output_pdb):
                    shutil.copy(predicted_files[0], output_pdb)
                return str(output_pdb)
            else:
                print("   Output PDB file not found")
                return None
        else:
            print(f"   ColabFold failed: {result.stderr[:300]}")
            return None
            
    except subprocess.TimeoutExpired:
        print("   ColabFold timed out (>1 hour)")
        return None
    except Exception as e:
        print(f"   Error: {e}")
        return None
    finally:
        if os.path.exists(tmp_fasta_path):
            os.unlink(tmp_fasta_path)


def predict_peptide_structure(sequence, output_dir="peptides", method="auto"):
    """
    Predict peptide 3D structure from amino acid sequence.
    
    Parameters:
    - sequence: Amino acid sequence (single letter code)
    - output_dir: Directory to save output files
    - method: "auto" (try ESMFold API, fall back to local), "esmfold", or "colabfold"
    
    Returns:
    - str: Path to predicted PDB file, or None if failed
    """
    sequence = sequence.upper().strip()
    
    # Validate sequence
    is_valid, error = validate_sequence(sequence)
    if not is_valid:
        print(f"   ✗ Invalid sequence: {error}")
        return None
    
    # Create output directory
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    output_pdb = output_path / f"{sequence}.pdb"
    
    # Check if already exists
    if output_pdb.exists():
        print(f"   ✓ Structure already exists: {output_pdb}")
        return str(output_pdb)
    
    print(f"   Length: {len(sequence)} amino acids")
    
    # Try prediction methods
    result = None
    
    if method in ["auto", "esmfold"]:
        result = predict_esmfold_api(sequence, output_pdb)
    
    if result is None and method in ["auto", "colabfold"]:
        result = predict_colabfold_local(sequence, output_pdb)
    
    if result:
        print(f"   ✓ Structure saved: {output_pdb}")
    else:
        print(f"   ✗ Prediction failed")
    
    return result


def read_peptide_list(filename):
    """
    Read peptide sequences from a text file.
    
    Parameters:
    - filename: Path to text file (one sequence per line)
    
    Returns:
    - list: List of peptide sequences
    """
    peptides = []
    
    try:
        with open(filename, 'r') as f:
            for line in f:
                line = line.strip()
                # Skip empty lines and comments
                if line and not line.startswith('#') and not line.startswith('>'):
                    # Handle FASTA format (sequence follows header)
                    peptides.append(line)
        
        return peptides
        
    except FileNotFoundError:
        print(f"Error: File not found: {filename}")
        return []
    except Exception as e:
        print(f"Error reading file: {e}")
        return []


def main():
    parser = argparse.ArgumentParser(
        description='Predict peptide 3D structure using ESMFold/AlphaFold API',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s peptide_prediction_list.txt
  %(prog)s --sequence GLAPYKLRPVAA
  %(prog)s peptide_prediction_list.txt --output-dir my_predictions
  %(prog)s peptide_prediction_list.txt --method colabfold
"""
    )
    
    parser.add_argument('input_file', nargs='?', 
                        help='Text file with peptide sequences (one per line)')
    parser.add_argument('--sequence', '-s', type=str,
                        help='Single peptide sequence to predict')
    parser.add_argument('--output-dir', '-o', type=str, default='peptides',
                        help='Output directory for PDB files (default: peptides)')
    parser.add_argument('--method', '-m', type=str, default='auto',
                        choices=['auto', 'esmfold', 'colabfold'],
                        help='Prediction method (default: auto)')
    parser.add_argument('--delay', '-d', type=float, default=2.0,
                        help='Delay between API requests in seconds (default: 2.0)')
    
    args = parser.parse_args()
    
    # Get sequences to process
    sequences = []
    
    if args.sequence:
        sequences = [args.sequence]
    elif args.input_file:
        sequences = read_peptide_list(args.input_file)
        if not sequences:
            print("No valid sequences found in input file")
            sys.exit(1)
    else:
        parser.print_help()
        print("\nError: Provide either an input file or --sequence argument")
        sys.exit(1)
    
    # Print header
    print("=" * 60)
    print("Peptide Structure Prediction")
    print("=" * 60)
    print(f"Sequences to process: {len(sequences)}")
    print(f"Output directory: {args.output_dir}")
    print(f"Method: {args.method}")
    print("=" * 60)
    
    # Process each sequence
    results = {}
    successful = 0
    failed = 0
    
    for i, seq in enumerate(sequences, 1):
        print(f"\n[{i}/{len(sequences)}] {seq}")
        
        result = predict_peptide_structure(seq, args.output_dir, args.method)
        results[seq] = result
        
        if result:
            successful += 1
        else:
            failed += 1
        
        # Add delay between requests to avoid rate limiting
        if i < len(sequences) and args.method in ['auto', 'esmfold']:
            time.sleep(args.delay)
    
    # Print summary
    print("\n" + "=" * 60)
    print("Summary")
    print("=" * 60)
    print(f"Total: {len(sequences)} | Successful: {successful} | Failed: {failed}")
    print()
    
    for seq, pdb_file in results.items():
        status = "✓" if pdb_file else "✗"
        result_str = pdb_file if pdb_file else "Failed"
        print(f"  {status} {seq}: {result_str}")
    
    print("=" * 60)
    
    # Exit with error code if any failed
    sys.exit(0 if failed == 0 else 1)


if __name__ == "__main__":
    main()

