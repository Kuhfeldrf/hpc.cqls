# HPC Computational Biology Workflow

This repository contains computational biology workflows for molecular dynamics simulations and protein-ligand docking.

## Directory Structure

### `docking/`
Contains Rosetta-based protein-ligand docking workflows:
- Ligand structures from AlphaFold predictions
- Docking scripts and protocols
- Post-processing tools for GROMACS compatibility

### `md/`
Contains GROMACS molecular dynamics simulation workflows:
- Structure preparation scripts
- MD simulation scripts
- Analysis tools

## Workflow Overview

1. **Docking**: Protein-ligand docking using Rosetta
   - Input: Receptor (3fxi_0001.pdb) + Ligands (alphafold_predictions/)
   - Output: Docked complexes ready for MD simulation

2. **Structure Preparation**: Prepare structures for MD
   - Convert docking outputs to GROMACS-compatible format
   - Generate topologies and system setup files

3. **MD Simulation**: Run production MD simulations
   - 100ns production runs
   - GPU-accelerated on GTX 1080 Ti

## Quick Start

### Docking
```bash
cd docking
sbatch dock_ligands_to_3fxi.sh
```

### Structure Preparation
```bash
cd md
./prepare_all_structures.sh
```

### MD Simulation
```bash
cd md
./submit_all_recligand.sh
```

## Requirements

- GROMACS 2024.5
- Rosetta (for docking)
- SLURM job scheduler
- GPU access (for MD simulations)

## Notes

- Large data files (trajectories, PDB files) are excluded from git
- Archive directories contain previous simulation results
- Job output files (*.out, *.err) are excluded from version control

