# HPC Computational Biology Workflow

This repository contains computational biology workflows for peptide-protein docking and molecular dynamics simulations. The workflow consists of two main pipelines:

1. **Docking Pipeline** (`docking/submit_docking_pipeline.sh`) - Predicts peptide structures and docks them to a receptor
2. **MD Pipeline** (`md/run_md_pipeline.sh`) - Runs molecular dynamics simulations on docked complexes

## Directory Structure

### `docking/`
Contains the peptide prediction and Rosetta-based docking workflow:
- `peptides/` - Predicted peptide structures (from ESMFold)
- `receptor/` - Receptor structure (single PDB file)
- `docked_complexes/` - Docking output and score files
- `submit_docking_pipeline.sh` - Main docking pipeline script (Slurm-based)
- `predict_peptide_structure.py` - ESMFold API structure prediction
- `relax_peptides_container.sh` - Rosetta relaxation
- `flexpep_dock_parallel.sh` - FlexPepDocking

### `md/`
Contains GROMACS molecular dynamics simulation workflow:
- Structure preparation scripts
- MD simulation scripts
- Analysis and plotting tools
- `run_md_pipeline.sh` - Main MD pipeline script

---

## Pipeline 1: Docking (`docking/`)

The docking pipeline predicts peptide structures, performs **global docking** to find binding poses, then **refines** with FlexPepDocking:

```
┌──────────────────────────────────────────────────────────────────┐
│  submit_docking_pipeline.sh <peptides.txt> [nrefine] [nglobal]   │
└──────────────────────────────────────────────────────────────────┘
                              │
        ┌─────────────────────┼─────────────────────┐
        │                     │                     │
        ▼                     ▼                     ▼
┌───────────────┐     ┌───────────────┐     ┌───────────────┐
│  GLAPYKLRPVAA │     │  LLFKDSAIGF   │     │    ...        │
└───────┬───────┘     └───────┬───────┘     └───────┬───────┘
        │                     │                     │
        ▼                     ▼                     ▼
┌─────────────────────────────────────────────────────────────────┐
│  STEP 1: Predict Structures (ESMFold API)                       │
│  Creates: peptides/GLAPYKLRPVAA.pdb, peptides/LLFKDSAIGF.pdb    │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│  STEP 2: Global Docking (Rosetta docking_protocol)              │
│  - Searches ENTIRE receptor surface                             │
│  - Randomizes peptide position                                  │
│  - Generates ~50 candidate poses per peptide                    │
│  Creates: global_docked/*_best_global.pdb                       │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│  STEP 3: FlexPepDocking Refinement (Rosetta)                    │
│  - Full peptide backbone flexibility                            │
│  - Side chain optimization                                      │
│  - Refines best global poses                                    │
│  Creates: docked_complexes/*_docked.pdb, *_scores.sc            │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│  STEP 4: Score Compilation & Best Structure Selection           │
│  Creates: docking_summary.csv                                   │
│  Copies: Best structure per peptide → md/ folder                │
└─────────────────────────────────────────────────────────────────┘
```

### Docking Pipeline Usage

```bash
cd docking

# Run complete docking pipeline (with global docking)
./submit_docking_pipeline.sh peptide_list.txt 10 50

# Arguments:
#   peptide_list.txt - Text file with peptide sequences (one per line)
#   10               - FlexPepDocking refinement structures per peptide (default: 10)
#   50               - Global docking structures per peptide (default: 50)

# Example peptide_list.txt:
#   GLAPYKLRPVAA
#   LLFKDSAIGF
#   YPFPGP
```

### Docking Pipeline Output
- `peptides/*.pdb` - Predicted peptide structures
- `global_docked/*_best_global.pdb` - Best global docking poses
- `global_docked/*_global_scores.sc` - Global docking scores
- `docked_complexes/*_docked.pdb` - Refined docked structures
- `docked_complexes/docking_summary_*.csv` - Compiled scores
- `../md/<peptide>_<receptor>.pdb` - Best structure for each peptide (ready for MD)

---

## Pipeline 2: Molecular Dynamics (`md/`)

The MD pipeline runs GROMACS simulations on docked complexes. It is orchestrated by `run_md_pipeline.sh` which uses SLURM job dependencies to chain all steps:

```
┌─────────────────────────────────────────────────────────────────┐
│            run_md_pipeline.sh STRUCT1 STRUCT2 STRUCT3 ...       │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│  STEP 1: PREPARATION (Array Job)                                │
│  ┌─────────┐ ┌─────────┐ ┌─────────┐ ┌─────────┐ ┌─────────┐   │
│  │ Task 0  │ │ Task 1  │ │ Task 2  │ │ Task 3  │ │ Task 4  │   │
│  │LKQVLLH..│ │GLAPYKL..│ │RVVRDPQ..│ │RPKLPLR..│ │RWKIFKK..│   │
│  └────┬────┘ └────┬────┘ └────┬────┘ └────┬────┘ └────┬────┘   │
└───────│──────────│──────────│──────────│──────────│──────────┘
        │          │          │          │          │
        ▼          ▼          ▼          ▼          ▼   (aftercorr)
┌─────────────────────────────────────────────────────────────────┐
│  STEP 2: PRODUCTION MD (5 Independent Jobs, ~12 hours each)     │
│  ┌─────────┐ ┌─────────┐ ┌─────────┐ ┌─────────┐ ┌─────────┐   │
│  │ MD Job  │ │ MD Job  │ │ MD Job  │ │ MD Job  │ │ MD Job  │   │
│  │  100ns  │ │  100ns  │ │  100ns  │ │  100ns  │ │  100ns  │   │
│  └────┬────┘ └────┬────┘ └────┬────┘ └────┬────┘ └────┬────┘   │
└───────│──────────│──────────│──────────│──────────│──────────┘
        │          │          │          │          │
        ▼          ▼          ▼          ▼          ▼   (afterok)
┌─────────────────────────────────────────────────────────────────┐
│  STEP 3: ANALYSIS (5 Independent Jobs)                          │
│  ┌─────────┐ ┌─────────┐ ┌─────────┐ ┌─────────┐ ┌─────────┐   │
│  │Analysis │ │Analysis │ │Analysis │ │Analysis │ │Analysis │   │
│  └────┬────┘ └────┬────┘ └────┬────┘ └────┬────┘ └────┬────┘   │
└───────│──────────│──────────│──────────│──────────│──────────┘
        └──────────┴──────────┴──────────┴──────────┘
                              │
                              ▼   (afterok: ALL analysis jobs)
┌─────────────────────────────────────────────────────────────────┐
│  STEP 4: PLOTTING (Single Job)                                  │
│  ┌───────────────────────────────────────────────────────────┐  │
│  │  Generate plots for all 5 simulations                     │  │
│  └───────────────────────────────────────────────────────────┘  │
└─────────────────────────────────────────────────────────────────┘
```

### MD Pipeline Usage

```bash
cd md

# Run the complete MD pipeline with structures as arguments
./run_md_pipeline.sh AVAVVKKGGSFQL_3fxi GLAPYKLRPVAA_3fxi LKKISQRYQKF_3fxi

# Test first (shows what would be submitted without actually submitting)
./run_md_pipeline.sh --dry-run AVAVVKKGGSFQL_3fxi GLAPYKLRPVAA_3fxi

# Skip preparation if structures are already prepared
./run_md_pipeline.sh --skip-prep AVAVVKKGGSFQL_3fxi GLAPYKLRPVAA_3fxi

# Skip both prep and MD (just run analysis and plotting)
./run_md_pipeline.sh --skip-prep --skip-md AVAVVKKGGSFQL_3fxi GLAPYKLRPVAA_3fxi

# Note: Structure names can be provided with or without the receptor suffix
# The script will automatically detect and use the correct PDB file
```

---

## Complete Workflow: Docking → MD

To run both pipelines sequentially:

```bash
# Step 1: Docking Pipeline
cd /home/kuhfeldr/hpc.cqls/docking
./submit_docking_pipeline.sh peptide_list.txt 10 50
# Monitor with: squeue -u $USER
# Best structures are automatically copied to md/ folder when complete

# Step 2: MD Pipeline (use structures from docking)
cd /home/kuhfeldr/hpc.cqls/md
./run_md_pipeline.sh GLAPYKLRPVAA_3fxi LLFKDSAIGF_3fxi
```

---

## Individual Steps (Manual Execution)

### Docking Steps
```bash
cd docking

# 1. Predict peptide structures (ESMFold API)
python3 predict_peptide_structure.py peptide_list.txt --output-dir peptides

# 2. Global docking (searches entire receptor surface)
sbatch global_dock_parallel.sh

# 3. FlexPepDocking refinement (refines best global poses)
sbatch flexpep_dock_parallel.sh
```

### MD Steps

#### 1. Structure Preparation
```bash
cd md
sbatch prepare_structures_parallel.sh
```

**`prepare_structures_parallel.sh`** - Prepares protein-peptide structures for MD simulation:
- **Purpose**: Parallel structure preparation using SLURM array jobs
- **Input**: PDB files (e.g., `STRUCTURE_NAME_3fxi.pdb`) and `.mdp` parameter files
- **Output**: Prepared structures in individual directories with `md.tpr` files
- **Steps Performed**:
  1. Topology generation (`pdb2gmx`) with chain separation (`-chainsep id_or_ter -merge no`)
  2. Simulation box creation (cubic, 1.2 nm buffer)
  3. Solvation with TIP3P water
  4. Ion addition (0.15 M NaCl)
  5. Energy minimization
  6. NVT equilibration (250 ps at 300 K)
  7. NPT equilibration (500 ps at 300 K, 1 bar)
  8. Production MD TPR generation
  9. Analysis index creation (`analysis.ndx` with Receptor/Peptide groups)
- **Usage**: Typically run via `run_full_pipeline.sh` as an array job, but can be run standalone
- **Requirements**: GROMACS container, GPU for equilibration steps

#### 2. MD Simulation
```bash
cd md/<structure_name>
sbatch ../submit_gromacs_gpu.sh
```

**`submit_gromacs_gpu.sh`** - Production MD simulation on GPU:
- **Purpose**: Runs 100 ns production MD simulation using GPU acceleration
- **Input**: Prepared structure directory with `md.tpr` file (from preparation step)
- **Output**: Trajectory (`md.xtc`), energy file (`md.edr`), log (`md.log`), final structure (`md.gro`)
- **Key Features**:
  - Multi-GPU support (4x L40S GPUs recommended)
  - Automatic checkpoint detection and continuation
  - Full GPU acceleration (NB, PME, bonded, update)
  - Checkpoint every 60 minutes for fault tolerance
  - Expected performance: ~50-80 ns/day with 4x L40S GPUs
- **Settings**: 300 K, 1 bar, 100 ns duration, includes all atoms in trajectory
- **Requirements**: GPU node (L40S or A30), GROMACS container
- **Note**: For 100 ns runs on normal partition (1-day limit), job will checkpoint and can be continued

#### 3. Analysis
```bash
cd md
./analyze_md.sh /path/to/simulation_dir
```

**`analyze_md.sh`** - Comprehensive MD trajectory analysis:
- **Purpose**: Analyzes MD trajectories to extract key metrics and generate XVG data files
- **Input**: Simulation directory containing `md.tpr` and trajectory file (`md.xtc`)
- **Output**: XVG files in `analysis_output/analysis_data/` directory
- **Analyses Performed**:
  1. **Receptor Backbone RMSD**: Receptor stability over time
  2. **Peptide RMSD (binding)**: Peptide displacement from initial docked pose (fit to receptor)
  3. **Peptide Backbone RMSD**: Peptide conformational flexibility (fit to itself)
  4. **RMSF (per residue)**: Flexibility analysis for receptor and peptide separately
  5. **Radius of Gyration**: Compactness measure for receptor and peptide
  6. **Hydrogen Bonds**: Intermolecular H-bonds between peptide and receptor
- **Key Features**:
  - Automatic receptor/peptide separation using custom index groups
  - PBC (Periodic Boundary Condition) correction for accurate analysis
  - Supports both moltype-based and residue-range selection methods
  - Filters outlier residues in peptide RMSF analysis
  - Works with GROMACS container (no GPU required for analysis)
- **Output Files**: All XVG files are saved in `analysis_output/analysis_data/<structure_name>/` and copied to `analysis_output/output_files/`
- **Requirements**: GROMACS container, trajectory file from production MD

#### 4. Plotting
```bash
cd md
# Plot a single simulation
module load intel-python/24.0.0
python3 plot_simulations.py /path/to/simulation_dir --peptide PEPTIDE_NAME

# Or use the plotting script (generates plots for all structures)
sbatch submit_plot_all.sh STRUCT1 STRUCT2 STRUCT3
```

**`plot_simulations.py`** - Generates publication-quality plots from GROMACS MD analysis data:
- **Input**: Simulation directory containing `analysis_output/` with XVG files from `analyze_md.sh`
- **Output**: PNG plots saved to `analysis_output/plots/`
- **Plots Generated**:
  - **RMSD Analysis**: Receptor backbone stability, peptide binding stability (vs receptor), and peptide internal conformational flexibility
  - **RMSF Analysis**: Per-residue root mean square fluctuation for receptor and peptide
  - **Radius of Gyration**: Time evolution and distribution for receptor and peptide
  - **Hydrogen Bonds**: Stability over time and distribution histogram
  - **Combined Summary**: 3x3 grid with all key metrics
- **Options**:
  - `--peptide NAME`: Specify peptide name for plot titles (defaults to directory name)
  - `--summary-only`: Generate only the combined summary plot
- **Requirements**: Python 3 with `numpy` and `matplotlib` (use `intel-python/24.0.0` module on ORCA)

## Estimated Timeline (100ns run using L40s GPUs on ORCA)

### Docking Pipeline
| Step | Duration | Notes |
|------|----------|-------|
| Structure Prediction | ~1 min | Per peptide (ESMFold API) |
| Global Docking | ~2-4 hours | Per peptide (searches entire surface) |
| FlexPepDocking Refinement | ~1-2 hours | Per peptide (refines best poses) |
| Score Compilation | ~1 min | Single job |

### MD Pipeline
| Step | Duration | Notes |
|------|----------|-------|
| Preparation | <1 hour | All structures in parallel |
| Production MD | ~12 hours | All structures in parallel (100ns each) |
| Analysis | < 1 hour | All structures in parallel |
| Plotting | < 1 min | Single job |

## Requirements

- GROMACS 2025.3 (via Apptainer container)
- Rosetta (via Apptainer container)
- SLURM job scheduler
- GPU access (NVIDIA A30 or L40s similar for MD simulations)
- Python 3 with matplotlib, numpy (for plotting)

## GROMACS Input Files (.mdp)

The MD simulation workflow uses several GROMACS parameter files (`.mdp` files) that define simulation parameters for each stage:

### `minim.mdp` - Energy Minimization
- **Purpose**: Initial energy minimization to remove steric clashes and bad contacts
- **Algorithm**: Steepest descent minimization
- **Key Parameters**:
  - Maximum force tolerance: 1000.0 kJ/mol/nm
  - Step size: 0.01 nm
  - Maximum steps: 50,000
- **Used in**: Step 5 of structure preparation (after ion addition)

### `ions.mdp` - Ion Addition Preprocessing
- **Purpose**: Energy minimization parameters for preprocessing before adding ions
- **Algorithm**: Steepest descent minimization (same as minim.mdp)
- **Key Parameters**: Same as minim.mdp
- **Used in**: Step 4 of structure preparation (before genion adds Na+/Cl- ions)

### `nvt.mdp` - NVT Equilibration (Constant Volume, Temperature)
- **Purpose**: Temperature equilibration at constant volume
- **Duration**: 250 ps (increased for docked protein-peptide complexes)
- **Key Parameters**:
  - Temperature: 300 K (V-rescale thermostat)
  - Position restraints: Applied to protein (`-DPOSRES`)
  - Velocity generation: Yes (Maxwell distribution at 300 K)
  - No pressure coupling
- **Used in**: Step 6 of structure preparation (after energy minimization)

### `npt.mdp` - NPT Equilibration (Constant Pressure, Temperature)
- **Purpose**: Pressure and density equilibration at constant temperature
- **Duration**: 500 ps (increased for docked protein-peptide complexes)
- **Key Parameters**:
  - Temperature: 300 K (V-rescale thermostat)
  - Pressure: 1.0 bar (Parrinello-Rahman barostat)
  - Position restraints: Applied to protein (`-DPOSRES`)
  - Continuation: Yes (from NVT)
- **Used in**: Step 7 of structure preparation (after NVT equilibration)

### `md.mdp` - Production MD Simulation
- **Purpose**: Production molecular dynamics simulation
- **Duration**: 100 ns
- **Key Parameters**:
  - Temperature: 300 K (V-rescale thermostat)
  - Pressure: 1.0 bar (Parrinello-Rahman barostat)
  - Trajectory output: Every 10 ps (includes all atoms including water)
  - Energy output: Every 10 ps
  - No position restraints (system is fully equilibrated)
  - Continuation: Yes (from NPT)
- **Used in**: Step 8 (production MD run)

**Simulation Sequence**: `minim.mdp` → `nvt.mdp` → `npt.mdp` → `md.mdp`

## Containers

- **GROMACS**: `/home/kuhfeldr/gromacs_v2025.3.sif`
- **Rosetta**: `/home/kuhfeldr/rosetta.sif`

## Notes

- Large data files (trajectories, PDB files) are excluded from git
- Archive directories contain previous simulation results
- Job output files (*.out, *.err) are excluded from version control
- Monitor jobs with: `squeue -u $USER`



## Gnerating Plots on ORCA
    module load intel-python/24.0.0
    python3 /home/kuhfeldr/hpc.cqls/md/plot_simulations.py /scratch/kuhfeldr-rktemp/GLAPYKLRPVAA_3fxi_100ns/