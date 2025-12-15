#!/usr/bin/env python3
"""
Generate publication-quality plots from GROMACS MD simulation data.
Core analyses: RMSD (backbone and peptide), RMSF, Radius of Gyration, and Hydrogen bonds

For extended analyses (DSSP, MinDist), use plot_simulations_extended.py
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import argparse
import sys


def parse_xvg(filename):
    """
    Parse GROMACS XVG file format.
    
    Returns:
    - time: Array of time values (in nanoseconds)
    - data: Array of data values (can be multi-column)
    - labels: List of y-axis labels if found
    """
    time = []
    data = []
    time_unit = None
    y_labels = []
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            
            if not line or line.startswith('#'):
                continue
            
            if line.startswith('@'):
                if 'xaxis' in line and 'label' in line:
                    if 'ps' in line.lower() or 'picosecond' in line.lower():
                        time_unit = 'ps'
                    elif 'ns' in line.lower() or 'nanosecond' in line.lower():
                        time_unit = 'ns'
                if line.startswith('@ s') and 'legend' in line:
                    parts = line.split('"')
                    if len(parts) >= 2:
                        y_labels.append(parts[1])
                continue
            
            try:
                parts = line.split()
                if len(parts) >= 2:
                    t = float(parts[0])
                    vals = [float(p) for p in parts[1:]]
                    time.append(t)
                    data.append(vals)
            except (ValueError, IndexError):
                continue
    
    time = np.array(time)
    data = np.array(data)
    
    if time_unit is None and len(time) > 0:
        max_time = np.max(time)
        time_unit = 'ps' if max_time > 1000 else 'ns'
    
    if time_unit == 'ps':
        time = time / 1000.0
    
    return time, data, y_labels


def find_xvg_file(base_dir, filename):
    """Find an XVG file by checking known locations from analyze_md.sh output."""
    base_dir = Path(base_dir)
    dir_name = base_dir.name
    
    search_paths = [
        base_dir / "analysis_output" / "analysis_data" / dir_name / filename,
        base_dir / "analysis_output" / "output_files" / f"{dir_name}_{filename}",
        base_dir / "analysis_output" / "analysis_data" / filename,
        base_dir / "analysis_output" / "output_files" / filename,
        base_dir / filename,
        base_dir / "analysis_data" / filename,
        base_dir / "analysis_data" / dir_name / filename,
        base_dir / "output_files" / filename,
        base_dir / "output_files" / f"{dir_name}_{filename}",
    ]
    
    for path in search_paths:
        if path.exists():
            return path
    
    return None


def get_output_dir(base_dir):
    """Determine appropriate output directory for plots."""
    base_dir = Path(base_dir)
    analysis_output = base_dir / "analysis_output"
    if analysis_output.exists():
        plots_dir = analysis_output / "plots"
        plots_dir.mkdir(exist_ok=True)
        return plots_dir
    return base_dir


def plot_hbond_analysis(time, hbonds, peptide_name, output_dir):
    """Create H-bond stability and distribution plots."""
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))
    
    mean_hbonds = np.mean(hbonds)
    std_hbonds = np.std(hbonds)
    min_hbonds = np.min(hbonds)
    max_hbonds = np.max(hbonds)
    
    # Plot 1: H-bond stability over time
    ax1.plot(time, hbonds, color='#9B59B6', linewidth=0.8, alpha=0.8)
    ax1.axhline(y=mean_hbonds, color='red', linestyle='--', linewidth=2, label=f'Mean: {mean_hbonds:.1f}')
    ax1.fill_between(time, mean_hbonds - std_hbonds, mean_hbonds + std_hbonds, 
                      color='#F39C12', alpha=0.3, label='±1 SD')
    
    stats_text = f'Mean: {mean_hbonds:.1f}\nStd: {std_hbonds:.1f}\nMin: {int(min_hbonds)}\nMax: {int(max_hbonds)}'
    ax1.text(0.02, 0.98, stats_text, transform=ax1.transAxes, 
             fontsize=10, verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    ax1.set_xlabel('Time (ns)', fontsize=12, fontweight='bold')
    ax1.set_ylabel('H-bonds', fontsize=12, fontweight='bold')
    ax1.set_title(f'{peptide_name} - Hydrogen Bond Stability Over Time', 
                  fontsize=14, fontweight='bold')
    ax1.legend(loc='upper right', fontsize=10)
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: H-bond distribution
    bins = np.arange(int(min_hbonds) - 0.5, int(max_hbonds) + 1.5, 1)
    ax2.hist(hbonds, bins=bins, color='#9B59B6', edgecolor='black', alpha=0.8)
    ax2.axvline(x=mean_hbonds, color='red', linestyle='--', linewidth=2, 
                label=f'Mean: {mean_hbonds:.1f}')
    
    ax2.set_xlabel('Number of Hydrogen Bonds', fontsize=12, fontweight='bold')
    ax2.set_ylabel('Frequency', fontsize=12, fontweight='bold')
    ax2.set_title(f'{peptide_name} - Hydrogen Bond Distribution', 
                  fontsize=14, fontweight='bold')
    ax2.legend(loc='upper right', fontsize=10)
    ax2.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    output_file = output_dir / f'{peptide_name}_hbond_analysis.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"✓ Saved: {output_file}")
    plt.close()


def plot_rmsd_analysis(time_backbone, rmsd_backbone, time_peptide, rmsd_peptide, 
                       time_peptide_int, rmsd_peptide_int, peptide_name, output_dir):
    """Create RMSD plots for backbone, peptide binding, and peptide internal."""
    fig, axes = plt.subplots(3, 1, figsize=(12, 14))
    ax1, ax2, ax3 = axes
    
    # Plot 1: Backbone RMSD
    mean_bb = np.mean(rmsd_backbone)
    std_bb = np.std(rmsd_backbone)
    min_bb = np.min(rmsd_backbone)
    max_bb = np.max(rmsd_backbone)
    
    ax1.plot(time_backbone, rmsd_backbone, color='#3498DB', linewidth=1.2, 
             label='Backbone RMSD')
    ax1.axhline(y=mean_bb, color='red', linestyle='--', linewidth=2, 
                label=f'Mean: {mean_bb:.3f} nm')
    ax1.fill_between(time_backbone, mean_bb - std_bb, mean_bb + std_bb, 
                      color='#F39C12', alpha=0.3, label='±1 SD')
    
    stats_text = f'Mean: {mean_bb:.3f} nm\nStd: {std_bb:.3f} nm\nMin: {min_bb:.3f} nm\nMax: {max_bb:.3f} nm'
    ax1.text(0.02, 0.98, stats_text, transform=ax1.transAxes, 
             fontsize=10, verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    ax1.set_xlabel('Time (ns)', fontsize=12, fontweight='bold')
    ax1.set_ylabel('RMSD (nm)', fontsize=12, fontweight='bold')
    ax1.set_title(f'{peptide_name} - Backbone RMSD (Receptor Stability)', 
                  fontsize=14, fontweight='bold')
    ax1.legend(loc='lower right', fontsize=10)
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Peptide RMSD (relative to receptor - binding stability)
    if rmsd_peptide is not None and time_peptide is not None:
        mean_pep = np.mean(rmsd_peptide)
        std_pep = np.std(rmsd_peptide)
        min_pep = np.min(rmsd_peptide)
        max_pep = np.max(rmsd_peptide)
        
        ax2.plot(time_peptide, rmsd_peptide, color='#2ECC71', linewidth=1.2, 
                 label='Peptide RMSD (vs receptor)')
        ax2.axhline(y=mean_pep, color='red', linestyle='--', linewidth=2, 
                    label=f'Mean: {mean_pep:.3f} nm')
        ax2.fill_between(time_peptide, mean_pep - std_pep, mean_pep + std_pep, 
                          color='#F39C12', alpha=0.3, label='±1 SD')
        
        stats_text = f'Mean: {mean_pep:.3f} nm\nStd: {std_pep:.3f} nm\nMin: {min_pep:.3f} nm\nMax: {max_pep:.3f} nm'
        ax2.text(0.02, 0.98, stats_text, transform=ax2.transAxes, 
                 fontsize=10, verticalalignment='top',
                 bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        ax2.set_xlabel('Time (ns)', fontsize=12, fontweight='bold')
        ax2.set_ylabel('RMSD (nm)', fontsize=12, fontweight='bold')
        ax2.set_title(f'{peptide_name} - Peptide RMSD (Binding Stability - fit to receptor)', 
                      fontsize=14, fontweight='bold')
        ax2.legend(loc='upper right', fontsize=10)
        ax2.grid(True, alpha=0.3)
    else:
        ax2.text(0.5, 0.5, 'Peptide binding RMSD data not available', 
                 transform=ax2.transAxes, ha='center', va='center', fontsize=14)
        ax2.set_axis_off()
    
    # Plot 3: Peptide Internal RMSD (conformational flexibility)
    if rmsd_peptide_int is not None and time_peptide_int is not None:
        mean_int = np.mean(rmsd_peptide_int)
        std_int = np.std(rmsd_peptide_int)
        min_int = np.min(rmsd_peptide_int)
        max_int = np.max(rmsd_peptide_int)
        
        ax3.plot(time_peptide_int, rmsd_peptide_int, color='#E67E22', linewidth=1.2, 
                 label='Peptide Internal RMSD')
        ax3.axhline(y=mean_int, color='red', linestyle='--', linewidth=2, 
                    label=f'Mean: {mean_int:.3f} nm')
        ax3.fill_between(time_peptide_int, mean_int - std_int, mean_int + std_int, 
                          color='#9B59B6', alpha=0.3, label='±1 SD')
        
        stats_text = f'Mean: {mean_int:.3f} nm\nStd: {std_int:.3f} nm\nMin: {min_int:.3f} nm\nMax: {max_int:.3f} nm'
        ax3.text(0.02, 0.98, stats_text, transform=ax3.transAxes, 
                 fontsize=10, verticalalignment='top',
                 bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        ax3.set_xlabel('Time (ns)', fontsize=12, fontweight='bold')
        ax3.set_ylabel('RMSD (nm)', fontsize=12, fontweight='bold')
        ax3.set_title(f'{peptide_name} - Peptide Internal RMSD (Conformational Flexibility - fit to self)', 
                      fontsize=14, fontweight='bold')
        ax3.legend(loc='upper right', fontsize=10)
        ax3.grid(True, alpha=0.3)
    else:
        ax3.text(0.5, 0.5, 'Peptide internal RMSD data not available', 
                 transform=ax3.transAxes, ha='center', va='center', fontsize=14)
        ax3.set_axis_off()
    
    plt.tight_layout()
    output_file = output_dir / f'{peptide_name}_rmsd_analysis.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"✓ Saved: {output_file}")
    plt.close()


def plot_rmsf_analysis(residues_receptor, rmsf_receptor, residues_peptide, rmsf_peptide, 
                       peptide_name, output_dir):
    """Create RMSF (flexibility) plots."""
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 10))
    
    if rmsf_receptor is not None:
        mean_rec = np.mean(rmsf_receptor)
        max_rec = np.max(rmsf_receptor)
        max_res = residues_receptor[np.argmax(rmsf_receptor)]
        
        ax1.bar(residues_receptor, rmsf_receptor, color='#3498DB', alpha=0.8, width=1.0)
        ax1.axhline(y=mean_rec, color='red', linestyle='--', linewidth=2, 
                    label=f'Mean: {mean_rec:.3f} nm')
        
        stats_text = f'Mean: {mean_rec:.3f} nm\nMax: {max_rec:.3f} nm (res {int(max_res)})'
        ax1.text(0.02, 0.98, stats_text, transform=ax1.transAxes, 
                 fontsize=10, verticalalignment='top',
                 bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        ax1.set_xlabel('Residue Number', fontsize=12, fontweight='bold')
        ax1.set_ylabel('RMSF (nm)', fontsize=12, fontweight='bold')
        ax1.set_title(f'{peptide_name} - Receptor RMSF per Residue', 
                      fontsize=14, fontweight='bold')
        ax1.legend(loc='upper right', fontsize=10)
        ax1.grid(True, alpha=0.3, axis='y')
    else:
        ax1.text(0.5, 0.5, 'Receptor RMSF data not available', 
                 transform=ax1.transAxes, ha='center', va='center', fontsize=14)
        ax1.set_axis_off()
    
    if rmsf_peptide is not None:
        mean_pep = np.mean(rmsf_peptide)
        max_pep = np.max(rmsf_peptide)
        max_res = residues_peptide[np.argmax(rmsf_peptide)]
        
        ax2.bar(residues_peptide, rmsf_peptide, color='#2ECC71', alpha=0.8, width=1.0)
        ax2.axhline(y=mean_pep, color='red', linestyle='--', linewidth=2, 
                    label=f'Mean: {mean_pep:.3f} nm')
        
        stats_text = f'Mean: {mean_pep:.3f} nm\nMax: {max_pep:.3f} nm (res {int(max_res)})'
        ax2.text(0.02, 0.98, stats_text, transform=ax2.transAxes, 
                 fontsize=10, verticalalignment='top',
                 bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        ax2.set_xlabel('Residue Number', fontsize=12, fontweight='bold')
        ax2.set_ylabel('RMSF (nm)', fontsize=12, fontweight='bold')
        ax2.set_title(f'{peptide_name} - Peptide RMSF per Residue', 
                      fontsize=14, fontweight='bold')
        ax2.legend(loc='upper right', fontsize=10)
        ax2.grid(True, alpha=0.3, axis='y')
    else:
        ax2.text(0.5, 0.5, 'Peptide RMSF data not available', 
                 transform=ax2.transAxes, ha='center', va='center', fontsize=14)
        ax2.set_axis_off()
    
    plt.tight_layout()
    output_file = output_dir / f'{peptide_name}_rmsf_analysis.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"✓ Saved: {output_file}")
    plt.close()


def plot_rog_analysis(time, rog, plot_title, output_dir):
    """Create radius of gyration plot."""
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))
    
    mean_rog = np.mean(rog)
    std_rog = np.std(rog)
    min_rog = np.min(rog)
    max_rog = np.max(rog)
    
    # Determine color based on title
    if 'Peptide' in plot_title:
        color = '#2ECC71'  # Green for peptide
    else:
        color = '#3498DB'  # Blue for receptor
    
    ax1.plot(time, rog, color=color, linewidth=1.0, alpha=0.8)
    ax1.axhline(y=mean_rog, color='red', linestyle='--', linewidth=2, 
                label=f'Mean: {mean_rog:.3f} nm')
    ax1.fill_between(time, mean_rog - std_rog, mean_rog + std_rog, 
                      color='#F39C12', alpha=0.3, label='±1 SD')
    
    stats_text = f'Mean: {mean_rog:.3f} nm\nStd: {std_rog:.3f} nm\nMin: {min_rog:.3f} nm\nMax: {max_rog:.3f} nm'
    ax1.text(0.02, 0.98, stats_text, transform=ax1.transAxes, 
             fontsize=10, verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    ax1.set_xlabel('Time (ns)', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Radius of Gyration (nm)', fontsize=12, fontweight='bold')
    ax1.set_title(f'{plot_title} - Radius of Gyration Over Time', 
                  fontsize=14, fontweight='bold')
    ax1.legend(loc='upper right', fontsize=10)
    ax1.grid(True, alpha=0.3)
    
    ax2.hist(rog, bins=50, color=color, edgecolor='black', alpha=0.8)
    ax2.axvline(x=mean_rog, color='red', linestyle='--', linewidth=2, 
                label=f'Mean: {mean_rog:.3f} nm')
    
    ax2.set_xlabel('Radius of Gyration (nm)', fontsize=12, fontweight='bold')
    ax2.set_ylabel('Frequency', fontsize=12, fontweight='bold')
    ax2.set_title(f'{plot_title} - Radius of Gyration Distribution', 
                  fontsize=14, fontweight='bold')
    ax2.legend(loc='upper right', fontsize=10)
    ax2.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    # Create safe filename from plot_title
    safe_title = plot_title.replace('(', '').replace(')', '').replace(' ', '_')
    output_file = output_dir / f'{safe_title}_rog_analysis.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"✓ Saved: {output_file}")
    plt.close()


def plot_combined_summary(data_dict, peptide_name, output_dir):
    """Create a combined summary plot with key metrics."""
    fig = plt.figure(figsize=(18, 16))
    
    gs = fig.add_gridspec(3, 3, hspace=0.3, wspace=0.25)
    
    plot_idx = 0
    positions = [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2), (2, 0), (2, 1), (2, 2)]
    
    # Backbone RMSD
    if 'rmsd_backbone' in data_dict:
        time, rmsd = data_dict['rmsd_backbone']
        ax = fig.add_subplot(gs[positions[plot_idx][0], positions[plot_idx][1]])
        ax.plot(time, rmsd, color='#3498DB', linewidth=0.8, alpha=0.8)
        ax.axhline(y=np.mean(rmsd), color='red', linestyle='--', linewidth=1.5)
        ax.set_xlabel('Time (ns)', fontsize=10)
        ax.set_ylabel('RMSD (nm)', fontsize=10)
        ax.set_title('Backbone RMSD\n(Receptor Stability)', fontsize=11, fontweight='bold')
        ax.grid(True, alpha=0.3)
        
        stats_text = f'Mean: {np.mean(rmsd):.3f} nm\nMax: {np.max(rmsd):.3f} nm'
        ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, fontsize=9,
                verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        plot_idx += 1
    
    # Peptide RMSD (binding)
    if 'rmsd_peptide' in data_dict:
        time, rmsd = data_dict['rmsd_peptide']
        ax = fig.add_subplot(gs[positions[plot_idx][0], positions[plot_idx][1]])
        ax.plot(time, rmsd, color='#27AE60', linewidth=0.8, alpha=0.8)
        ax.axhline(y=np.mean(rmsd), color='red', linestyle='--', linewidth=1.5)
        ax.set_xlabel('Time (ns)', fontsize=10)
        ax.set_ylabel('RMSD (nm)', fontsize=10)
        ax.set_title('Peptide RMSD\n(Binding Stability)', fontsize=11, fontweight='bold')
        ax.grid(True, alpha=0.3)
        
        stats_text = f'Mean: {np.mean(rmsd):.3f} nm\nMax: {np.max(rmsd):.3f} nm'
        ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, fontsize=9,
                verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        plot_idx += 1
    
    # Peptide Internal RMSD
    if 'rmsd_peptide_internal' in data_dict:
        time, rmsd = data_dict['rmsd_peptide_internal']
        ax = fig.add_subplot(gs[positions[plot_idx][0], positions[plot_idx][1]])
        ax.plot(time, rmsd, color='#E67E22', linewidth=0.8, alpha=0.8)
        ax.axhline(y=np.mean(rmsd), color='red', linestyle='--', linewidth=1.5)
        ax.set_xlabel('Time (ns)', fontsize=10)
        ax.set_ylabel('RMSD (nm)', fontsize=10)
        ax.set_title('Peptide Internal RMSD\n(Conformational Flex)', fontsize=11, fontweight='bold')
        ax.grid(True, alpha=0.3)
        
        stats_text = f'Mean: {np.mean(rmsd):.3f} nm\nMax: {np.max(rmsd):.3f} nm'
        ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, fontsize=9,
                verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        plot_idx += 1
    
    # Hydrogen Bonds over time
    if 'hbonds' in data_dict:
        time, hbonds = data_dict['hbonds']
        ax = fig.add_subplot(gs[positions[plot_idx][0], positions[plot_idx][1]])
        ax.plot(time, hbonds, color='#9B59B6', linewidth=0.8, alpha=0.8)
        ax.axhline(y=np.mean(hbonds), color='red', linestyle='--', linewidth=1.5)
        ax.set_xlabel('Time (ns)', fontsize=10)
        ax.set_ylabel('H-bonds', fontsize=10)
        ax.set_title('Hydrogen Bonds\n(Over Time)', fontsize=11, fontweight='bold')
        ax.grid(True, alpha=0.3)
        
        stats_text = f'Mean: {np.mean(hbonds):.1f}\nMax: {int(np.max(hbonds))}'
        ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, fontsize=9,
                verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        plot_idx += 1
    
    # H-bond distribution
    if 'hbonds' in data_dict:
        time, hbonds = data_dict['hbonds']
        ax = fig.add_subplot(gs[positions[plot_idx][0], positions[plot_idx][1]])
        bins = np.arange(int(np.min(hbonds)) - 0.5, int(np.max(hbonds)) + 1.5, 1)
        ax.hist(hbonds, bins=bins, color='#9B59B6', edgecolor='black', alpha=0.8)
        ax.axvline(x=np.mean(hbonds), color='red', linestyle='--', linewidth=1.5)
        ax.set_xlabel('Number of H-bonds', fontsize=10)
        ax.set_ylabel('Frequency', fontsize=10)
        ax.set_title('H-bond Distribution', fontsize=11, fontweight='bold')
        ax.grid(True, alpha=0.3, axis='y')
        plot_idx += 1
    
    # Radius of Gyration - Receptor
    if 'rog_receptor' in data_dict and plot_idx < len(positions):
        time, rog = data_dict['rog_receptor']
        ax = fig.add_subplot(gs[positions[plot_idx][0], positions[plot_idx][1]])
        ax.plot(time, rog, color='#3498DB', linewidth=0.8, alpha=0.8)
        ax.axhline(y=np.mean(rog), color='red', linestyle='--', linewidth=1.5)
        ax.set_xlabel('Time (ns)', fontsize=10)
        ax.set_ylabel('RoG (nm)', fontsize=10)
        ax.set_title('Receptor RoG', fontsize=11, fontweight='bold')
        ax.grid(True, alpha=0.3)
        
        stats_text = f'Mean: {np.mean(rog):.3f} nm\nStd: {np.std(rog):.3f} nm'
        ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, fontsize=9,
                verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        plot_idx += 1
    
    # Radius of Gyration - Peptide
    if 'rog_peptide' in data_dict and plot_idx < len(positions):
        time, rog = data_dict['rog_peptide']
        ax = fig.add_subplot(gs[positions[plot_idx][0], positions[plot_idx][1]])
        ax.plot(time, rog, color='#2ECC71', linewidth=0.8, alpha=0.8)
        ax.axhline(y=np.mean(rog), color='red', linestyle='--', linewidth=1.5)
        ax.set_xlabel('Time (ns)', fontsize=10)
        ax.set_ylabel('RoG (nm)', fontsize=10)
        ax.set_title('Peptide RoG', fontsize=11, fontweight='bold')
        ax.grid(True, alpha=0.3)
        
        stats_text = f'Mean: {np.mean(rog):.3f} nm\nStd: {np.std(rog):.3f} nm'
        ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, fontsize=9,
                verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        plot_idx += 1
    
    # RMSF receptor (simplified)
    if 'rmsf_receptor' in data_dict and plot_idx < len(positions):
        residues, rmsf = data_dict['rmsf_receptor']
        ax = fig.add_subplot(gs[positions[plot_idx][0], positions[plot_idx][1]])
        ax.bar(residues, rmsf, color='#3498DB', alpha=0.8, width=1.0)
        ax.axhline(y=np.mean(rmsf), color='red', linestyle='--', linewidth=1.5)
        ax.set_xlabel('Residue', fontsize=10)
        ax.set_ylabel('RMSF (nm)', fontsize=10)
        ax.set_title('Receptor RMSF', fontsize=11, fontweight='bold')
        ax.grid(True, alpha=0.3, axis='y')
        
        stats_text = f'Mean: {np.mean(rmsf):.3f} nm\nMax: {np.max(rmsf):.3f} nm'
        ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, fontsize=9,
                verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        plot_idx += 1
    
    # RMSF peptide (simplified)
    if 'rmsf_peptide' in data_dict and plot_idx < len(positions):
        residues, rmsf = data_dict['rmsf_peptide']
        ax = fig.add_subplot(gs[positions[plot_idx][0], positions[plot_idx][1]])
        ax.bar(residues, rmsf, color='#2ECC71', alpha=0.8, width=1.0)
        ax.axhline(y=np.mean(rmsf), color='red', linestyle='--', linewidth=1.5)
        ax.set_xlabel('Residue', fontsize=10)
        ax.set_ylabel('RMSF (nm)', fontsize=10)
        ax.set_title('Peptide RMSF', fontsize=11, fontweight='bold')
        ax.grid(True, alpha=0.3, axis='y')
        
        stats_text = f'Mean: {np.mean(rmsf):.3f} nm\nMax: {np.max(rmsf):.3f} nm'
        ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, fontsize=9,
                verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        plot_idx += 1
    
    fig.suptitle(f'{peptide_name} - MD Simulation Analysis Summary', 
                 fontsize=14, fontweight='bold', y=0.98)
    
    output_file = output_dir / f'{peptide_name}_summary.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"✓ Saved: {output_file}")
    plt.close()


def main():
    parser = argparse.ArgumentParser(
        description='Generate MD analysis plots (RMSD, RMSF, RoG, H-bonds)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s /path/to/simulation_dir
  %(prog)s /path/to/simulation_dir --peptide "My Peptide"
  %(prog)s /path/to/simulation_dir --summary-only

For extended analyses (DSSP, MinDist), run:
  python plot_simulations_extended.py /path/to/simulation_dir
"""
    )
    parser.add_argument('input_dir', type=str, help='Directory containing simulation data')
    parser.add_argument('--peptide', type=str, default='', help='Peptide name for plot titles')
    parser.add_argument('--summary-only', action='store_true', help='Generate only summary plot')
    args = parser.parse_args()
    
    input_dir = Path(args.input_dir)
    if not input_dir.exists():
        print(f"Error: Directory '{input_dir}' does not exist!")
        sys.exit(1)
    
    peptide_name = args.peptide if args.peptide else input_dir.name
    output_dir = get_output_dir(input_dir)
    
    print("=" * 60)
    print("MD Simulation Plot Generator")
    print("(RMSD, RMSF, RoG, H-bonds)")
    print("=" * 60)
    print(f"Input directory: {input_dir}")
    print(f"Output directory: {output_dir}")
    print(f"Peptide name: {peptide_name}")
    print("=" * 60 + "\n")
    
    # Find XVG files
    print("Searching for data files...")
    
    data_dict = {}
    plots_created = []
    
    # Load RMSD receptor backbone
    rmsd_bb_path = find_xvg_file(input_dir, 'rmsd_receptor_backbone.xvg')
    time_bb, rmsd_bb = None, None
    if rmsd_bb_path:
        print(f"  ✓ Found receptor backbone RMSD: {rmsd_bb_path.relative_to(input_dir)}")
        time_bb, data_bb, _ = parse_xvg(rmsd_bb_path)
        rmsd_bb = data_bb[:, 0] if len(data_bb.shape) > 1 else data_bb
        data_dict['rmsd_backbone'] = (time_bb, rmsd_bb)
    else:
        print("  ✗ rmsd_receptor_backbone.xvg not found")
    
    # Load RMSD peptide (binding - relative to receptor)
    # Try both naming conventions: rmsd_peptide_binding.xvg (new) and rmsd_peptide_backbone.xvg (old)
    rmsd_pep_path = find_xvg_file(input_dir, 'rmsd_peptide_binding.xvg')
    if not rmsd_pep_path:
        rmsd_pep_path = find_xvg_file(input_dir, 'rmsd_peptide_backbone.xvg')
    time_pep, rmsd_pep = None, None
    if rmsd_pep_path:
        print(f"  ✓ Found peptide RMSD (binding): {rmsd_pep_path.relative_to(input_dir)}")
        time_pep, data_pep, _ = parse_xvg(rmsd_pep_path)
        rmsd_pep = data_pep[:, 0] if len(data_pep.shape) > 1 else data_pep
        data_dict['rmsd_peptide'] = (time_pep, rmsd_pep)
    else:
        print("  ✗ rmsd_peptide_binding.xvg (or rmsd_peptide_backbone.xvg) not found")
    
    # Load RMSD peptide internal (conformational flexibility)
    rmsd_pep_int_path = find_xvg_file(input_dir, 'rmsd_peptide_internal.xvg')
    time_pep_int, rmsd_pep_int = None, None
    if rmsd_pep_int_path:
        print(f"  ✓ Found peptide internal RMSD: {rmsd_pep_int_path.relative_to(input_dir)}")
        time_pep_int, data_pep_int, _ = parse_xvg(rmsd_pep_int_path)
        rmsd_pep_int = data_pep_int[:, 0] if len(data_pep_int.shape) > 1 else data_pep_int
        data_dict['rmsd_peptide_internal'] = (time_pep_int, rmsd_pep_int)
    else:
        print("  ✗ rmsd_peptide_internal.xvg not found")
    
    # Load H-bonds
    hbond_path = find_xvg_file(input_dir, 'hbonds.xvg')
    time_hb, hbonds = None, None
    if hbond_path:
        print(f"  ✓ Found hydrogen bonds: {hbond_path.relative_to(input_dir)}")
        time_hb, data_hb, _ = parse_xvg(hbond_path)
        hbonds = data_hb[:, 0] if len(data_hb.shape) > 1 else data_hb
        data_dict['hbonds'] = (time_hb, hbonds)
    else:
        print("  ✗ hbonds.xvg not found")
    
    # Load Radius of Gyration - Receptor
    rog_rec_path = find_xvg_file(input_dir, 'rog_receptor.xvg')
    time_rog_rec, rog_rec = None, None
    if rog_rec_path:
        print(f"  ✓ Found receptor radius of gyration: {rog_rec_path.relative_to(input_dir)}")
        time_rog_rec, data_rog_rec, _ = parse_xvg(rog_rec_path)
        rog_rec = data_rog_rec[:, 0] if len(data_rog_rec.shape) > 1 else data_rog_rec
        data_dict['rog_receptor'] = (time_rog_rec, rog_rec)
    else:
        print("  ✗ rog_receptor.xvg not found")
    
    # Load Radius of Gyration - Peptide
    rog_pep_path = find_xvg_file(input_dir, 'rog_peptide.xvg')
    time_rog_pep, rog_pep = None, None
    if rog_pep_path:
        print(f"  ✓ Found peptide radius of gyration: {rog_pep_path.relative_to(input_dir)}")
        time_rog_pep, data_rog_pep, _ = parse_xvg(rog_pep_path)
        rog_pep = data_rog_pep[:, 0] if len(data_rog_pep.shape) > 1 else data_rog_pep
        data_dict['rog_peptide'] = (time_rog_pep, rog_pep)
    else:
        print("  ✗ rog_peptide.xvg not found")
    
    # Load RMSF receptor
    rmsf_rec_path = find_xvg_file(input_dir, 'rmsf_receptor.xvg')
    res_rmsf_rec, rmsf_rec = None, None
    if rmsf_rec_path:
        print(f"  ✓ Found receptor RMSF: {rmsf_rec_path.relative_to(input_dir)}")
        res_rmsf_rec, data_rmsf_rec, _ = parse_xvg(rmsf_rec_path)
        rmsf_rec = data_rmsf_rec[:, 0] if len(data_rmsf_rec.shape) > 1 else data_rmsf_rec
        data_dict['rmsf_receptor'] = (res_rmsf_rec, rmsf_rec)
    else:
        print("  ✗ rmsf_receptor.xvg not found")
    
    # Load RMSF peptide
    rmsf_pep_path = find_xvg_file(input_dir, 'rmsf_peptide.xvg')
    res_rmsf_pep, rmsf_pep = None, None
    if rmsf_pep_path:
        print(f"  ✓ Found peptide RMSF: {rmsf_pep_path.relative_to(input_dir)}")
        res_rmsf_pep, data_rmsf_pep, _ = parse_xvg(rmsf_pep_path)
        rmsf_pep = data_rmsf_pep[:, 0] if len(data_rmsf_pep.shape) > 1 else data_rmsf_pep
        data_dict['rmsf_peptide'] = (res_rmsf_pep, rmsf_pep)
    else:
        print("  ✗ rmsf_peptide.xvg not found")
    
    print("\n" + "=" * 60)
    print("Generating plots...")
    print("=" * 60 + "\n")
    
    # Generate RMSD plots
    if rmsd_bb is not None and not args.summary_only:
        try:
            plot_rmsd_analysis(time_bb, rmsd_bb, time_pep, rmsd_pep, 
                             time_pep_int, rmsd_pep_int, peptide_name, output_dir)
            plots_created.append("RMSD analysis")
        except Exception as e:
            print(f"✗ Error creating RMSD plot: {e}")
    
    # Generate H-bond plots
    if hbonds is not None and not args.summary_only:
        try:
            plot_hbond_analysis(time_hb, hbonds, peptide_name, output_dir)
            plots_created.append("H-bond analysis")
        except Exception as e:
            print(f"✗ Error creating H-bond plot: {e}")
    
    # Generate RMSF plots
    if (rmsf_rec is not None or rmsf_pep is not None) and not args.summary_only:
        try:
            plot_rmsf_analysis(res_rmsf_rec, rmsf_rec, res_rmsf_pep, rmsf_pep, peptide_name, output_dir)
            plots_created.append("RMSF analysis")
        except Exception as e:
            print(f"✗ Error creating RMSF plot: {e}")
    
    # Generate RoG plots - Receptor
    if rog_rec is not None and not args.summary_only:
        try:
            plot_rog_analysis(time_rog_rec, rog_rec, f"{peptide_name} (Receptor)", output_dir)
            plots_created.append("Receptor radius of gyration")
        except Exception as e:
            print(f"✗ Error creating receptor RoG plot: {e}")
    
    # Generate RoG plots - Peptide
    if rog_pep is not None and not args.summary_only:
        try:
            plot_rog_analysis(time_rog_pep, rog_pep, f"{peptide_name} (Peptide)", output_dir)
            plots_created.append("Peptide radius of gyration")
        except Exception as e:
            print(f"✗ Error creating peptide RoG plot: {e}")
    
    # Generate combined summary
    if data_dict:
        try:
            plot_combined_summary(data_dict, peptide_name, output_dir)
            plots_created.append("Combined summary")
        except Exception as e:
            print(f"✗ Error creating summary plot: {e}")
    
    print("\n" + "=" * 60)
    print("Summary")
    print("=" * 60)
    if plots_created:
        print(f"✓ Successfully generated {len(plots_created)} plot(s):")
        for plot in plots_created:
            print(f"  - {plot}")
        print(f"\nPlots saved to: {output_dir}")
    else:
        print("✗ No plots were generated.")
        print("\nExpected files from analyze_md.sh:")
        print("  - rmsd_receptor_backbone.xvg (receptor stability)")
        print("  - rmsd_peptide_backbone.xvg (binding stability)")
        print("  - rmsd_peptide_internal.xvg (conformational flexibility)")
        print("  - rmsf_receptor.xvg (receptor flexibility)")
        print("  - rmsf_peptide.xvg (peptide flexibility)")
        print("  - rog_receptor.xvg (receptor radius of gyration)")
        print("  - rog_peptide.xvg (peptide radius of gyration)")
        print("  - hbonds.xvg (hydrogen bonds)")
    
    print("\nFor extended analyses (DSSP, MinDist), run:")
    print(f"  python plot_simulations_extended.py {input_dir}")
    print("=" * 60)


if __name__ == "__main__":
    main()
