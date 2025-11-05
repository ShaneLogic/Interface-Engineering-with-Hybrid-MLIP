#!/usr/bin/env python3
"""
fix_interface_layers.py

This script adds Selective Dynamics constraints to a heterojunction POSCAR file,
relaxing the two closest layers at the interface (one from each material) and fixing all other layers.

The script:
1. Reads a combined heterojunction POSCAR file
2. Identifies the interface position (boundary between two materials)
3. Finds the two closest atomic layers to the interface (one from each side)
4. Adds Selective Dynamics flags: T T T (free/relax) for interface layers, F F F (fixed) for others
5. Writes a new POSCAR file with Selective Dynamics

Author: Generated for heterojunction interface relaxation
Requirements: pymatgen, numpy
"""

import argparse
import numpy as np
from pymatgen.core import Structure
from pymatgen.io.vasp import Poscar


def identify_interface_z(structure, method='density_gap'):
    """
    Identify the z-coordinate of the interface between two materials.
    
    Args:
        structure: pymatgen Structure object
        method: Method to identify interface ('density_gap' or 'median')
    
    Returns:
        float: z-coordinate of the interface
    """
    z_coords = structure.cart_coords[:, 2]
    
    if method == 'density_gap':
        # Find z-coordinate with minimum atomic density (gap between materials)
        z_min, z_max = z_coords.min(), z_coords.max()
        z_range = z_max - z_min
        n_bins = max(50, len(structure) // 10)  # Adaptive binning
        bins = np.linspace(z_min, z_max, n_bins)
        hist, bin_edges = np.histogram(z_coords, bins=bins)
        
        # Find bin with minimum density (gap)
        min_density_idx = np.argmin(hist)
        interface_z = (bin_edges[min_density_idx] + bin_edges[min_density_idx + 1]) / 2.0
        
    elif method == 'median':
        # Use median as interface (simple but may not be accurate)
        interface_z = np.median(z_coords)
    
    else:
        # Default: find the largest gap in z-coordinates
        sorted_z = np.sort(z_coords)
        gaps = np.diff(sorted_z)
        max_gap_idx = np.argmax(gaps)
        interface_z = (sorted_z[max_gap_idx] + sorted_z[max_gap_idx + 1]) / 2.0
    
    return interface_z


def find_interface_layers(structure, interface_z, n_layers_per_side=1, layer_thickness=2.0):
    """
    Find atomic layers closest to the interface on both sides.
    Selects n_layers_per_side distinct atomic layers from each side (TiO2 and perovskite).
    
    Args:
        structure: pymatgen Structure object
        interface_z: z-coordinate of the interface
        n_layers_per_side: Number of layers to relax on each side (default: 1)
        layer_thickness: Thickness threshold for grouping atoms into layers (Angstroms)
    
    Returns:
        list: List of site indices to relax (at interface layers)
    """
    z_coords = structure.cart_coords[:, 2]
    
    # Separate atoms into bottom (substrate, e.g., TiO2) and top (film, e.g., perovskite)
    bottom_indices = [i for i, z in enumerate(z_coords) if z < interface_z]
    top_indices = [i for i, z in enumerate(z_coords) if z >= interface_z]
    
    if len(bottom_indices) == 0 or len(top_indices) == 0:
        raise ValueError("Cannot identify interface: all atoms on one side")
    
    # Find z-coordinates for bottom and top materials
    bottom_z = [z_coords[i] for i in bottom_indices]
    top_z = [z_coords[i] for i in top_indices]
    
    # Round z-coordinates to identify distinct atomic layers
    # Use reasonable precision (0.01 Å) to group atoms that are essentially at the same z-level
    z_precision = 0.01  # 0.01 Å precision for layer identification
    
    bottom_z_rounded = [round(z, 2) for z in bottom_z]  # Round to 0.01 Å
    top_z_rounded = [round(z, 2) for z in top_z]
    
    # Get unique z-levels for each side (these represent distinct atomic layers)
    bottom_unique_z = sorted(set(bottom_z_rounded), reverse=True)  # Highest z first (closest to interface)
    top_unique_z = sorted(set(top_z_rounded))  # Lowest z first (closest to interface)
    
    # Select n_layers_per_side unique z-levels closest to interface
    # Each z-level represents a distinct atomic layer
    selected_bottom_z_levels = bottom_unique_z[:n_layers_per_side] if len(bottom_unique_z) >= n_layers_per_side else bottom_unique_z
    selected_top_z_levels = top_unique_z[:n_layers_per_side] if len(top_unique_z) >= n_layers_per_side else top_unique_z
    
    # Collect atoms in selected layers (these will be relaxed)
    # Use z_precision (0.01 Å) to match atoms to exact z-levels, not layer_thickness
    # layer_thickness is only used for documentation/clustering purposes, not for selection
    relaxed_indices = []
    
    # Bottom side (TiO2): select atoms in the closest n_layers_per_side layers
    for selected_z in selected_bottom_z_levels:
        for i in bottom_indices:
            z = z_coords[i]
            # Check if atom's z (rounded) matches the selected z-level exactly
            if abs(round(z, 2) - selected_z) < z_precision:
                relaxed_indices.append(i)
    
    # Top side (perovskite): select atoms in the closest n_layers_per_side layers
    for selected_z in selected_top_z_levels:
        for i in top_indices:
            z = z_coords[i]
            # Check if atom's z (rounded) matches the selected z-level exactly
            if abs(round(z, 2) - selected_z) < z_precision:
                relaxed_indices.append(i)
    
    # Remove duplicates
    relaxed_indices = list(set(relaxed_indices))
    
    return relaxed_indices


def add_selective_dynamics(poscar_file, output_file, n_layers_per_side=1, 
                          layer_thickness=2.0, interface_method='density_gap'):
    """
    Add Selective Dynamics to a POSCAR file, relaxing interface layers and fixing others.
    
    Args:
        poscar_file: Input POSCAR file path
        output_file: Output POSCAR file path
        n_layers_per_side: Number of layers to relax on each side of interface
        layer_thickness: Thickness threshold for defining a layer (Angstroms)
        interface_method: Method to identify interface ('density_gap', 'median', or 'max_gap')
    """
    print(f"Reading POSCAR file: {poscar_file}")
    structure = Structure.from_file(poscar_file)
    
    print(f"Structure contains {len(structure)} atoms")
    print(f"Lattice parameters: {structure.lattice.abc}")
    
    # Identify interface position
    interface_z = identify_interface_z(structure, method=interface_method)
    print(f"Identified interface at z = {interface_z:.4f} Å")
    
    # Find z-coordinate range
    z_coords = structure.cart_coords[:, 2]
    z_min, z_max = z_coords.min(), z_coords.max()
    print(f"Z-coordinate range: {z_min:.4f} - {z_max:.4f} Å")
    
    # Find interface layers to relax
    relaxed_indices = find_interface_layers(structure, interface_z, 
                                           n_layers_per_side=n_layers_per_side,
                                           layer_thickness=layer_thickness)
    
    print(f"Found {len(relaxed_indices)} atoms in interface layers to relax")
    
    # Create Selective Dynamics flags
    # T T T = free (relax), F F F = fixed
    selective_dynamics = []
    for i in range(len(structure)):
        if i in relaxed_indices:
            selective_dynamics.append([True, True, True])  # Free (relax interface layers)
        else:
            selective_dynamics.append([False, False, False])  # Fixed (all other layers)
    
    # Add selective dynamics to structure
    structure.add_site_property("selective_dynamics", selective_dynamics)
    
    # Write POSCAR with Selective Dynamics
    poscar = Poscar(structure)
    poscar.write_file(output_file)
    
    print(f"✓ Written POSCAR with Selective Dynamics to: {output_file}")
    print(f"  Relaxed atoms: {len(relaxed_indices)} (interface layers)")
    print(f"  Fixed atoms: {len(structure) - len(relaxed_indices)} (all other layers)")
    
    # Print some statistics
    relaxed_z_coords = [z_coords[i] for i in relaxed_indices]
    if len(relaxed_z_coords) > 0:
        print(f"  Relaxed layer z-range: {min(relaxed_z_coords):.4f} - {max(relaxed_z_coords):.4f} Å")


def main():
    parser = argparse.ArgumentParser(
        description="Add Selective Dynamics to heterojunction POSCAR, relaxing interface layers and fixing others"
    )
    parser.add_argument(
        "input_file",
        help="Input POSCAR file (combined heterojunction structure)"
    )
    parser.add_argument(
        "-o", "--output",
        dest="output_file",
        default=None,
        help="Output POSCAR file (default: input_file with '_relaxed' suffix)"
    )
    parser.add_argument(
        "-n", "--n_layers",
        dest="n_layers",
        type=int,
        default=1,
        help="Number of layers to relax on each side of interface (default: 1)"
    )
    parser.add_argument(
        "-t", "--thickness",
        dest="layer_thickness",
        type=float,
        default=2.0,
        help="Layer thickness threshold in Angstroms (default: 2.0)"
    )
    parser.add_argument(
        "-m", "--method",
        dest="interface_method",
        choices=['density_gap', 'median', 'max_gap'],
        default='density_gap',
        help="Method to identify interface position (default: density_gap)"
    )
    
    args = parser.parse_args()
    
    # Generate output filename if not provided
    if args.output_file is None:
        if args.input_file.endswith('.vasp') or args.input_file.endswith('.POSCAR'):
            args.output_file = args.input_file.replace('.vasp', '_relaxed.vasp').replace('.POSCAR', '_relaxed.POSCAR')
        else:
            args.output_file = args.input_file + '_relaxed'
    
    try:
        add_selective_dynamics(
            args.input_file,
            args.output_file,
            n_layers_per_side=args.n_layers,
            layer_thickness=args.layer_thickness,
            interface_method=args.interface_method
        )
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
        return 1
    
    return 0


if __name__ == "__main__":
    exit(main())

