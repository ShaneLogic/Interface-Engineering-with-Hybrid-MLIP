#!/usr/bin/env python3
"""
Stack two VASP-format slab files with a specified gap between them
"""

import numpy as np
import sys
import os


def read_vasp(filename):
    """Read VASP format file"""
    try:
        with open(filename, 'r') as f:
            lines = f.readlines()
        
        if len(lines) < 8:
            raise ValueError(f"File {filename} appears to be incomplete (less than 8 lines)")
        
        # Read title line
        title = lines[0].strip()
        
        # Read scale factor
        scale = float(lines[1].strip())
        
        # Read lattice vectors
        lattice = np.array([[float(x) for x in lines[i].split()] for i in range(2, 5)])
        
        # Read element types
        elements = lines[5].split()
        
        # Read number of atoms for each element
        n_atoms = [int(x) for x in lines[6].split()]
        
        # Check coordinate format
        coord_type = lines[7].strip()
        
        # Read atomic coordinates
        coords = []
        idx = 8
        for n in n_atoms:
            for _ in range(n):
                if idx >= len(lines):
                    raise ValueError(f"File {filename} has insufficient coordinate lines")
                coords.append([float(x) for x in lines[idx].split()[:3]])
                idx += 1
        
        coords = np.array(coords)
        
        return {
            'title': title,
            'scale': scale,
            'lattice': lattice,
            'elements': elements,
            'n_atoms': n_atoms,
            'coord_type': coord_type,
            'coords': coords
        }
    except FileNotFoundError:
        raise FileNotFoundError(f"File {filename} not found")
    except (ValueError, IndexError) as e:
        raise ValueError(f"Error reading VASP file {filename}: {e}")


def write_vasp(filename, data):
    """Write VASP format file"""
    with open(filename, 'w') as f:
        # Write title
        f.write(f"{data['title']}\n")
        
        # Write scale factor
        f.write(f"   {data['scale']:.15f}\n")
        
        # Write lattice vectors
        for vec in data['lattice']:
            f.write(f"  {vec[0]:18.15f} {vec[1]:18.15f} {vec[2]:18.15f}\n")
        
        # Write element types
        f.write("  " + "   ".join(data['elements']) + "\n")
        
        # Write number of atoms for each element
        f.write("  " + "  ".join([str(n) for n in data['n_atoms']]) + "\n")
        
        # Write coordinate type
        f.write(f"{data['coord_type']}\n")
        
        # Write atomic coordinates
        for coord in data['coords']:
            f.write(f"  {coord[0]:18.15f} {coord[1]:18.15f} {coord[2]:18.15f}\n")
        
        # Write empty line and velocities
        f.write("\n")
        total_atoms = sum(data['n_atoms'])
        if 'velocities' in data:
            for vel in data['velocities']:
                f.write(f"  {vel[0]:18.8E} {vel[1]:18.8E} {vel[2]:18.8E}\n")
        else:
            # Write zero velocities
            for _ in range(total_atoms):
                f.write("  0.00000000E+00  0.00000000E+00  0.00000000E+00\n")


def get_z_cartesian(coords, coord_type, c_vector):
    """
    Convert z-coordinates to Cartesian
    
    Parameters:
        coords: Array of coordinates (N, 3)
        coord_type: 'Direct' or 'Cartesian'
        c_vector: z-direction lattice vector length
    
    Returns:
        z_cartesian: Array of z-coordinates in Cartesian (N,)
    """
    if coord_type == 'Direct':
        return coords[:, 2] * c_vector
    else:
        return coords[:, 2]


def transform_coordinates(coords, coord_type, z_min_old, z_max_old, 
                          z_min_new_frac, z_max_new_frac, new_c):
    """
    Transform z-coordinates from old range to new fractional range
    
    Parameters:
        coords: Array of coordinates (N, 3)
        coord_type: 'Direct' or 'Cartesian'
        z_min_old: Minimum z in old coordinate system
        z_max_old: Maximum z in old coordinate system
        z_min_new_frac: New minimum z in fractional coordinates
        z_max_new_frac: New maximum z in fractional coordinates
        new_c: New lattice vector z-length
    
    Returns:
        new_coords: Transformed coordinates (N, 3)
    """
    new_coords = coords.copy()
    
    if coord_type == 'Direct':
        # Fractional coordinate mode
        if z_max_old != z_min_old:
            new_coords[:, 2] = z_min_new_frac + (coords[:, 2] - z_min_old) * \
                              (z_max_new_frac - z_min_new_frac) / (z_max_old - z_min_old)
        else:
            new_coords[:, 2] = z_min_new_frac
    else:
        # Cartesian coordinate mode - convert to fractional
        if z_max_old != z_min_old:
            # Map from [z_min_old, z_max_old] to [z_min_new_frac*new_c, z_max_new_frac*new_c] in Cartesian
            # Then convert to fractional
            z_new_cart_min = z_min_new_frac * new_c
            z_new_cart_max = z_max_new_frac * new_c
            z_cart = coords[:, 2]
            z_new_cart = z_new_cart_min + (z_cart - z_min_old) * \
                        (z_new_cart_max - z_new_cart_min) / (z_max_old - z_min_old)
            new_coords[:, 2] = z_new_cart / new_c
        else:
            new_coords[:, 2] = z_min_new_frac
    
    return new_coords


def merge_elements(elements1, n_atoms1, elements2, n_atoms2):
    """
    Merge element lists from two slabs, combining counts for same elements
    
    Parameters:
        elements1, n_atoms1: Elements and counts from first slab
        elements2, n_atoms2: Elements and counts from second slab
    
    Returns:
        merged_elements: List of unique elements
        merged_n_atoms: List of atom counts for each element
    """
    # Start with slab1's elements
    merged_elements = list(elements1)
    merged_n_atoms = list(n_atoms1)
    
    # Add slab2's elements, combining counts if element already exists
    for elem, n in zip(elements2, n_atoms2):
        if elem in merged_elements:
            idx = merged_elements.index(elem)
            merged_n_atoms[idx] += n
        else:
            merged_elements.append(elem)
            merged_n_atoms.append(n)
    
    return merged_elements, merged_n_atoms


def generate_output_filename(file1, file2, output_dir='relax_slab'):
    """
    Generate output filename with duplicate handling
    
    Parameters:
        file1, file2: Input file paths
        output_dir: Output directory
    
    Returns:
        output_file: Full path to output file
    """
    base1 = os.path.splitext(os.path.basename(file1))[0]
    base2 = os.path.splitext(os.path.basename(file2))[0]
    
    # Remove common suffixes to create cleaner names
    base1_clean = base1.replace('_slab_relax', '').replace('_relax', '')
    base2_clean = base2.replace('_slab_relax', '').replace('_relax', '')
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Generate base filename
    base_filename = f"{base1_clean}@{base2_clean}_stacked.vasp"
    output_file = os.path.join(output_dir, base_filename)
    
    # Check if file exists and add number suffix if needed
    if os.path.exists(output_file):
        counter = 1
        while True:
            base_name = f"{base1_clean}@{base2_clean}_stacked_{counter}.vasp"
            output_file = os.path.join(output_dir, base_name)
            if not os.path.exists(output_file):
                break
            counter += 1
        print(f"Warning: File {base_filename} already exists, using {os.path.basename(output_file)} instead")
    
    return output_file


def stack_slabs(file1, file2, gap=3.0, vacuum=20.0, output_file=None):
    """
    Stack two slab files with vacuum layers
    
    Parameters:
        file1: Path to first slab file (bottom)
        file2: Path to second slab file (top)
        gap: Gap between two slabs in Angstroms (default: 3.0)
        vacuum: Total vacuum layer thickness in Angstroms (default: 20.0)
                This will be split equally between bottom and top
        output_file: Output filename, if None will be auto-generated
    
    Returns:
        output_file: Path to the generated output file
    """
    # Read both files
    print(f"Reading {file1}...")
    slab1 = read_vasp(file1)
    
    print(f"Reading {file2}...")
    slab2 = read_vasp(file2)
    
    # Check if lattice vectors are compatible (xy plane should be the same)
    lattice1_xy = slab1['lattice'][:2, :2]
    lattice2_xy = slab2['lattice'][:2, :2]
    
    if not np.allclose(lattice1_xy, lattice2_xy, atol=1e-6):
        print("Warning: xy-plane lattice vectors of two slabs are not identical, will use first slab's lattice vectors")
    
    # Get z-direction lattice vector lengths
    c1 = np.linalg.norm(slab1['lattice'][2])
    c2 = np.linalg.norm(slab2['lattice'][2])
    
    # Convert to Cartesian coordinates (z-direction only)
    z1_cart = get_z_cartesian(slab1['coords'], slab1['coord_type'], c1)
    z2_cart = get_z_cartesian(slab2['coords'], slab2['coord_type'], c2)
    
    # Find z-coordinate ranges and thicknesses
    z1_min, z1_max = np.min(z1_cart), np.max(z1_cart)
    z2_min, z2_max = np.min(z2_cart), np.max(z2_cart)
    
    slab1_thickness = z1_max - z1_min
    slab2_thickness = z2_max - z2_min
    
    print(f"\nSlab 1 z-range: {z1_min:.3f} - {z1_max:.3f} Å (thickness: {slab1_thickness:.3f} Å)")
    print(f"Slab 2 z-range: {z2_min:.3f} - {z2_max:.3f} Å (thickness: {slab2_thickness:.3f} Å)")
    
    # Split total vacuum equally between bottom and top
    bottom_vacuum = vacuum / 2.0
    top_vacuum = vacuum / 2.0
    
    # Calculate new z-direction lattice vector length
    new_c = bottom_vacuum + slab1_thickness + gap + slab2_thickness + top_vacuum
    
    print(f"Total vacuum thickness: {vacuum:.3f} Å (split equally: {bottom_vacuum:.3f} Å bottom + {top_vacuum:.3f} Å top)")
    print(f"Total z-direction length: {new_c:.3f} Å")
    
    # Transform coordinates
    # Slab1: from [z1_min, z1_max] to [bottom_vacuum, bottom_vacuum + slab1_thickness]
    z1_min_frac = bottom_vacuum / new_c
    z1_max_frac = (bottom_vacuum + slab1_thickness) / new_c
    
    if slab1['coord_type'] == 'Direct':
        z1_frac_min = np.min(slab1['coords'][:, 2])
        z1_frac_max = np.max(slab1['coords'][:, 2])
    else:
        z1_frac_min = z1_min
        z1_frac_max = z1_max
    
    slab1_new_coords = transform_coordinates(
        slab1['coords'], slab1['coord_type'],
        z1_frac_min, z1_frac_max,
        z1_min_frac, z1_max_frac, new_c
    )
    
    # Slab2: from [z2_min, z2_max] to [bottom_vacuum + slab1_thickness + gap, ...]
    z2_min_frac = (bottom_vacuum + slab1_thickness + gap) / new_c
    z2_max_frac = (bottom_vacuum + slab1_thickness + gap + slab2_thickness) / new_c
    
    if slab2['coord_type'] == 'Direct':
        z2_frac_min = np.min(slab2['coords'][:, 2])
        z2_frac_max = np.max(slab2['coords'][:, 2])
    else:
        z2_frac_min = z2_min
        z2_frac_max = z2_max
    
    slab2_new_coords = transform_coordinates(
        slab2['coords'], slab2['coord_type'],
        z2_frac_min, z2_frac_max,
        z2_min_frac, z2_max_frac, new_c
    )
    
    # Merge elements and coordinates
    merged_elements, merged_n_atoms = merge_elements(
        slab1['elements'], slab1['n_atoms'],
        slab2['elements'], slab2['n_atoms']
    )
    merged_coords = np.vstack([slab1_new_coords, slab2_new_coords])
    
    # Create new lattice vectors
    new_lattice = slab1['lattice'].copy()
    new_lattice[2] = new_lattice[2] / np.linalg.norm(new_lattice[2]) * new_c
    
    # Generate output filename
    if output_file is None:
        output_file = generate_output_filename(file1, file2)
    else:
        # Handle output file path
        if not os.path.dirname(output_file):
            output_dir = 'relax_slab'
            os.makedirs(output_dir, exist_ok=True)
            output_file = os.path.join(output_dir, output_file)
        else:
            output_dir = os.path.dirname(output_file)
            os.makedirs(output_dir, exist_ok=True)
    
    # Create output data structure
    output_data = {
        'title': f"{slab1['title'].split()[0]}@{slab2['title'].split()[0]}",
        'scale': slab1['scale'],
        'lattice': new_lattice,
        'elements': merged_elements,
        'n_atoms': merged_n_atoms,
        'coord_type': slab1['coord_type'],
        'coords': merged_coords
    }
    
    # Write output file
    print(f"\nWriting {output_file}...")
    write_vasp(output_file, output_data)
    
    print(f"\nDone!")
    print(f"Total atoms: {sum(merged_n_atoms)}")
    print(f"New lattice z-direction length: {new_c:.3f} Å")
    print(f"  - Bottom vacuum: {bottom_vacuum:.3f} Å")
    print(f"  - Slab 1: {slab1_thickness:.3f} Å")
    print(f"  - Gap: {gap:.3f} Å")
    print(f"  - Slab 2: {slab2_thickness:.3f} Å")
    print(f"  - Top vacuum: {top_vacuum:.3f} Å")
    print(f"  - Total vacuum: {vacuum:.3f} Å (bottom + top)")
    
    return output_file


def main():
    """Main entry point for command-line usage"""
    if len(sys.argv) < 3:
        print("Usage: python stack_slabs.py <slab1_file> <slab2_file> [gap_in_angstrom] [vacuum_in_angstrom] [output_file]")
        print("Example: python stack_slabs.py relax_slab/fapbi3_slab_relax.vasp relax_slab/tio2_fa_slab_relax.vasp 3.0 20.0")
        sys.exit(1)
    
    file1 = sys.argv[1]
    file2 = sys.argv[2]
    gap = float(sys.argv[3]) if len(sys.argv) > 3 else 3.0
    vacuum = float(sys.argv[4]) if len(sys.argv) > 4 else 20.0
    output_file = sys.argv[5] if len(sys.argv) > 5 else None
    
    # Validate input files
    if not os.path.exists(file1):
        print(f"Error: File {file1} does not exist")
        sys.exit(1)
    
    if not os.path.exists(file2):
        print(f"Error: File {file2} does not exist")
        sys.exit(1)
    
    # Validate parameters
    if gap < 0:
        print("Error: Gap must be non-negative")
        sys.exit(1)
    
    if vacuum < 0:
        print("Error: Vacuum thickness must be non-negative")
        sys.exit(1)
    
    try:
        stack_slabs(file1, file2, gap, vacuum, output_file)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
