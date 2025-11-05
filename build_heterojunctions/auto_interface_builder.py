#!/usr/bin/env python3
"""
auto_interface_builder.py

Given two bulk structures (CIF/POSCAR), this script:
 - finds/optimizes primitive cells,
 - searches for low-strain commensurate interface supercells (in-plane),
 - builds slabs for given Miller indices,
 - applies strain to one or both materials per user choice,
 - stacks them with an initial separation and vacuum,
 - writes POSCARs for each slab (strained) and the combined interface.

Author: Xuan-Yan (worked example tuned for MAPbI3 / TiO2 workflows)
Requirements: pymatgen, numpy
"""

import argparse
import math
import numpy as np
from pymatgen.core import Structure, Lattice
from pymatgen.core.surface import SlabGenerator
from pymatgen.transformations.standard_transformations import SupercellTransformation
#from pymatgen.analysis.interfaces import InterfaceMatcher
from pymatgen.analysis.interfaces.coherent_interfaces import CoherentInterfaceBuilder
from pymatgen.io.vasp import Poscar

def load_structure(path):
    s = Structure.from_file(path)
    return s

def get_primitive(struct):
    """Return primitive cell if found, else the original structure."""
    try:
        prim = struct.get_primitive_structure()
        # sometimes get_primitive returns same; ensure lattice volume not pathologically changed
        return prim
    except Exception as e:
        print("Warning: get_primitive_structure failed, returning original. Err:", e)
        return struct

def build_slab(struct, miller, min_slab_size, vacuum, center_slab=True):
    """Construct a slab using SlabGenerator. Returns the first termination slab."""
    sg = SlabGenerator(struct, miller_index=tuple(map(int, miller)), 
                       min_slab_size=min_slab_size, min_vacuum_size=vacuum, center_slab=center_slab)
    slabs = sg.get_slabs()
    if len(slabs) == 0:
        raise RuntimeError(f"No slabs generated for miller {miller}")
    return slabs[0]

def search_matches(structA, structB, miller_a, miller_b, tol=0.03, max_area=500, max_match=20):
    """
    Use CoherentInterfaceBuilder to find candidate in-plane supercells that match.
    Requires bulk structures and Miller indices, not slabs.
    Returns a list of match dicts with strains and transforms.
    """
    # CoherentInterfaceBuilder needs bulk structures and Miller indices
    # structA is film, structB is substrate (can be swapped if needed)
    from pymatgen.analysis.interfaces.zsl import ZSLGenerator
    
    # Create ZSLGenerator - try with max_area parameter, fallback to default
    try:
        zslgen = ZSLGenerator(max_area=max_area)
    except TypeError:
        # If max_area is not a valid parameter, use default ZSLGenerator
        zslgen = ZSLGenerator()
    
    # CoherentInterfaceBuilder signature: (substrate_structure, film_structure, 
    #                                       film_miller, substrate_miller, zslgen)
    matcher = CoherentInterfaceBuilder(
        substrate_structure=structB,  # B is substrate
        film_structure=structA,       # A is film
        film_miller=miller_a,
        substrate_miller=miller_b,
        zslgen=zslgen
    )
    
    # Get matches from zsl_matches attribute (not get_matches method)
    zsl_matches = matcher.zsl_matches
    if not zsl_matches or len(zsl_matches) == 0:
        return []
    
    # Convert ZSLMatch objects to summary dicts
    summary = []
    for m in zsl_matches[:max_match]:
        # Filter by area only (ZSLMatch already represents valid matches)
        if m.match_area > max_area:
            continue
        
        # Get transformation matrices (2x2 for in-plane)
        film_trans = np.array(m.film_transformation)
        substrate_trans = np.array(m.substrate_transformation) if hasattr(m, 'substrate_transformation') else np.eye(2)
        
        # Calculate approximate strain vectors from transformation matrices
        # These represent the relative deformation needed for matching
        strain_A = np.array([film_trans[0,0] - 1.0, film_trans[1,1] - 1.0, film_trans[0,1]])
        strain_B = np.array([substrate_trans[0,0] - 1.0, substrate_trans[1,1] - 1.0, substrate_trans[0,1]])
        
        # Convert 2x2 transformation matrices to 3x3 supercell matrices
        # Extend 2x2 to 3x3 by adding identity for z-direction
        def extend_matrix_2d_to_3d(mat_2d):
            mat_3d = np.eye(3)
            mat_3d[0:2, 0:2] = mat_2d
            # Round to nearest integer for supercell matrix
            return np.round(mat_3d).astype(int)
        
        supercell_a = extend_matrix_2d_to_3d(film_trans)
        supercell_b = extend_matrix_2d_to_3d(substrate_trans)
        
        summary.append({
            "strain_A": strain_A,
            "strain_B": strain_B,
            "area": m.match_area,
            "supercell_a": supercell_a,
            "supercell_b": supercell_b,
            "tilt": 0.0,  # ZSLMatch doesn't have tilt info directly
            "match_obj": m
        })
    
    return summary

def compute_inplane_norms(latt):
    """Return norms of the first two lattice vectors (in-plane)."""
    mat = np.array(latt.matrix)
    a = np.linalg.norm(mat[0])
    b = np.linalg.norm(mat[1])
    return a, b

def apply_inplane_strain(slab, scale_x, scale_y):
    """
    Return a copy of slab with in-plane lattice vectors scaled by scale_x, scale_y.
    Preserves the c vector magnitude and its relative orientation to maintain slab thickness.
    """
    new_lat = slab.lattice.matrix.copy()
    # Scale in-plane vectors
    new_lat[0] = new_lat[0] * scale_x
    new_lat[1] = new_lat[1] * scale_y
    # Keep c vector unchanged to preserve slab thickness
    # The c vector will be adjusted later in align_and_stack_ordered if needed
    new_lattice = Lattice(new_lat)
    s2 = Structure(new_lattice, slab.species, slab.frac_coords)
    return s2

def align_and_stack_ordered(slab_bottom, slab_top, separation=3.2, vacuum=20.0):
    """
    Align and stack two slabs with proper lattice matching.
    Uses the original slab lattices directly without modification to preserve structure integrity.
    """
    # Copy slabs to avoid modifying originals
    bottom = slab_bottom.copy()
    top = slab_top.copy()
    
    # Align bottom slab: ensure min z = 0
    bottom_cart = bottom.cart_coords
    minz_b = bottom_cart[:,2].min()
    bottom.translate_sites(list(range(len(bottom))), [0,0, -minz_b])
    
    # Compute bottom thickness (max z)
    bottom_cart = bottom.cart_coords
    maxz_b = bottom_cart[:,2].max()
    
    # Align top slab: translate so minz_t sits at maxz_b + separation
    top_cart = top.cart_coords
    minz_t = top_cart[:,2].min()
    dz = maxz_b + separation - minz_t
    top.translate_sites(list(range(len(top))), [0,0, dz])
    
    # Determine combined cell z dimension
    maxz_top = top.cart_coords[:,2].max()
    total_z = maxz_top + vacuum
    
    # Use bottom's lattice as the reference (it should match top after straining)
    # Get in-plane vectors from bottom
    bottom_lat = bottom.lattice.matrix
    a_vec = bottom_lat[0]
    b_vec = bottom_lat[1]
    
    # Verify top's in-plane vectors match (should be true after straining)
    top_lat = top.lattice.matrix
    if np.linalg.norm(a_vec - top_lat[0]) > 0.01 or np.linalg.norm(b_vec - top_lat[1]) > 0.01:
        print("Warning: In-plane lattice vectors don't match exactly. Using bottom slab vectors.")
        # Create new lattice for top using bottom's in-plane vectors but preserving top's c
        new_top_lat = top_lat.copy()
        new_top_lat[0] = a_vec
        new_top_lat[1] = b_vec
        top_lattice_new = Lattice(new_top_lat)
        # Recalculate fractional coordinates with new lattice
        top_cart_coords = top.cart_coords
        top = Structure(top_lattice_new, top.species, 
                       [top_lattice_new.get_fractional_coords(coord) for coord in top_cart_coords],
                       coords_are_cartesian=False)
    
    # Calculate surface normal from in-plane vectors
    ab_normal = np.cross(a_vec, b_vec)
    ab_normal_norm = np.linalg.norm(ab_normal)
    if ab_normal_norm < 1e-10:
        raise ValueError("In-plane vectors are parallel or degenerate")
    ab_normal_unit = ab_normal / ab_normal_norm
    
    # Create c vector perpendicular to ab plane with correct magnitude
    c_vec = ab_normal_unit * total_z
    
    # Build combined lattice using bottom's in-plane vectors and perpendicular c
    combined_lat_matrix = np.vstack([a_vec, b_vec, c_vec])
    combined_lattice = Lattice(combined_lat_matrix)
    
    # Create combined structure
    combined = Structure(combined_lattice, [], [])
    
    # Add bottom slab sites - use cartesian coordinates then convert to fractional
    for site in bottom:
        frac_coords = combined_lattice.get_fractional_coords(site.coords)
        # Wrap fractional coordinates to [0, 1)
        frac_coords = frac_coords % 1.0
        combined.append(site.species_string, frac_coords, coords_are_cartesian=False)
    
    # Add top slab sites - use cartesian coordinates then convert to fractional
    for site in top:
        frac_coords = combined_lattice.get_fractional_coords(site.coords)
        # Wrap fractional coordinates to [0, 1)
        frac_coords = frac_coords % 1.0
        combined.append(site.species_string, frac_coords, coords_are_cartesian=False)
    
    # Sort sites by z-coordinate for better ordering
    combined.sort(key=lambda site: site.frac_coords[2])
    
    return bottom, top, combined

def write_poscars(bottom, top, combined, prefix="interface"):
    Poscar(bottom).write_file(f"{prefix}_bottom_strained.vasp")
    Poscar(top).write_file(f"{prefix}_top_strained.vasp")
    Poscar(combined).write_file(f"{prefix}_combined.vasp")
    print("Wrote POSCARs:", f"{prefix}_bottom_strained.vasp", f"{prefix}_top_strained.vasp", f"{prefix}_combined.vasp")

def estimate_atoms_for_match(structA, structB, miller_a, miller_b, supercell_a, supercell_b, 
                             slab_thickness_a, slab_thickness_b):
    """Estimate total number of atoms for a given match."""
    # Build temporary slabs to estimate atom count
    try:
        slabA_temp = build_slab(structA, miller_a, slab_thickness_a, 0, center_slab=False)
        slabB_temp = build_slab(structB, miller_b, slab_thickness_b, 0, center_slab=False)
        
        # Apply supercell transformations
        transA = [[int(supercell_a[0,0]), int(supercell_a[0,1]), int(supercell_a[0,2])],
                  [int(supercell_a[1,0]), int(supercell_a[1,1]), int(supercell_a[1,2])],
                  [int(supercell_a[2,0]), int(supercell_a[2,1]), int(supercell_a[2,2])]]
        transB = [[int(supercell_b[0,0]), int(supercell_b[0,1]), int(supercell_b[0,2])],
                  [int(supercell_b[1,0]), int(supercell_b[1,1]), int(supercell_b[1,2])],
                  [int(supercell_b[2,0]), int(supercell_b[2,1]), int(supercell_b[2,2])]]
        
        slabA_super = SupercellTransformation(transA).apply_transformation(slabA_temp)
        slabB_super = SupercellTransformation(transB).apply_transformation(slabB_temp)
        
        total_atoms = len(slabA_super) + len(slabB_super)
        return total_atoms
    except:
        # If estimation fails, return a large number
        return 10000

def build_interface_from_builder(structA, structB, miller_a, miller_b, slab_thickness_a, slab_thickness_b, 
                                  vacuum, sep, max_area, zslgen, max_atoms=400):
    """
    Build ordered interface using CoherentInterfaceBuilder's get_interfaces method.
    This ensures proper alignment and high symmetry.
    Automatically adjusts thickness if atom count exceeds max_atoms.
    """
    from pymatgen.analysis.interfaces.zsl import ZSLGenerator
    
    builder = CoherentInterfaceBuilder(
        substrate_structure=structB,
        film_structure=structA,
        film_miller=miller_a,
        substrate_miller=miller_b,
        zslgen=zslgen
    )
    
    if not builder.zsl_matches or len(builder.zsl_matches) == 0:
        return None, None, None
    
    # Get terminations
    terminations = builder.terminations
    if len(terminations) == 0:
        return None, None, None
    
    # Try different terminations and thicknesses to find one within atom limit
    best_interface = None
    best_atoms = float('inf')
    best_thickness_a = slab_thickness_a
    best_thickness_b = slab_thickness_b
    best_termination = None
    
    # Try reducing thickness if needed - start with smaller thickness
    # Calculate initial atom count to determine starting factor
    initial_interface = None
    try:
        initial_interfaces = list(builder.get_interfaces(
            termination=terminations[0],
            gap=sep,
            vacuum_over_film=vacuum,
            film_thickness=slab_thickness_a,
            substrate_thickness=slab_thickness_b,
            in_layers=False
        ))
        if len(initial_interfaces) > 0:
            initial_interface = initial_interfaces[0]
            initial_atoms = len(initial_interface)
            if initial_atoms > max_atoms:
                # Calculate aggressive reduction factor
                thickness_factor = (max_atoms / initial_atoms) * 0.8  # 80% to be safe
            else:
                thickness_factor = 1.0
        else:
            thickness_factor = 0.5  # Start with 50% if we can't estimate
    except:
        thickness_factor = 0.5
    
    max_attempts = 8
    min_thickness_a = 3.0  # Minimum thickness in Angstroms
    min_thickness_b = 3.0
    
    for attempt in range(max_attempts):
        current_thickness_a = max(slab_thickness_a * thickness_factor, min_thickness_a)
        current_thickness_b = max(slab_thickness_b * thickness_factor, min_thickness_b)
        
        for termination in terminations:
            try:
                interfaces = list(builder.get_interfaces(
                    termination=termination,
                    gap=sep,
                    vacuum_over_film=vacuum,
                    film_thickness=current_thickness_a,
                    substrate_thickness=current_thickness_b,
                    in_layers=False  # Use Angstroms
                ))
                
                if len(interfaces) > 0:
                    interface = interfaces[0]
                    n_atoms = len(interface)
                    
                    if n_atoms <= max_atoms and n_atoms < best_atoms:
                        best_interface = interface
                        best_atoms = n_atoms
                        best_thickness_a = current_thickness_a
                        best_thickness_b = current_thickness_b
                        best_termination = termination
            except:
                continue
        
        if best_interface is not None:
            break
        
        # Reduce thickness more aggressively for next attempt
        thickness_factor *= 0.6
    
    if best_interface is None:
        print(f"Warning: Could not find interface with <={max_atoms} atoms. Trying smallest available...")
        # Try with very small thickness to find minimum possible
        for termination in terminations:
            for thin_factor in [0.3, 0.2, 0.15, 0.1]:
                try:
                    thin_a = max(slab_thickness_a * thin_factor, min_thickness_a)
                    thin_b = max(slab_thickness_b * thin_factor, min_thickness_b)
                    interfaces = list(builder.get_interfaces(
                        termination=termination,
                        gap=sep,
                        vacuum_over_film=vacuum,
                        film_thickness=thin_a,
                        substrate_thickness=thin_b,
                        in_layers=False
                    ))
                    if len(interfaces) > 0:
                        interface = interfaces[0]
                        n_atoms = len(interface)
                        if n_atoms < best_atoms:
                            best_interface = interface
                            best_atoms = n_atoms
                            best_thickness_a = thin_a
                            best_thickness_b = thin_b
                            best_termination = termination
                            if n_atoms <= max_atoms:
                                break
                except:
                    continue
            if best_interface is not None and best_atoms <= max_atoms:
                break
        
        if best_interface is None:
            # Last resort: use first available
            termination = terminations[0]
            interfaces = list(builder.get_interfaces(
                termination=termination,
                gap=sep,
                vacuum_over_film=vacuum,
                film_thickness=slab_thickness_a,
                substrate_thickness=slab_thickness_b,
                in_layers=False
            ))
            if len(interfaces) > 0:
                best_interface = interfaces[0]
                best_atoms = len(best_interface)
                best_termination = termination
            else:
                return None, None, None
    
    print(f"Using termination: {best_termination}")
    print(f"Interface contains {best_atoms} atoms (target: <={max_atoms})")
    if best_atoms > max_atoms:
        print(f"Warning: Atom count ({best_atoms}) exceeds limit ({max_atoms})")
    
    interface = best_interface
    
    # The interface structure from CoherentInterfaceBuilder.get_interfaces is already complete
    # It contains both slabs properly aligned
    combined = interface
    
    # Extract bottom and top slabs from interface for separate POSCAR files
    # Find the interface plane (approximate)
    z_coords = combined.cart_coords[:, 2]
    # Use median or density gap method for better interface identification
    interface_z = np.median(z_coords)
    
    bottom_sites = [i for i, site in enumerate(combined) if site.coords[2] < interface_z]
    top_sites = [i for i, site in enumerate(combined) if site.coords[2] >= interface_z]
    
    # Create separate structures using the same lattice as combined
    # This preserves the structure integrity
    bottom = Structure(combined.lattice, 
                      [combined[i].species_string for i in bottom_sites],
                      [combined[i].frac_coords for i in bottom_sites],
                      coords_are_cartesian=False)
    
    top = Structure(combined.lattice,
                   [combined[i].species_string for i in top_sites],
                   [combined[i].frac_coords for i in top_sites],
                   coords_are_cartesian=False)
    
    return bottom, top, combined

def auto_run(a_file, b_file, miller_a, miller_b, slab_thickness_a, slab_thickness_b, vacuum, sep, tol, max_area, strain_target, use_builder_interface=False, max_atoms=400):
    """Main orchestrator."""
    print("Loading structures...")
    A = load_structure(a_file)
    B = load_structure(b_file)
    print("Optimizing / getting primitive cells...")
    Aprim = get_primitive(A)
    Bprim = get_primitive(B)
    print("Primitive cell A lattice:", Aprim.lattice.abc)
    print("Primitive cell B lattice:", Bprim.lattice.abc)

    # Option to use CoherentInterfaceBuilder's get_interfaces for ordered interface
    if use_builder_interface:
        print("Building ordered interface using CoherentInterfaceBuilder.get_interfaces...")
        from pymatgen.analysis.interfaces.zsl import ZSLGenerator
        try:
            zslgen = ZSLGenerator(max_area=max_area)
        except TypeError:
            zslgen = ZSLGenerator()
        
        bottom, top, combined = build_interface_from_builder(
            Aprim, Bprim, miller_a, miller_b,
            slab_thickness_a, slab_thickness_b, vacuum, sep, max_area, zslgen, max_atoms
        )
        
        if combined is not None:
            print("Successfully built ordered interface using CoherentInterfaceBuilder.")
            write_poscars(bottom, top, combined, prefix="auto_interface")
            return
    
    # Search for matches using bulk structures and Miller indices
    print("Searching for low-strain matches (this may take a moment)...")
    matches = search_matches(Aprim, Bprim, miller_a, miller_b, tol=tol, max_area=max_area)
    
    if len(matches) == 0:
        raise RuntimeError("No matches found within tolerance. Increase tol or max_area.")
    
    # Build slabs after finding matches (needed for final structure construction)
    print("Building slabs (unstrained) for final structure...")
    slabA = build_slab(Aprim, miller_a, slab_thickness_a, vacuum)
    slabB = build_slab(Bprim, miller_b, slab_thickness_b, vacuum)
    
    # First, estimate atom counts for all matches with reduced thickness if needed
    # Calculate optimal thickness based on best match atom density
    def estimate_optimal_thickness(match, target_atoms):
        """Estimate optimal thickness to achieve target atom count."""
        est_atoms_full = estimate_atoms_for_match(
            Aprim, Bprim, miller_a, miller_b,
            match["supercell_a"], match["supercell_b"],
            slab_thickness_a, slab_thickness_b
        )
        
        if est_atoms_full <= target_atoms:
            return slab_thickness_a, slab_thickness_b
        
        # Calculate reduction factor
        reduction = (target_atoms / est_atoms_full) * 0.75  # 75% safety margin
        optimal_a = max(slab_thickness_a * reduction, 2.5)
        optimal_b = max(slab_thickness_b * reduction, 2.5)
        
        return optimal_a, optimal_b
    
    # Score matches prioritizing atom count, then strain
    # For DFT efficiency, atom count is more important than perfect strain matching
    def match_score(m):
        sa = np.array(m["strain_A"]) if m["strain_A"] is not None else np.array([0, 0, 0])
        sb = np.array(m["strain_B"]) if m["strain_B"] is not None else np.array([0, 0, 0])
        strain_magnitude = np.linalg.norm(sa) + np.linalg.norm(sb)
        
        # Estimate atom count with optimal thickness for this match
        opt_thick_a, opt_thick_b = estimate_optimal_thickness(m, max_atoms)
        est_atoms = estimate_atoms_for_match(
            Aprim, Bprim, miller_a, miller_b,
            m["supercell_a"], m["supercell_b"],
            opt_thick_a, opt_thick_b
        )
        
        # Calculate supercell volume to estimate atom density
        det_a = abs(np.linalg.det(m["supercell_a"]))
        det_b = abs(np.linalg.det(m["supercell_b"]))
        supercell_factor = max(det_a, det_b)  # Larger supercell = more atoms
        
        # Prioritize atom count heavily - if over limit, heavily penalize
        if est_atoms > max_atoms:
            # Very heavy penalty for exceeding limit
            atom_penalty = (est_atoms - max_atoms) * 2.0 + (est_atoms / max_atoms - 1.0) * 10.0
        else:
            # Small bonus for being well under limit
            atom_penalty = -(max_atoms - est_atoms) / max_atoms * 0.1
        
        # Area penalty (smaller area generally means smaller supercell)
        area_penalty = m["area"] / 500.0
        
        # Strain penalty (less important than atom count)
        strain_penalty = strain_magnitude * 0.1
        
        # Primary score: atom count first, then strain
        primary_score = atom_penalty + strain_penalty + area_penalty
        
        return primary_score, est_atoms, strain_magnitude
    
    # Score all matches
    matches_scored = [(m, match_score(m)) for m in matches]
    
    # Sort primarily by atom count (within limit), then by total score
    matches_scored.sort(key=lambda x: (
        0 if x[1][1] <= max_atoms else 1,  # Within limit first
        x[1][1],  # Then by atom count
        x[1][0]   # Then by score
    ))
    
    # Prioritize matches within atom limit
    valid_matches = [m for m, (score, atoms, strain) in matches_scored if atoms <= max_atoms]
    
    if len(valid_matches) > 0:
        matches_sorted = valid_matches
        print(f"Found {len(valid_matches)} matches within atom limit ({max_atoms})")
        # Among valid matches, prefer smallest atom count
        matches_sorted.sort(key=lambda m: match_score(m)[1])
    else:
        # If no matches within limit, find the one that's closest to limit
        print(f"Warning: No matches within atom limit ({max_atoms}). Finding closest match...")
        matches_sorted = [m for m, _ in matches_scored]
        # Sort by how close to limit (prefer smaller)
        matches_sorted.sort(key=lambda m: match_score(m)[1])
    
    best = matches_sorted[0]
    (best_score, best_atoms, best_strain) = match_score(best)
    
    # Calculate optimal thickness for selected match
    optimal_thick_a, optimal_thick_b = estimate_optimal_thickness(best, max_atoms)
    
    print(f"Selected match with estimated {best_atoms} atoms (strain: {best_strain:.4f})")
    print(f"Optimal thickness for this match: A={optimal_thick_a:.2f} Å, B={optimal_thick_b:.2f} Å")
    
    # Use optimal thickness if significantly different
    if abs(optimal_thick_a - slab_thickness_a) > 0.5 or abs(optimal_thick_b - slab_thickness_b) > 0.5:
        print(f"Adjusting thickness to meet atom limit:")
        print(f"  A: {slab_thickness_a:.2f} -> {optimal_thick_a:.2f} Å")
        print(f"  B: {slab_thickness_b:.2f} -> {optimal_thick_b:.2f} Å")
        slab_thickness_a = optimal_thick_a
        slab_thickness_b = optimal_thick_b
    
    mobj = best["match_obj"]
    print("Best match summary:")
    print(" area:", best["area"])
    print(" supercell A matrix:\n", best["supercell_a"])
    print(" supercell B matrix:\n", best["supercell_b"])
    print(" strain on A (approx):", best["strain_A"])
    print(" strain on B (approx):", best["strain_B"])
    # Use the InterfaceMatch object to build the supercell structures (pymatgen has helper)
    # NOTE: CoherentInterfaceBuilder provides methods to get transformed slabs via match_obj.construct_interface? 
    # For portability, we'll manually apply supercell transformations using matrices.

    # Apply supercell transform matrices to slabs (they are 3x3 integer matrices)
    supA = np.array(best["supercell_a"])
    supB = np.array(best["supercell_b"])
    # Convert int matrices into 3x3 lists
    transA = [[int(supA[0,0]), int(supA[0,1]), int(supA[0,2])],
              [int(supA[1,0]), int(supA[1,1]), int(supA[1,2])],
              [int(supA[2,0]), int(supA[2,1]), int(supA[2,2])]]
    transB = [[int(supB[0,0]), int(supB[0,1]), int(supB[0,2])],
              [int(supB[1,0]), int(supB[1,1]), int(supB[1,2])],
              [int(supB[2,0]), int(supB[2,1]), int(supB[2,2])]]
    slabA_super = SupercellTransformation(transA).apply_transformation(slabA)
    slabB_super = SupercellTransformation(transB).apply_transformation(slabB)

    # Compute in-plane lattice norms to determine actual strains
    aA_pre = compute_inplane_norms(slabA_super.lattice)[0:2]
    aB_pre = compute_inplane_norms(slabB_super.lattice)[0:2]
    print("In-plane sizes after supercell (A):", aA_pre)
    print("In-plane sizes after supercell (B):", aB_pre)

    # Decide which to strain based on user input: 'A', 'B', or 'both'
    if strain_target.lower() == 'a':
        scale_x = aB_pre[0] / aA_pre[0]
        scale_y = aB_pre[1] / aA_pre[1]
        print(f"Straining A to match B with factors x={scale_x:.6f}, y={scale_y:.6f}")
        slabA_strained = apply_inplane_strain(slabA_super, scale_x, scale_y)
        slabB_strained = slabB_super.copy()
    elif strain_target.lower() == 'b':
        scale_x = aA_pre[0] / aB_pre[0]
        scale_y = aA_pre[1] / aB_pre[1]
        print(f"Straining B to match A with factors x={scale_x:.6f}, y={scale_y:.6f}")
        slabB_strained = apply_inplane_strain(slabB_super, scale_x, scale_y)
        slabA_strained = slabA_super.copy()
    elif strain_target.lower() == 'both':
        # apply sqrt of ratio to both (split strain evenly)
        scale_x = math.sqrt(aB_pre[0] / aA_pre[0])
        scale_y = math.sqrt(aB_pre[1] / aA_pre[1])
        print(f"Splitting strain: scale factors applied to A: {scale_x:.6f},{scale_y:.6f} and to B: {1.0/scale_x:.6f},{1.0/scale_y:.6f}")
        slabA_strained = apply_inplane_strain(slabA_super, scale_x, scale_y)
        slabB_strained = apply_inplane_strain(slabB_super, 1.0/scale_x, 1.0/scale_y)
    else:
        raise ValueError("strain_target must be one of 'A', 'B', 'both'")

    # Ensure in-plane lattice vectors match exactly before stacking
    # This is critical for ordered, high-symmetry interfaces
    bottom_lat = slabB_strained.lattice.matrix
    top_lat = slabA_strained.lattice.matrix
    
    # Verify vectors match (should be true after straining, but double-check)
    if np.linalg.norm(bottom_lat[0] - top_lat[0]) > 0.01 or np.linalg.norm(bottom_lat[1] - top_lat[1]) > 0.01:
        print("Warning: In-plane vectors don't match. Adjusting top slab lattice.")
        # Preserve top's c vector while matching in-plane vectors
        new_top_lat = top_lat.copy()
        new_top_lat[0] = bottom_lat[0]
        new_top_lat[1] = bottom_lat[1]
        # Keep original c vector to preserve slab structure
        slabA_strained = Structure(Lattice(new_top_lat), slabA_strained.species, slabA_strained.frac_coords)
    
    # Iteratively optimize thickness to meet atom limit
    # This is critical for DFT efficiency
    current_thickness_a = slab_thickness_a
    current_thickness_b = slab_thickness_b
    min_thickness = 2.5  # Minimum 2.5 Angstroms (very thin, but physically reasonable)
    max_iterations = 10
    
    for iteration in range(max_iterations):
        n_atoms_current = len(slabA_super) + len(slabB_super)
        
        if n_atoms_current <= max_atoms:
            break
        
        print(f"\nIteration {iteration + 1}: Current atom count ({n_atoms_current}) exceeds limit ({max_atoms})")
        
        # Calculate aggressive reduction factor
        # Use 0.6-0.7 factor to ensure we get well under limit
        reduction_factor = (max_atoms / n_atoms_current) * 0.65  # 65% to be safe
        
        # Reduce thickness proportionally
        new_thickness_a = max(current_thickness_a * reduction_factor, min_thickness)
        new_thickness_b = max(current_thickness_b * reduction_factor, min_thickness)
        
        # Don't reduce if already at minimum
        if new_thickness_a >= current_thickness_a - 0.1 and new_thickness_b >= current_thickness_b - 0.1:
            print("Warning: Cannot reduce thickness further. Atom count may exceed limit.")
            break
        
        print(f"Reducing thickness: A: {current_thickness_a:.2f} -> {new_thickness_a:.2f} Å")
        print(f"                   B: {current_thickness_b:.2f} -> {new_thickness_b:.2f} Å")
        
        current_thickness_a = new_thickness_a
        current_thickness_b = new_thickness_b
        
        # Rebuild slabs with reduced thickness
        slabA = build_slab(Aprim, miller_a, current_thickness_a, vacuum)
        slabB = build_slab(Bprim, miller_b, current_thickness_b, vacuum)
        slabA_super = SupercellTransformation(transA).apply_transformation(slabA)
        slabB_super = SupercellTransformation(transB).apply_transformation(slabB)
        
        # Reapply strain
        aA_pre = compute_inplane_norms(slabA_super.lattice)[0:2]
        aB_pre = compute_inplane_norms(slabB_super.lattice)[0:2]
        
        if strain_target.lower() == 'a':
            scale_x = aB_pre[0] / aA_pre[0]
            scale_y = aB_pre[1] / aA_pre[1]
            slabA_strained = apply_inplane_strain(slabA_super, scale_x, scale_y)
            slabB_strained = slabB_super.copy()
        elif strain_target.lower() == 'b':
            scale_x = aA_pre[0] / aB_pre[0]
            scale_y = aA_pre[1] / aB_pre[1]
            slabB_strained = apply_inplane_strain(slabB_super, scale_x, scale_y)
            slabA_strained = slabA_super.copy()
        elif strain_target.lower() == 'both':
            scale_x = math.sqrt(aB_pre[0] / aA_pre[0])
            scale_y = math.sqrt(aB_pre[1] / aA_pre[1])
            slabA_strained = apply_inplane_strain(slabA_super, scale_x, scale_y)
            slabB_strained = apply_inplane_strain(slabB_super, 1.0/scale_x, 1.0/scale_y)
        
        # Update for next iteration
        slabA_super = slabA_strained if strain_target.lower() != 'b' else slabA_super
        slabB_super = slabB_strained if strain_target.lower() != 'a' else slabB_super
    
    final_atoms_before_stack = len(slabA_strained) + len(slabB_strained)
    print(f"\nFinal atom count before stacking: {final_atoms_before_stack} (target: <={max_atoms})")
    
    # align and stack with ordered arrangement (put B as bottom, A as top)
    bottom, top, combined = align_and_stack_ordered(slabB_strained, slabA_strained, separation=sep, vacuum=vacuum)
    
    # Report final atom count
    final_atoms = len(combined)
    print(f"\nFinal structure contains {final_atoms} atoms (limit: {max_atoms})")
    if final_atoms > max_atoms:
        print(f"⚠ WARNING: Atom count ({final_atoms}) exceeds limit ({max_atoms})!")
        print(f"  Consider reducing max_area or using thinner initial slabs.")
    else:
        print(f"✓ Atom count is within limit ({(max_atoms - final_atoms) / max_atoms * 100:.1f}% under limit)")
    
    # Verify final structure has matching in-plane lattice
    final_bottom_lat = bottom.lattice.matrix
    final_top_lat = top.lattice.matrix
    final_combined_lat = combined.lattice.matrix
    
    print("\nFinal lattice verification:")
    print("Bottom in-plane vectors:", final_bottom_lat[0], final_bottom_lat[1])
    print("Top in-plane vectors:", final_top_lat[0], final_top_lat[1])
    print("Combined in-plane vectors:", final_combined_lat[0], final_combined_lat[1])
    
    write_poscars(bottom, top, combined, prefix="auto_interface")

    # print final in-plane match results
    final_a = np.linalg.norm(bottom.lattice.matrix[0]), np.linalg.norm(bottom.lattice.matrix[1])
    final_b = np.linalg.norm(top.lattice.matrix[0]), np.linalg.norm(top.lattice.matrix[1])
    print("Final in-plane lattice sizes (bottom):", final_a)
    print("Final in-plane lattice sizes (top):", final_b)
    # estimate percent strains applied to original prims (informational)
    # NOTE: more careful calculation can be added

if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Auto build matched heterojunction interface from two bulks.")
    p.add_argument("--a", dest="afile", required=True, help="File for material A (CIF/POSCAR)")
    p.add_argument("--b", dest="bfile", required=True, help="File for material B (CIF/POSCAR)")
    p.add_argument("--miller_a", dest="miller_a", default="0,0,1", help="Miller for A, e.g. 0,0,1")
    p.add_argument("--miller_b", dest="miller_b", default="1,0,1", help="Miller for B")
    p.add_argument("--slab_thickness_a", dest="stka", type=float, default=20.0, help="Slab thickness for A (Å)")
    p.add_argument("--slab_thickness_b", dest="stkb", type=float, default=12.0, help="Slab thickness for B (Å) or trilayer estimate")
    p.add_argument("--vacuum", dest="vac", type=float, default=20.0, help="Vacuum size (Å)")
    p.add_argument("--sep", dest="sep", type=float, default=3.2, help="Initial separation between surfaces (Å)")
    p.add_argument("--tol", dest="tol", type=float, default=0.03, help="Matching strain tolerance (fraction, e.g. 0.03)")
    p.add_argument("--max_area", dest="max_area", type=float, default=800.0, help="Maximum interface area (Å^2) for searching matches")
    p.add_argument("--strain_target", dest="strain_target", default="A", choices=["A","B","both"], help="Which material to strain: A, B, or both (split)")
    p.add_argument("--use_builder_interface", dest="use_builder", action="store_true", help="Use CoherentInterfaceBuilder.get_interfaces for ordered interface (recommended)")
    p.add_argument("--max_atoms", dest="max_atoms", type=int, default=400, help="Maximum number of atoms in final structure (default: 400)")
    args = p.parse_args()

    miller_a = tuple(map(int, args.miller_a.split(',')))
    miller_b = tuple(map(int, args.miller_b.split(',')))

    auto_run(args.afile, args.bfile, miller_a, miller_b,
             args.stka, args.stkb, args.vac, args.sep,
             args.tol, args.max_area, args.strain_target, args.use_builder, args.max_atoms)
