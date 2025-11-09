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
from pathlib import Path
from typing import Dict, Iterable, List, Set, Tuple

import numpy as np
from pymatgen.core import Lattice, Structure
from pymatgen.core.surface import SlabGenerator
from pymatgen.transformations.standard_transformations import SupercellTransformation
# from pymatgen.analysis.interfaces import InterfaceMatcher
from pymatgen.analysis.interfaces.coherent_interfaces import CoherentInterfaceBuilder, Interface
from pymatgen.io.vasp import Poscar


def _extend_matrix_2d_to_3d(mat_2d: np.ndarray) -> np.ndarray:
    """
    Promote a 2x2 transformation matrix to a 3x3 supercell matrix by keeping the
    original in-plane mapping and leaving the out-of-plane direction unchanged.
    """
    mat_3d = np.eye(3)
    mat_3d[0:2, 0:2] = np.array(mat_2d)
    return np.rint(mat_3d).astype(int)


def _matrix_key(film_trans: np.ndarray, sub_trans: np.ndarray) -> Tuple[int, ...]:
    """Create a hashable key describing a pair of 2x2 transformation matrices."""
    film_int = tuple(np.rint(film_trans).astype(int).flatten().tolist())
    sub_int = tuple(np.rint(sub_trans).astype(int).flatten().tolist())
    return film_int + sub_int


def _angle_between(vec_a: np.ndarray, vec_b: np.ndarray) -> float:
    """Return the angle between two vectors in degrees."""
    dot_val = np.dot(vec_a, vec_b)
    norms = np.linalg.norm(vec_a) * np.linalg.norm(vec_b)
    if norms < 1e-8:
        return 0.0
    cos_theta = np.clip(dot_val / norms, -1.0, 1.0)
    return float(np.degrees(np.arccos(cos_theta)))


def _compute_match_metrics(match_obj) -> Dict[str, np.ndarray | float]:
    """Extract length/angle mismatch information from a ZSLMatch."""
    film_vectors = np.array(match_obj.film_vectors)
    sub_vectors = np.array(match_obj.substrate_vectors)

    film_lengths = np.linalg.norm(film_vectors, axis=1)
    sub_lengths = np.linalg.norm(sub_vectors, axis=1)
    avg_lengths = (film_lengths + sub_lengths) / 2.0
    # Relative mismatch per in-plane vector (signed)
    length_components = np.divide(
        (film_lengths - sub_lengths),
        np.where(avg_lengths < 1e-8, 1.0, avg_lengths),
    )

    film_angle = _angle_between(film_vectors[0], film_vectors[1])
    sub_angle = _angle_between(sub_vectors[0], sub_vectors[1])
    angle_mismatch = abs(film_angle - sub_angle)

    return {
        "film_lengths": film_lengths,
        "substrate_lengths": sub_lengths,
        "length_components": length_components,
        "max_length_mismatch": float(np.max(np.abs(length_components))),
        "film_angle": film_angle,
        "substrate_angle": sub_angle,
        "angle_mismatch": angle_mismatch,
    }


def _generate_zsl_parameter_grid(
    base_tol: float,
    max_area: float,
    allow_bidirectional: bool = True,
) -> List[Dict[str, float | bool]]:
    """
    Generate a prioritized list of parameter dictionaries for ZSLGenerator searches.
    Parameters with smaller area, tighter tolerances, and bidirectional matching are
    considered first to encourage compact supercells.
    """
    base_tol = max(0.01, min(base_tol, 0.15))
    length_candidates = {
        base_tol,
        min(0.15, base_tol * 1.35),
        min(0.15, base_tol * 1.7),
        min(0.15, base_tol + 0.02),
        min(0.15, 0.12 if base_tol < 0.12 else base_tol),
    }
    length_values = sorted(length_candidates)

    area_candidates = {
        max_area,
        max_area * 1.3,
        max_area * 1.7,
        max_area * 2.2,
    }
    area_values = sorted({float(max(10.0, a)) for a in area_candidates})

    bidi_options: Iterable[bool] = [True, False] if allow_bidirectional else [False]

    param_grid: List[Dict[str, float | bool]] = []
    for length_tol in length_values:
        angle_tol = min(0.1, max(0.015, length_tol / 3.0))
        for area in area_values:
            for bidi in bidi_options:
                param_grid.append(
                    {
                        "max_area": float(area),
                        "max_length_tol": float(length_tol),
                        "max_angle_tol": float(angle_tol),
                        "bidirectional": bool(bidi),
                    }
                )

    param_grid.sort(key=lambda p: (p["max_area"], p["max_length_tol"], 0 if p["bidirectional"] else 1))
    return param_grid


def _safe_structure_label(path_str: str) -> str:
    """
    Derive a readable label from an input structure path by stripping directories,
    removing the file suffix, and replacing spaces with underscores.
    """
    path = Path(path_str)
    label = path.stem  # Removes the final suffix such as .cif, .vasp, etc.
    return label.replace(" ", "_")


def _reorder_structure_for_poscar(struct: Structure) -> Structure:
    """
    Group sites by species so that POSCAR headers list each element once in the
    order they first appear. Within each species block, maintain ascending
    fractional z for readability.
    """
    struct_copy = struct.copy()
    species_order: List[str] = []
    for site in struct_copy:
        symbol = site.species_string
        if symbol not in species_order:
            species_order.append(symbol)
    order_index = {symbol: idx for idx, symbol in enumerate(species_order)}
    struct_copy.sort(key=lambda site: (order_index[site.species_string], site.frac_coords[2]))
    return struct_copy

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

def search_matches(
    structA,
    structB,
    miller_a,
    miller_b,
    tol: float = 0.03,
    max_area: float = 500,
    max_match: int = 20,
    allow_bidirectional: bool = True,
):
    """
    Use CoherentInterfaceBuilder to find candidate in-plane supercells that match.
    Requires bulk structures and Miller indices, not slabs.
    Returns a list of match dicts with strains and transforms.
    """
    from pymatgen.analysis.interfaces.zsl import ZSLGenerator

    param_grid = _generate_zsl_parameter_grid(tol, max_area, allow_bidirectional=allow_bidirectional)
    matches: List[dict] = []
    seen_keys: Set[Tuple[int, ...]] = set()

    for params in param_grid:
        try:
            zslgen = ZSLGenerator(
                max_area=params["max_area"],
                max_length_tol=params["max_length_tol"],
                max_angle_tol=params["max_angle_tol"],
                bidirectional=params["bidirectional"],
            )
        except TypeError:
            # Legacy pymatgen versions may not support all keyword arguments.
            zslgen = ZSLGenerator(max_area=params["max_area"])

        matcher = CoherentInterfaceBuilder(
            substrate_structure=structB,
            film_structure=structA,
            film_miller=miller_a,
            substrate_miller=miller_b,
            zslgen=zslgen,
        )

        zsl_matches = matcher.zsl_matches or []
        if len(zsl_matches) == 0:
            continue

        for idx, match_obj in enumerate(zsl_matches[: max_match or None]):
            film_trans = np.array(match_obj.film_transformation)
            sub_trans = (
                np.array(match_obj.substrate_transformation)
                if hasattr(match_obj, "substrate_transformation")
                else np.eye(2)
            )
            key = _matrix_key(film_trans, sub_trans)
            if key in seen_keys:
                continue

            metrics = _compute_match_metrics(match_obj)
            summary = {
                "area": float(match_obj.match_area),
                "film_transformation": film_trans,
                "substrate_transformation": sub_trans,
                "supercell_a": _extend_matrix_2d_to_3d(film_trans),
                "supercell_b": _extend_matrix_2d_to_3d(sub_trans),
                "det_a": float(abs(np.linalg.det(film_trans))),
                "det_b": float(abs(np.linalg.det(sub_trans))),
                "length_components": metrics["length_components"],
                "length_mismatch": metrics["max_length_mismatch"],
                "angle_mismatch": metrics["angle_mismatch"],
                "angles": (metrics["film_angle"], metrics["substrate_angle"]),
                "film_lengths": metrics["film_lengths"],
                "substrate_lengths": metrics["substrate_lengths"],
                "zsl_params": params,
                "priority": (
                    params["max_area"],
                    params["max_length_tol"],
                    0 if params["bidirectional"] else 1,
                    idx,
                ),
                "match_index": idx,
            }
            matches.append(summary)
            seen_keys.add(key)

    matches.sort(key=lambda m: m["priority"])
    if max_match:
        matches = matches[:max_match]
    return matches

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
    Returns bottom, top, and combined structures, where bottom and top are properly separated.
    """
    # Copy slabs to avoid modifying originals
    bottom = slab_bottom.copy()
    top = slab_top.copy()
    
    # Use bottom's lattice as the reference (it should match top after straining)
    # Get in-plane vectors from bottom to determine the interface normal.
    bottom_lat = bottom.lattice.matrix
    a_vec = bottom_lat[0]
    b_vec = bottom_lat[1]

    # Calculate surface normal from in-plane vectors
    ab_normal = np.cross(a_vec, b_vec)
    ab_normal_norm = np.linalg.norm(ab_normal)
    if ab_normal_norm < 1e-10:
        raise ValueError("In-plane vectors are parallel or degenerate")
    ab_normal_unit = ab_normal / ab_normal_norm
    # Ensure the surface normal points in the same general direction as the original c vector.
    if np.dot(ab_normal_unit, bottom_lat[2]) < 0:
        ab_normal_unit *= -1.0

    # Align bottom slab along the interface normal so the minimum projection is 0.
    if len(bottom) > 0:
        bottom_proj = np.dot(bottom.cart_coords, ab_normal_unit)
        min_proj_bottom = bottom_proj.min()
        bottom.translate_sites(
            list(range(len(bottom))),
            -ab_normal_unit * min_proj_bottom,
        )
        bottom_proj = np.dot(bottom.cart_coords, ab_normal_unit)
        max_proj_bottom = bottom_proj.max()
    else:
        max_proj_bottom = 0.0

    # Align top slab so its minimum projection is at the desired separation.
    if len(top) > 0:
        top_proj = np.dot(top.cart_coords, ab_normal_unit)
        min_proj_top = top_proj.min()
        top.translate_sites(
            list(range(len(top))),
            -ab_normal_unit * min_proj_top,
        )
        top.translate_sites(
            list(range(len(top))),
            ab_normal_unit * (max_proj_bottom + separation),
        )
        top_proj = np.dot(top.cart_coords, ab_normal_unit)
        desired_min = max_proj_bottom + separation + 0.1
        min_proj_top_post = top_proj.min()
        if min_proj_top_post < desired_min - 1e-4:
            adjust = desired_min - min_proj_top_post
            top.translate_sites(list(range(len(top))), ab_normal_unit * adjust)
            top_proj = np.dot(top.cart_coords, ab_normal_unit)
        max_proj_top = top_proj.max()
    else:
        max_proj_top = max_proj_bottom + separation

    # Determine combined cell length along the interface normal and construct the new lattice.
    total_z = max_proj_top + vacuum
    
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
    # Create c vector perpendicular to ab plane with correct magnitude
    c_vec = ab_normal_unit * total_z
    
    # Build combined lattice using bottom's in-plane vectors and perpendicular c
    combined_lat_matrix = np.vstack([a_vec, b_vec, c_vec])
    combined_lattice = Lattice(combined_lat_matrix)
    
    # Create structures with the combined lattice for bottom and top
    # This ensures they use the same lattice as the combined structure
    # Convert bottom slab to use combined lattice
    bottom_cart_coords = bottom.cart_coords
    bottom_frac_coords = [combined_lattice.get_fractional_coords(coord) for coord in bottom_cart_coords]
    # Wrap fractional coordinates to [0, 1)
    bottom_frac_coords = [(fc % 1.0) for fc in bottom_frac_coords]
    bottom_separated = Structure(combined_lattice,
                                bottom.species,
                                bottom_frac_coords,
                                coords_are_cartesian=False)
    
    # Convert top slab to use combined lattice
    top_cart_coords = top.cart_coords
    top_frac_coords = [combined_lattice.get_fractional_coords(coord) for coord in top_cart_coords]
    # Wrap fractional coordinates to [0, 1)
    top_frac_coords = [(fc % 1.0) for fc in top_frac_coords]
    top_separated = Structure(combined_lattice,
                            top.species,
                            top_frac_coords,
                            coords_are_cartesian=False)
    
    # Create combined structure by adding all sites
    combined = Structure(combined_lattice, [], [])
    
    # Add bottom slab sites first
    for site in bottom_separated:
        combined.append(site.species_string, site.frac_coords, coords_are_cartesian=False)
    
    # Add top slab sites
    for site in top_separated:
        combined.append(site.species_string, site.frac_coords, coords_are_cartesian=False)
    
    # Sort combined by z-coordinate for better ordering
    combined.sort(key=lambda site: site.frac_coords[2])
    
    # Verify separation: ensure bottom and top are properly separated
    if len(bottom_separated) > 0 and len(top_separated) > 0:
        c_unit = ab_normal_unit
        bottom_proj = np.dot(bottom_separated.cart_coords, c_unit)
        top_proj = np.dot(top_separated.cart_coords, c_unit)
        bottom_proj_max = bottom_proj.max()
        top_proj_min = top_proj.min()
        gap_size = top_proj_min - bottom_proj_max
        if gap_size < 0:
            print(
                f"Warning: Overlap detected along interface normal! "
                f"bottom_max={bottom_proj_max:.2f} Å, top_min={top_proj_min:.2f} Å"
            )
        elif gap_size < 0.5:
            print(f"Warning: Very small gap ({gap_size:.2f} Å) along interface normal.")
    
    return bottom_separated, top_separated, combined

def write_poscars(bottom, top, combined, output_dir: Path, base_name: str):
    output_dir.mkdir(parents=True, exist_ok=True)
    bottom_for_io = _reorder_structure_for_poscar(bottom)
    top_for_io = _reorder_structure_for_poscar(top)
    combined_for_io = _reorder_structure_for_poscar(combined)

    bottom_path = output_dir / f"{base_name}_bottom_strained.vasp"
    top_path = output_dir / f"{base_name}_top_strained.vasp"
    combined_path = output_dir / f"{base_name}.vasp"

    Poscar(bottom_for_io, sort_structure=False).write_file(bottom_path)
    Poscar(top_for_io, sort_structure=False).write_file(top_path)
    Poscar(combined_for_io, sort_structure=False).write_file(combined_path)

    print(
        "Wrote POSCARs:",
        bottom_path,
        top_path,
        combined_path,
    )

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

def build_interface_from_builder(
    structA,
    structB,
    miller_a,
    miller_b,
    slab_thickness_a,
    slab_thickness_b,
    vacuum,
    sep,
    zsl_param_candidates,
    target_match=None,
    max_atoms: int = 400,
):
    """
    Attempt to construct an ordered interface using CoherentInterfaceBuilder. The search
    iterates over candidate ZSL parameter sets and slab thickness scaling factors until
    an interface within the atom limit is located. If none satisfy the limit, the
    smallest interface found is returned.
    """
    from pymatgen.analysis.interfaces.zsl import ZSLGenerator

    if not zsl_param_candidates:
        zsl_param_candidates = [
            {
                "max_area": 800.0,
                "max_length_tol": 0.08,
                "max_angle_tol": 0.03,
                "bidirectional": True,
            }
        ]

    thickness_factors = [1.0, 0.85, 0.7, 0.55, 0.45, 0.35, 0.3]
    min_thickness = 3.0

    target_key = None
    if target_match:
        target_key = _matrix_key(
            np.array(target_match["film_transformation"]),
            np.array(target_match["substrate_transformation"]),
        )

    best_candidate = None
    best_over_limit = None

    for factor in thickness_factors:
        current_thickness_a = max(slab_thickness_a * factor, min_thickness)
        current_thickness_b = max(slab_thickness_b * factor, min_thickness)

        for params in zsl_param_candidates:
            try:
                zslgen = ZSLGenerator(
                    max_area=params.get("max_area", 800.0),
                    max_length_tol=params.get("max_length_tol", 0.08),
                    max_angle_tol=params.get("max_angle_tol", 0.03),
                    bidirectional=params.get("bidirectional", True),
                )
            except TypeError:
                zslgen = ZSLGenerator(max_area=params.get("max_area", 800.0))

            builder = CoherentInterfaceBuilder(
                substrate_structure=structB,
                film_structure=structA,
                film_miller=miller_a,
                substrate_miller=miller_b,
                zslgen=zslgen,
            )

            if not builder.zsl_matches:
                continue

            terminations = builder.terminations
            if not terminations:
                continue

            matches = builder.zsl_matches
            for termination in terminations:
                try:
                    interfaces = list(
                        builder.get_interfaces(
                            termination=termination,
                            gap=sep,
                            vacuum_over_film=vacuum,
                            film_thickness=current_thickness_a,
                            substrate_thickness=current_thickness_b,
                            in_layers=False,
                        )
                    )
                except Exception:
                    continue

                for idx, interface in enumerate(interfaces):
                    if idx >= len(matches):
                        break

                    match_obj = matches[idx]
                    match_key = _matrix_key(
                        np.array(match_obj.film_transformation),
                        np.array(match_obj.substrate_transformation),
                    )
                    if target_key and match_key != target_key:
                        continue

                    n_atoms = len(interface)
                    candidate_info = {
                        "interface": interface,
                        "termination": termination,
                        "params": params,
                        "thickness_a": current_thickness_a,
                        "thickness_b": current_thickness_b,
                        "atoms": n_atoms,
                        "match_key": match_key,
                    }

                    if n_atoms <= max_atoms:
                        if best_candidate is None or n_atoms < best_candidate["atoms"]:
                            best_candidate = candidate_info
                    else:
                        if best_over_limit is None or n_atoms < best_over_limit["atoms"]:
                            best_over_limit = candidate_info

        if best_candidate is not None:
            break

    chosen = best_candidate or best_over_limit
    if chosen is None:
        return None, None, None

    interface: "Interface" = chosen["interface"]
    print(f"Using termination: {chosen['termination']}")
    print(f"Interface contains {chosen['atoms']} atoms (target: <= {max_atoms})")
    if chosen["atoms"] > max_atoms:
        print(f"Warning: Atom count ({chosen['atoms']}) exceeds limit ({max_atoms})")
    print(
        f"ZSL parameters: {chosen['params']} | film thickness={chosen['thickness_a']:.2f} Å, "
        f"substrate thickness={chosen['thickness_b']:.2f} Å"
    )

    combined = interface.copy()
    combined.sort(key=lambda site: site.frac_coords[2])

    bottom = Structure.from_sites(interface.substrate_sites)
    top = Structure.from_sites(interface.film_sites)

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

    a_label = _safe_structure_label(a_file)
    b_label = _safe_structure_label(b_file)
    hetero_base_name = f"{a_label}@{b_label}"
    hetero_output_dir = Path("structures") / "heterojunctions"

    # Search for matches using bulk structures and Miller indices
    print("Searching for low-strain matches (this may take a moment)...")
    matches = search_matches(
        Aprim,
        Bprim,
        miller_a,
        miller_b,
        tol=tol,
        max_area=max_area,
        allow_bidirectional=True,
    )
    
    if len(matches) == 0:
        raise RuntimeError("No matches found within tolerance. Increase tol or max_area.")
    
    # Collect unique ZSL generator parameter sets observed in the matches.
    zsl_param_list: List[Dict[str, float | bool]] = []
    param_keys_seen: Set[Tuple[Tuple[str, float | bool], ...]] = set()
    for match in matches:
        params = match.get("zsl_params")
        if not params:
            continue
        key = tuple(sorted(params.items()))
        if key in param_keys_seen:
            continue
        param_keys_seen.add(key)
        zsl_param_list.append(params)
    
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
    
    # Score matches prioritizing atom count, then strain/area metrics
    def match_score(m):
        opt_thick_a, opt_thick_b = estimate_optimal_thickness(m, max_atoms)
        est_atoms = estimate_atoms_for_match(
            Aprim,
            Bprim,
            miller_a,
            miller_b,
            m["supercell_a"],
            m["supercell_b"],
            opt_thick_a,
            opt_thick_b,
        )

        length_metric = float(m.get("length_mismatch", 0.0))
        angle_metric = float(m.get("angle_mismatch", 0.0))
        area_metric = float(m.get("area", 0.0))
        supercell_factor = float(max(m.get("det_a", 0.0), m.get("det_b", 0.0)))

        if est_atoms > max_atoms:
            atom_penalty = (est_atoms - max_atoms) * 4.0 + (est_atoms / max_atoms - 1.0) * 20.0
        else:
            atom_penalty = -((max_atoms - est_atoms) / max_atoms) * 5.0

        strain_penalty = length_metric * 200.0 + angle_metric * 2.0
        area_penalty = (area_metric / max(max_area, 1.0)) * 5.0
        supercell_penalty = supercell_factor

        total_penalty = atom_penalty + strain_penalty + area_penalty + supercell_penalty

        m["_cached_opt_thickness"] = (opt_thick_a, opt_thick_b)
        m["_cached_est_atoms"] = est_atoms
        m["_cached_strain_penalty"] = strain_penalty

        return total_penalty, est_atoms, strain_penalty

    matches_with_scores = [(idx, m, match_score(m)) for idx, m in enumerate(matches)]

    matches_with_scores.sort(
        key=lambda item: (
            0 if item[2][1] <= max_atoms else 1,
            item[2][1],
            item[1].get("priority", (item[0],)),
            item[2][0],
        )
    )

    valid_matches = [item for item in matches_with_scores if item[2][1] <= max_atoms]

    if len(valid_matches) > 0:
        ordered_items = valid_matches
        print(f"Found {len(valid_matches)} matches within atom limit ({max_atoms})")
    else:
        ordered_items = matches_with_scores
        print(f"Warning: No matches within atom limit ({max_atoms}). Selecting closest candidate...")

    best_item = ordered_items[0]
    best = best_item[1]
    (best_score, best_atoms, best_strain_metric) = best_item[2]
    
    # Calculate optimal thickness for selected match
    optimal_thick_a, optimal_thick_b = best.get("_cached_opt_thickness", estimate_optimal_thickness(best, max_atoms))
    
    print(f"Selected match with estimated {best_atoms} atoms (strain metric: {best_strain_metric:.4f})")
    print(f"Optimal thickness for this match: A={optimal_thick_a:.2f} Å, B={optimal_thick_b:.2f} Å")
    
    # Use optimal thickness if significantly different
    if abs(optimal_thick_a - slab_thickness_a) > 0.5 or abs(optimal_thick_b - slab_thickness_b) > 0.5:
        print(f"Adjusting thickness to meet atom limit:")
        print(f"  A: {slab_thickness_a:.2f} -> {optimal_thick_a:.2f} Å")
        print(f"  B: {slab_thickness_b:.2f} -> {optimal_thick_b:.2f} Å")
        slab_thickness_a = optimal_thick_a
        slab_thickness_b = optimal_thick_b

    if use_builder_interface:
        print("Attempting ordered interface construction via CoherentInterfaceBuilder...")
        builder_param_candidates = zsl_param_list or [
            {
                "max_area": float(max_area),
                "max_length_tol": float(max(tol, 0.03)),
                "max_angle_tol": float(min(0.08, max(tol / 2.0, 0.02))),
                "bidirectional": True,
            }
        ]

        bottom, top, combined = build_interface_from_builder(
            Aprim,
            Bprim,
            miller_a,
            miller_b,
            slab_thickness_a,
            slab_thickness_b,
            vacuum,
            sep,
            builder_param_candidates,
            target_match=best,
            max_atoms=max_atoms,
        )

        if combined is not None:
            print("Successfully built ordered interface using CoherentInterfaceBuilder.")
            write_poscars(bottom, top, combined, hetero_output_dir, hetero_base_name)
            return
        else:
            print("Ordered interface generation failed; falling back to manual stacking workflow.")
    
    print("Best match summary:")
    print(" area:", best["area"])
    print(" supercell A matrix:\n", best["supercell_a"])
    print(" supercell B matrix:\n", best["supercell_b"])
    if "det_a" in best and "det_b" in best:
        print(f" determinants det(A)={best['det_a']:.2f}, det(B)={best['det_b']:.2f}")
    if "length_components" in best:
        print(
            " length mismatch (fractional per in-plane vector):",
            np.array2string(np.array(best["length_components"]), precision=5),
        )
    if "angle_mismatch" in best:
        angles = best.get("angles", (None, None))
        if angles[0] is not None and angles[1] is not None:
            print(
                f" angles (film/substrate): {angles[0]:.3f}° / {angles[1]:.3f}° | mismatch {best['angle_mismatch']:.3f}°"
            )
        else:
            print(f" angle mismatch: {best['angle_mismatch']:.3f}°")
    if "zsl_params" in best:
        print(" ZSL parameter set:", best["zsl_params"])
    # Use the InterfaceMatch data to build the supercell structures manually.

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
    
    write_poscars(bottom, top, combined, hetero_output_dir, hetero_base_name)

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
