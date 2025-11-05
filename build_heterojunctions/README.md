# Auto Interface Builder - Usage Guide

This guide provides examples of how to use `auto_interface_builder.py` to build heterojunction interfaces between different materials.

## Basic Usage

```bash
python build_heterojunctions/auto_interface_builder.py \
    --a <structure_A_file> \
    --b <structure_B_file> \
    --miller_a <h,k,l> \
    --miller_b <h,k,l>
```

## Required Arguments

- `--a`: Path to structure file A (CIF/POSCAR format)
- `--b`: Path to structure file B (CIF/POSCAR format)
- `--miller_a`: Miller indices for material A (e.g., `0,0,1`)
- `--miller_b`: Miller indices for material B (e.g., `1,0,1`)

## Optional Arguments

- `--slab_thickness_a`: Slab thickness for A in Å (default: 20.0)
- `--slab_thickness_b`: Slab thickness for B in Å (default: 12.0)
- `--vacuum`: Vacuum size in Å (default: 20.0)
- `--sep`: Initial separation between surfaces in Å (default: 3.2)
- `--tol`: Matching strain tolerance (default: 0.03)
- `--max_area`: Maximum interface area in Å² for searching matches (default: 800.0)
- `--strain_target`: Which material to strain: `A`, `B`, or `both` (default: `A`)
- `--use_builder_interface`: Use CoherentInterfaceBuilder method (recommended for ordered interfaces)
- `--max_atoms`: Maximum number of atoms in final structure (default: 400)

## Example 1: Perovskite/TiO2 Interface (Basic)

```bash
python build_heterojunctions/auto_interface_builder.py \
    --a structures/perovskites/H6PbCI3N.cif \
    --b structures/etl/TiO2.cif \
    --miller_a 0,0,1 \
    --miller_b 1,0,1 \
    --vacuum 20 \
    --sep 3.2 \
    --strain_target A
```

## Example 2: Perovskite/TiO2 with Atom Limit Optimization

Limit the total number of atoms to 300 for efficient DFT calculations:

```bash
python build_heterojunctions/auto_interface_builder.py \
    --a structures/perovskites/H6PbCI3N.cif \
    --b structures/etl/TiO2.cif \
    --miller_a 0,0,1 \
    --miller_b 1,0,1 \
    --vacuum 20 \
    --sep 3.2 \
    --strain_target A \
    --max_atoms 300 \
    --use_builder_interface
```

## Example 3: Different Miller Indices

Try different surface orientations for better lattice matching:

```bash
# Using (1,0,0) surfaces
python build_heterojunctions/auto_interface_builder.py \
    --a structures/perovskites/H6PbCI3N.cif \
    --b structures/etl/TiO2.cif \
    --miller_a 1,0,0 \
    --miller_b 1,0,0 \
    --vacuum 20 \
    --sep 3.2 \
    --strain_target A \
    --max_atoms 300

# Using (1,1,0) surfaces
python build_heterojunctions/auto_interface_builder.py \
    --a structures/perovskites/H6PbCI3N.cif \
    --b structures/etl/TiO2.cif \
    --miller_a 1,1,0 \
    --miller_b 1,1,0 \
    --vacuum 20 \
    --sep 3.2 \
    --strain_target both \
    --max_atoms 300
```

## Example 4: Thinner Slabs for Smaller Systems

Use thinner slabs to reduce atom count:

```bash
python build_heterojunctions/auto_interface_builder.py \
    --a structures/perovskites/H6PbCI3N.cif \
    --b structures/etl/TiO2.cif \
    --miller_a 0,0,1 \
    --miller_b 1,0,1 \
    --slab_thickness_a 10.0 \
    --slab_thickness_b 6.0 \
    --vacuum 20 \
    --sep 3.2 \
    --strain_target A \
    --max_atoms 250
```

## Example 5: Larger Tolerance for More Matches

Increase tolerance to find more matching possibilities:

```bash
python build_heterojunctions/auto_interface_builder.py \
    --a structures/perovskites/H6PbCI3N.cif \
    --b structures/etl/TiO2.cif \
    --miller_a 0,0,1 \
    --miller_b 1,0,1 \
    --tol 0.05 \
    --max_area 1500 \
    --vacuum 20 \
    --sep 3.2 \
    --strain_target A \
    --max_atoms 300
```

## Example 6: Split Strain Between Both Materials

Apply strain to both materials evenly:

```bash
python build_heterojunctions/auto_interface_builder.py \
    --a structures/perovskites/H6PbCI3N.cif \
    --b structures/etl/TiO2.cif \
    --miller_a 0,0,1 \
    --miller_b 1,0,1 \
    --vacuum 20 \
    --sep 3.2 \
    --strain_target both \
    --max_atoms 300 \
    --use_builder_interface
```

## Example 7: Different Material Combinations

### MAPbI3 / SnO2 Interface
```bash
python build_heterojunctions/auto_interface_builder.py \
    --a structures/perovskites/MAPbI3.cif \
    --b structures/etl/SnO2.cif \
    --miller_a 0,0,1 \
    --miller_b 1,0,0 \
    --vacuum 20 \
    --sep 3.5 \
    --strain_target A \
    --max_atoms 350
```

### Perovskite / Perovskite Interface
```bash
python build_heterojunctions/auto_interface_builder.py \
    --a structures/perovskites/CsPbI3.cif \
    --b structures/perovskites/MAPbI3.cif \
    --miller_a 0,0,1 \
    --miller_b 0,0,1 \
    --vacuum 20 \
    --sep 3.0 \
    --strain_target both \
    --max_atoms 400
```

## Example 8: Very Small System for Quick Testing

Build a minimal system for rapid testing:

```bash
python build_heterojunctions/auto_interface_builder.py \
    --a structures/perovskites/H6PbCI3N.cif \
    --b structures/etl/TiO2.cif \
    --miller_a 0,0,1 \
    --miller_b 1,0,1 \
    --slab_thickness_a 5.0 \
    --slab_thickness_b 3.0 \
    --vacuum 15 \
    --sep 3.0 \
    --max_area 500 \
    --tol 0.04 \
    --strain_target A \
    --max_atoms 200 \
    --use_builder_interface
```

## Example 9: Large Area Interface

Build a larger interface for better statistics:

```bash
python build_heterojunctions/auto_interface_builder.py \
    --a structures/perovskites/H6PbCI3N.cif \
    --b structures/etl/TiO2.cif \
    --miller_a 0,0,1 \
    --miller_b 1,0,1 \
    --slab_thickness_a 15.0 \
    --slab_thickness_b 10.0 \
    --vacuum 25 \
    --sep 3.5 \
    --max_area 2000 \
    --tol 0.03 \
    --strain_target A \
    --max_atoms 500
```

## Output Files

The script generates three POSCAR files:

1. `auto_interface_bottom_strained.vasp` - Bottom slab (substrate) structure
2. `auto_interface_top_strained.vasp` - Top slab (film) structure  
3. `auto_interface_combined.vasp` - Combined heterojunction structure

## Tips for Optimization

1. **Atom Count**: Use `--max_atoms` to control system size. The script will automatically reduce slab thickness if needed.

2. **Miller Indices**: Try different Miller indices combinations to find better lattice matches with smaller supercells.

3. **Tolerance**: Increase `--tol` to find more matches, but this may increase strain.

4. **Builder Interface**: Use `--use_builder_interface` for better ordered interfaces with high symmetry.

5. **Strain Target**: 
   - `A`: Strain material A to match B (good when A is more flexible)
   - `B`: Strain material B to match A (good when B is more flexible)
   - `both`: Split strain evenly (good for balanced systems)

6. **Max Area**: Reducing `--max_area` helps find smaller supercells, reducing atom count.

## Troubleshooting

- **No matches found**: Increase `--tol` or `--max_area`
- **Too many atoms**: Reduce `--max_atoms`, `--slab_thickness_a/b`, or `--max_area`
- **Poor matching**: Try different `--miller_a` and `--miller_b` combinations
- **Interface quality**: Use `--use_builder_interface` for better ordered structures

## Next Steps

After generating the interface, use `fix_interface_layers.py` to apply selective dynamics for relaxation:

```bash
python build_heterojunctions/fix_interface_layers.py \
    --poscar auto_interface_combined.vasp \
    --n_layers 2 \
    --output auto_interface_combined_relaxed.vasp
```

