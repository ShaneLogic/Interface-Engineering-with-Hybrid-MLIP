# Auto Interface Builder - Usage Guide

This guide provides examples of how to use `auto_interface_builder.py` and `fix_interface_layers.py` to build heterojunction interfaces with selective dynamics for DFT relaxation calculations.

## Overview

The workflow consists of two steps:
1. **Generate interface structure** using `auto_interface_builder.py`
2. **Add selective dynamics** using `fix_interface_layers.py` to relax interface layers and fix bulk layers

## Step 1: Generate Interface Structure

### Basic Usage

```bash
python build_heterojunctions/auto_interface_builder.py \
    --a structures/perovskites/H6PbCI3N.cif \
    --b structures/etl/TiO2.cif \
    --miller_a 0,0,1 \
    --miller_b 1,0,1 \
    --vacuum 20 \
    --sep 3.2 \
    --strain_target A \
    --use_builder_interface \
    --max_atoms 300
```

### Required Arguments

- `--a`: Path to structure file A (CIF/POSCAR format)
- `--b`: Path to structure file B (CIF/POSCAR format)
- `--miller_a`: Miller indices for material A (e.g., `0,0,1`)
- `--miller_b`: Miller indices for material B (e.g., `1,0,1`)

### Optional Arguments

- `--slab_thickness_a`: Slab thickness for A in Å (default: 20.0)
- `--slab_thickness_b`: Slab thickness for B in Å (default: 12.0)
- `--vacuum`: Vacuum size in Å (default: 20.0)
- `--sep`: Initial separation between surfaces in Å (default: 3.2)
- `--tol`: Matching strain tolerance (default: 0.03)
- `--max_area`: Maximum interface area in Å² for searching matches (default: 800.0)
- `--strain_target`: Which material to strain: `A`, `B`, or `both` (default: `A`)
- `--use_builder_interface`: Use CoherentInterfaceBuilder method (recommended for ordered interfaces)
- `--max_atoms`: Maximum number of atoms in final structure (default: 400)

### Output Files from Step 1

The script generates three POSCAR files:

1. `auto_interface_bottom_strained.vasp` - Bottom slab (substrate) structure
2. `auto_interface_top_strained.vasp` - Top slab (film) structure  
3. `auto_interface_combined.vasp` - Combined heterojunction structure (input for Step 2)

## Step 2: Add Selective Dynamics

### Basic Usage

```bash
python build_heterojunctions/fix_interface_layers.py \
    auto_interface_combined.vasp \
    -o auto_interface_combined_relaxed.vasp \
    -n 2 \
    -t 2.0
```

### Parameters for fix_interface_layers.py

- `input_file`: Input POSCAR file (combined heterojunction from Step 1)
- `-o, --output`: Output POSCAR file (default: input_file with '_relaxed' suffix)
- `-n, --n_layers`: Number of layers to relax on each side of interface (default: 1)
- `-t, --thickness`: Layer thickness threshold in Angstroms (default: 2.0)
- `-m, --method`: Interface identification method: `density_gap`, `median`, or `max_gap` (default: `density_gap`)

### Selective Dynamics Format

The output POSCAR file contains:

```
Selective dynamics
direct
  0.1234  0.5678  0.9012  T T T  Ti4+
  0.2345  0.6789  0.0123  F F F  O2-
  ...
```

- `T T T`: Atom is free to relax (interface layers)
- `F F F`: Atom is fixed (bulk layers)

## Complete Workflow Examples

### Example 1: Standard Workflow (Recommended)

**Step 1**: Generate interface structure

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

**Step 2**: Add selective dynamics

```bash
python build_heterojunctions/fix_interface_layers.py \
    auto_interface_combined.vasp \
    -o auto_interface_combined_relaxed.vasp \
    -n 2 \
    -t 2.0 \
    -m density_gap
```

This generates `auto_interface_combined_relaxed.vasp` ready for DFT calculations.

### Example 2: Relax More Interface Layers

To relax 3 layers on each side of the interface:

```bash
# Step 1: Generate interface
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

# Step 2: Add selective dynamics with 3 layers
python build_heterojunctions/fix_interface_layers.py \
    auto_interface_combined.vasp \
    -o auto_interface_combined_relaxed.vasp \
    -n 3 \
    -t 2.5
```

### Example 3: Smaller System for Quick Testing

```bash
# Step 1: Generate small interface
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

# Step 2: Add selective dynamics
python build_heterojunctions/fix_interface_layers.py \
    auto_interface_combined.vasp \
    -n 1 \
    -t 2.0
```

### Example 4: Different Strain Strategy

```bash
# Step 1: Generate interface with split strain
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

# Step 2: Add selective dynamics
python build_heterojunctions/fix_interface_layers.py \
    auto_interface_combined.vasp \
    -n 2 \
    -t 2.0
```

### Example 5: Using Different Interface Identification Method

```bash
# Step 1: Generate interface
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

# Step 2: Use median method for interface identification
python build_heterojunctions/fix_interface_layers.py \
    auto_interface_combined.vasp \
    -n 2 \
    -t 2.0 \
    -m median
```

## Understanding Selective Dynamics

The selective dynamics feature:

1. **Identifies interface**: Automatically finds the boundary between two materials using density gap, median, or max gap method
2. **Selects interface layers**: Chooses the N closest atomic layers on each side of the interface based on z-coordinates
3. **Applies constraints**:
   - Interface layers: `T T T` (free to relax in all directions)
   - Bulk layers: `F F F` (fixed positions)

This is ideal for DFT calculations where you want to:
- Relax the interface to find equilibrium structure
- Keep bulk layers fixed to maintain bulk properties
- Reduce computational cost by limiting degrees of freedom

## Tips for Optimization

1. **Use `--use_builder_interface`**: Recommended for generating ordered, high-symmetry interfaces
2. **Atom Count Control**: Use `--max_atoms` to limit system size. The script automatically adjusts thickness if needed
3. **Selective Dynamics**: Use `-n 2` (2 layers per side) for most DFT calculations
4. **Layer Thickness**: Adjust `-t` (typically 2.0-3.0 Å) based on your material's layer spacing
5. **Strain Target**: 
   - `A`: Strain material A to match B (good when A is more flexible)
   - `B`: Strain material B to match A (good when B is more flexible)
   - `both`: Split strain evenly (good for balanced systems)

## Troubleshooting

- **No matches found**: Increase `--tol` or `--max_area`
- **Too many atoms**: Reduce `--max_atoms`, `--slab_thickness_a/b`, or `--max_area`
- **Poor matching**: Try different `--miller_a` and `--miller_b` combinations
- **Interface layers not relaxed**: Check `-n` and `-t` parameters in `fix_interface_layers.py`
- **Interface identification failed**: Try different `-m` method (`density_gap`, `median`, or `max_gap`)

## Recommended Workflow

**For most DFT calculations, use this two-step workflow:**

```bash
# Step 1: Generate interface
python build_heterojunctions/auto_interface_builder.py \
    --a <structure_A> \
    --b <structure_B> \
    --miller_a <h,k,l> \
    --miller_b <h,k,l> \
    --use_builder_interface \
    --max_atoms 300

# Step 2: Add selective dynamics
python build_heterojunctions/fix_interface_layers.py \
    auto_interface_combined.vasp \
    -o auto_interface_combined_relaxed.vasp \
    -n 2 \
    -t 2.0
```

This generates a structure ready for DFT relaxation calculations with interface layers relaxed and bulk layers fixed.

## Step 3: Stack Pre-relaxed Slabs (Alternative Workflow)

If you have pre-relaxed slab structures and want to stack them with a specific gap, you can use `stack_slabs.py`:

### Basic Usage

```bash
python build_heterojunctions/stack_slabs.py <slab1_file> <slab2_file> [gap_in_angstrom] [vacuum_in_angstrom] [output_file]
```

### Parameters

- `slab1_file`: Path to first slab file (bottom slab)
- `slab2_file`: Path to second slab file (top slab)
- `gap_in_angstrom`: Gap between two slabs in Angstroms (default: 3.0)
- `vacuum_in_angstrom`: Total vacuum layer thickness in Angstroms (default: 20.0), split equally between bottom and top
- `output_file`: Output filename (optional), if not specified will be auto-generated and saved to `relax_slab/` folder

### Examples

**Example 1**: Stack two pre-relaxed slabs with 3 Å gap and 20 Å vacuum

```bash
python build_heterojunctions/stack_slabs.py \
    relax_slab/fapbi3_slab_relax.vasp \
    relax_slab/tio2_fa_slab_relax.vasp \
    3.0 \
    20.0
```

This will create `relax_slab/fapbi3@tio2_fa_stacked.vasp` with structure:
- Bottom vacuum (10 Å, half of total 20 Å)
- Slab 1 (fapbi3)
- Gap (3 Å)
- Slab 2 (tio2_fa)
- Top vacuum (10 Å, half of total 20 Å)

**Example 2**: Stack slabs with custom gap and vacuum

```bash
python build_heterojunctions/stack_slabs.py \
    relax_slab/fapbi3_slab_relax.vasp \
    relax_slab/tio2_fa_slab_relax.vasp \
    5.0 \
    15.0
```

**Example 3**: Stack with different spacing and custom output name

```bash
python build_heterojunctions/stack_slabs.py \
    relax_slab/mapbi3_slab_relax.vasp \
    relax_slab/tio2_ma_slab_relax.vasp \
    3.2 \
    20.0 \
    custom_name.vasp
```

This will create `relax_slab/custom_name.vasp` (still saved in relax_slab folder)

### How It Works

1. Reads two VASP-format slab files
2. Calculates the thickness of each slab in the z-direction
3. Creates a structure with the following layout:
   - Bottom vacuum layer
   - Slab 1 (bottom slab)
   - Gap between slabs
   - Slab 2 (top slab)
   - Top vacuum layer
4. Adjusts z-coordinates of both slabs to fit in the new structure
5. Merges atomic coordinates and element information from both slabs
6. Updates the z-component of lattice vectors to accommodate both slabs, gap, and vacuum layers
7. Outputs the combined VASP file

### Output Information

The script displays:
- z-coordinate range and thickness for each slab
- Total number of atoms
- New lattice z-direction length
- Gap between the two slabs

### Notes

- Both slabs should have the same (or compatible) xy-plane lattice vectors
- Output file uses Direct coordinates (fractional coordinates)
- The script automatically merges elements of the same type
- Output file includes velocity information (initialized to zero)
- **Output location**: By default, stacked structures are saved to the `relax_slab/` folder
- **File naming**: Output files are named as `{slab1}@{slab2}_stacked.vasp` (e.g., `fapbi3@tio2_fa_stacked.vasp`)
- **Duplicate handling**: If a file with the same name already exists, a number suffix will be added (e.g., `fapbi3@tio2_fa_stacked_1.vasp`)
- **Vacuum layers**: The structure includes vacuum layers at both bottom and top. The specified vacuum thickness is the total (default: 20 Å), which is split equally between bottom and top (10 Å each) to ensure proper isolation for DFT calculations
- **Structure layout**: The final structure follows: [bottom vacuum] → [slab1] → [gap] → [slab2] → [top vacuum]

### Use Cases

This tool is useful when:
- You have pre-relaxed slab structures from separate calculations
- You want to create heterojunctions by stacking relaxed slabs
- You need precise control over the gap between slabs
- You want to combine slabs without going through the full interface builder workflow
