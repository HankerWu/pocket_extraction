# Protein-Ligand Pocket Extraction

## Overview
This project provides Python scripts for extracting ligand binding pockets and ligands from Protein Data Bank (PDB) files. The scripts utilize the Biopython and RDKit libraries to process molecular structures efficiently.

## Features
- **extract_pocket.py**: Extracts the binding pocket around a ligand in a PDB file using different methods.
- **extract_ligand.py**: Extracts ligands from a PDB file using multiple selection criteria.
- **selection.py**: Implements selection logic for identifying pockets and ligands.
- **data_utils.py**: Provides utility functions for handling molecular structures.

## Dependencies
The following Python packages are required:
- `biopython`
- `rdkit`
- `numpy`

You can install them using:
```sh
pip install biopython rdkit numpy
```

## Usage

### Extracting a Pocket
```sh
python extract_pocket.py input.pdb --ligand_file ligand.sdf --radius 10.0 -o pocket.pdb
```
**Extraction Methods:**
1. **Ligand Center Mode:**
   - Specify the ligand center manually using `--ligand_center x y z`.
   - Example:
     ```sh
     python extract_pocket.py input.pdb --ligand_center 10.0 20.0 30.0 --radius 10.0 -o pocket.pdb
     ```
   - This mode does not require a ligand file and relies on a manually provided center.
2. **Full Ligand Atom Coordinates Mode:**
   - Extracts pocket residues based on the full ligand structure using `--ligand_file`.
   - Example:
     ```sh
     python extract_pocket.py input.pdb --ligand_file ligand.sdf --radius 10.0 -o pocket.pdb
     ```
   - Supported ligand file formats:
     - **SDF (`.sdf`)**: Structure-data file format commonly used for small molecules.
     - **MOL2 (`.mol2`)**: Tripos MOL2 format containing molecular structure and properties.
     - **PDB (`.pdb`)**: PDB format containing atomic coordinates.
   - If the ligand file is provided, the center is automatically computed from ligand atom coordinates.
3. **Pocket Search Radius:**
   - The `--radius` argument determines how far from the ligand the pocket residues will be extracted.
   - Default is `10.0 Å`, but this can be adjusted based on binding site size.

### Extracting Ligands
```sh
python extract_ligand.py input.pdb -o ligand.pdb --ligand_names ATP
```
**Extraction Methods:**
1. **Extract Specific Ligands by Name:**
   - Use `--ligand_names` to specify one or more residue names corresponding to ligands in the PDB file.
   - Example (extract ATP and NAD as a single ligand file):
     ```sh
     python extract_ligand.py input.pdb -o ligand.pdb --ligand_names ATP NAD
     ```
   - If `--multi_ligand` is **not** specified, all residues in `--ligand_names` will be merged into a single output file.
2. **Extract Multiple Ligands Separately:**
   - Use `--multi_ligand` to extract each specified ligand as a separate file.
   - Example (extract ATP and NAD into separate files):
     ```sh
     python extract_ligand.py input.pdb -o ligand.pdb --ligand_names ATP NAD --multi_ligand
     ```
3. **Extract All HETATM Ligands:**
   - If `--ligand_names` is **not** specified, the script extracts all HETATM residues (non-standard residues).
   - If `--multi_ligand` is **not** specified, all extracted ligands will be merged into one file.
   - If `--multi_ligand` **is** specified, ligands will be split into separate files based on their residue names.
   - Example (extract all ligands into one file):
     ```sh
     python extract_ligand.py input.pdb -o ligand.pdb
     ```
   - Example (extract all ligands into separate files based on residue name):
     ```sh
     python extract_ligand.py input.pdb -o ligand.pdb --multi_ligand
     ```
4. **Filter by Model and Chain:**
   - Extract ligands from a specific model and chain using `--model_id` and `--chain_id`.
   - Example:
     ```sh
     python extract_ligand.py input.pdb -o ligand.pdb --model_id 0 --chain_id A
     ```

### **Extract Ligands and Corresponding Pockets**
The `extract_ligand_and_pocket.py` script allows you to extract both ligands and their corresponding binding pockets in one step.

```sh
python extract_ligand_and_pocket.py input.pdb -l ligand.pdb -p pocket.pdb --ligand_names ATP --radius 10.0
```

#### **Extraction Modes:**
1. **Single Pocket for All Ligands (Default Mode)**
   - If multiple ligands are found, they will be merged into a single `ligand.pdb` file.
   - The pocket will be extracted based on the center of all ligand atoms combined.
   - Example:
     ```sh
     python extract_ligand_and_pocket.py input.pdb -l ligand.pdb -p pocket.pdb --ligand_names ATP NAD --radius 10.0
     ```

2. **Multiple Ligands and Individual Pockets (`--multi_ligand`)**
   - If `--multi_ligand` is specified, each ligand will be extracted into a separate file.
   - Each pocket will be extracted based on the center of its corresponding ligand.
   - Example:
     ```sh
     python extract_ligand_and_pocket.py input.pdb -l ligand.pdb -p pocket.pdb --ligand_names ATP NAD --multi_ligand --radius 10.0
     ```

3. **Extract All Ligands and Their Pockets**
   - If `--ligand_names` is not provided, the script will extract all HETATM ligands.
   - Each ligand will have a corresponding pocket extracted.
   - Example:
     ```sh
     python extract_ligand_and_pocket.py input.pdb -l ligand.pdb -p pocket.pdb --multi_ligand --radius 10.0
     ```
