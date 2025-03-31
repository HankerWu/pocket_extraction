# Pocket Extraction

**Pocket Extraction** is a Python package built on **Biopython** for extracting ligands and binding pockets from structural biology files (PDB/mmCIF). It supports high-throughput screening as well as detailed structural analyses.

---

## Key Features ✨

- **Binding Pocket Extraction**  
  Extract pockets around ligands using either:
  - An existing ligand file, or  
  - Manually specified coordinates  
  *(Adjust the search radius as needed.)*

- **Ligand Extraction**  
  Retrieve ligands by specifying names (single/multiple) or by automatically processing all non-solvent HETATM residues.

- **Flexible I/O Support**  
  - **Input:** PDB, mmCIF  
  - **Output:** PDB (default), mmCIF

- **Advanced Filtering & Batch Processing**  
  Filter by model ID, chain ID, or ligand names; process multiple ligands/pockets in one command.

---

## Installation

Install via pip:

```bash
pip install pocket_extraction
```

---

## Usage Examples

The package provides both CLI and Python interfaces. Below are the consolidated examples that cover the most common use cases.

### 1. Extracting Binding Pockets

**CLI:**
```bash
# Using an existing ligand file:
extract_pocket input.pdb -o pocket.cif --ligand_file ligand.pdb --radius 12.5

# Using manual coordinates (specify ligand center):
extract_pocket input.cif -o pocket.pdb --ligand_center 10.0 20.0 30.0 --radius 10.0
```

**Python:**
```python
from pocket_extraction import extract_pocket, get_ligand_coords

# Option A: Using an existing ligand file
ligand_coords = get_ligand_coords("ligand.pdb")
extract_pocket("input.pdb", "pocket.pdb", ligand_coords=ligand_coords, radius=12.5)

# Option B: Using manually provided coordinates
extract_pocket("input.cif", "pocket.cif", ligand_center=[10.0, 20.0, 30.0], radius=10.0)
```

*Note: Choose the option that matches your data source. The only difference is in providing either a ligand file (or coordinates obtained from it) or manual coordinates.*

---

### 2. Extracting Ligands

**CLI:**
```bash
# For a specific ligand or multiple ligands:
extract_ligand input.pdb -o output_path --ligands NAD         # Single ligand
extract_ligand input.cif -o output_dir --ligands ATP NAD --multi  # Multiple ligands, each saved separately
```

**Python:**
```python
from pocket_extraction import extract_ligand

# Example for a specific ligand with optional filtering parameters:
extract_ligand("input.pdb", "nad.pdb", ligand_names=["NAD"], model_id=0, chain_id="A")

# Example for multiple ligands:
extract_ligand("input.cif", "output_dir/", ligand_names=["ATP", "NAD"], multi_mode=True)
```

---

### 3. Combined Extraction of Ligands and Pockets

Use the combined function for simultaneous extraction. Adjust parameters based on your workflow:

**CLI:**
```bash
# Example 1: Merged multi-residue ligand with a unified pocket
extract_ligand_and_pocket input.pdb -l ligand.pdb -p pocket.pdb --ligand_names HIS ARG --model_id 0 --chain_id E --radius 12.0

# Example 2: Separate files for each ligand and pocket:
extract_ligand_and_pocket input.pdb -l ligands/ -p pockets/ --ligands ATP NAD --multi --radius 10.0

# Example 3: Automatic extraction of all non-solvent ligands and pockets:
extract_ligand_and_pocket input.pdb -l auto_ligands/ -p auto_pockets/ --multi --radius 10.0
```

**Python:**
```python
from pocket_extraction import extract_ligand_and_pocket

# Option: Customize parameters for your workflow.
extract_ligand_and_pocket(
    pdb_file="input.pdb",      # or pdb_path for automatic extraction
    ligand_path="ligand.pdb",  # or a directory (e.g., "ligands/")
    pocket_path="pocket.pdb",  # or a directory (e.g., "pockets/")
    ligand_names=["ATP", "NAD"],   # omit or adjust as needed
    model_id=0,                # optional filtering
    chain_id="E",              # optional filtering
    multi_mode=True,           # set to True for separate files
    radius=12.0                # adjust the search radius
)
```

*Note: The Python interface uses similar parameters to the CLI. The function adapts based on whether single files or directories are provided, and whether multi_mode is enabled.*

---

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

---

## Author

**Hanker Wu**  
📧 GitHub: [HankerWu](https://github.com/HankerWu/pocket_extraction)  
💬 *For bug reports or feature requests, please open a GitHub issue.*