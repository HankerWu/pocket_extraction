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

## Command-line Arguments

Different command-line tools support different arguments. Below is a breakdown of supported options per tool.

### **`extract_pocket`**
- `--ligand_file <file>`: Specify a ligand file to extract a binding pocket around.
- `--ligand_center <x y z>`: Provide manual coordinates as the center of the pocket.
- `--radius <value>`: Define the search radius for the pocket extraction.
- `--model_id <id>`: Filter by model ID.
- `--chain_id <id>`: Filter by chain ID.
- `--exclude <residue(s)>`: Exclude specific residues (e.g., `HOH` for water molecules).
- `--quiet`: Suppress output messages for silent execution.

### **`extract_ligand`**
- `--ligands <name(s)>`: Specify one or more ligand names for extraction.
- `--multi`: Save each ligand separately when extracting multiple.
- `--model_id <id>`: Filter by model ID.
- `--chain_id <id>`: Filter by chain ID.
- `--exclude <residue(s)>`: Exclude specific residues (e.g., `HOH` for water molecules).
- `--quiet`: Suppress output messages for silent execution.

### **`extract_ligand_and_pocket`**
- `--ligands <name(s)>`: Specify one or more ligand names for extraction.
- `--multi`: Save each ligand and pocket separately when extracting multiple.
- `--model_id <id>`: Filter by model ID.
- `--chain_id <id>`: Filter by chain ID.
- `--radius <value>`: Define the search radius for pocket extraction.
- `--exclude <residue(s)>`: Exclude specific residues (e.g., `HOH` for water molecules).
- `--quiet`: Suppress output messages for silent execution.

---

## Usage Examples

The package provides both CLI and Python interfaces. Below are examples of common use cases.

### 1. Extracting Binding Pockets

#### CLI:
```bash
# Using an existing ligand file:
extract_pocket input.pdb -o pocket.cif --ligand_file ligand.pdb --radius 12.5 --quiet

# Using manual coordinates (specify ligand center):
extract_pocket input.cif -o pocket.pdb --ligand_center 10.0 20.0 30.0 --radius 10.0 --exclude HOH
```

#### Python:
```python
from pocket_extraction import extract_pocket, get_ligand_coords

# Option A: Using an existing ligand file
ligand_coords = get_ligand_coords("ligand.pdb")
extract_pocket("input.pdb", "pocket.pdb", ligand_coords=ligand_coords, radius=12.5, quiet=True)

# Option B: Using manually provided coordinates
extract_pocket("input.cif", "pocket.cif", ligand_center=[10.0, 20.0, 30.0], radius=10.0, exclude=["HOH"])
```

---

### 2. Extracting Ligands

#### CLI:
```bash
# Extract a specific ligand or multiple ligands:
extract_ligand input.pdb -o output_path --ligands NAD --quiet        # Single ligand
extract_ligand input.cif -o output_dir --ligands ATP NAD --multi --exclude HOH  # Multiple ligands, each saved separately
```

#### Python:
```python
from pocket_extraction import extract_ligand

# Example for a specific ligand with optional filtering:
extract_ligand("input.pdb", "nad.pdb", ligand_names=["NAD"], model_id=0, chain_id="A", quiet=True)

# Example for multiple ligands:
extract_ligand("input.cif", "output_dir/", ligand_names=["ATP", "NAD"], multi_mode=True, exclude=["HOH"])
```

---

### 3. Combined Extraction of Ligands and Pockets

Use the combined function for simultaneous extraction:

#### CLI:
```bash
# Example 1: Merged multi-residue ligand with a unified pocket
extract_ligand_and_pocket input.pdb -l ligand.pdb -p pocket.pdb --ligands HIS ARG --model_id 0 --chain_id E --radius 12.0 --quiet

# Example 2: Separate files for each ligand and pocket:
extract_ligand_and_pocket input.pdb -l ligands/ -p pockets/ --ligands ATP NAD --multi --radius 10.0 --exclude HOH

# Example 3: Automatic extraction of all non-solvent ligands and pockets:
extract_ligand_and_pocket input.pdb -l auto_ligands/ -p auto_pockets/ --multi --radius 10.0 --quiet
```

#### Python:
```python
from pocket_extraction import extract_ligand_and_pocket

extract_ligand_and_pocket(
    pdb_file="input.pdb",      # or pdb_path for automatic extraction
    ligand_path="ligand.pdb",  # or a directory (e.g., "ligands/")
    pocket_path="pocket.pdb",  # or a directory (e.g., "pockets/")
    ligand_names=["ATP", "NAD"],   # omit or adjust as needed
    model_id=0,                # optional filtering
    chain_id="E",              # optional filtering
    multi_mode=True,           # set to True for separate files
    radius=12.0,               # adjust the search radius
    quiet=True,                # suppress output messages
    exclude=["HOH"]          # exclude specific ligands or chains
)
```

---

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

---

## Author

**Hanker Wu**  
📧 GitHub: [HankerWu](https://github.com/HankerWu/pocket_extraction)  
💬 *For bug reports or feature requests, please open a GitHub issue.*

