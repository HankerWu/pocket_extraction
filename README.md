# Pocket Extraction

**Pocket Extraction** is a Python package for extracting ligands and binding pockets from PDB files. It combines the power of **Biopython** and **RDKit** to provide flexible and efficient molecular structure processing.

## Features ✨
- **Extract Binding Pockets**: Identify pockets around ligands using coordinates, ligand files, or custom radii.
- **Extract Ligands**: Retrieve ligands by name, multiple ligands separately, or all HETATM residues (excluding solvents/ions).
- **Multi-Format Support**:  
  • **Input**: PDB, SDF, MOL2 (ligand files).  
  • **Output**: PDB (default), SDF, MOL2.
- **Advanced Filtering**: Select by model ID, chain ID, or ligand names.
- **Batch Processing**: Extract individual pockets for multiple ligands in one command.

---

## Installation

```bash
pip install pocket_extraction
```

---

## Quick Start 🚀

### 1. Extract Ligand and Its Pocket (CLI)
```bash
extract_ligand_and_pocket input.pdb \
  -l ligand.pdb \
  -p pocket.pdb \
  --ligand_names ATP \
  --radius 10.0
```

### 2. Extract All Ligands with Individual Pockets (CLI)
```bash
extract_ligand_and_pocket input.pdb \
  -l ligands/ \
  -p pockets/ \
  --multi_ligand
```

### 3. Python API Example
```python
from pocket_extraction import extract_ligand, extract_pocket

# Extract ligand "HEM" from Chain B
extract_ligand("input.pdb", "heme.pdb", 
              ligand_names=["HEM"], chain_id="B")

# Extract pocket around a manually defined center
extract_pocket("input.pdb", "pocket.pdb",
              ligand_center=[15.3, 24.7, 32.1], 
              radius=8.5)
```

---

## Usage Guide

### 🔍 Extracting Binding Pockets

#### Method 1: Ligand File (SDF/MOL2/PDB)
_Calculate pocket from ligand structure_  
**CLI**:
```bash
extract_pocket input.pdb --ligand_file ligand.sdf -o pocket.pdb --radius 12.5
```

**Python**:
```python
from pocket_extraction import extract_pocket, get_ligand_coords

ligand_coords = get_ligand_coords("ligand.mol2")
extract_pocket("input.pdb", "pocket.pdb", 
              ligand_coords=ligand_coords, 
              radius=12.5)
```

#### Method 2: Manual Coordinates
_Specify exact pocket center_  
**CLI**:
```bash
extract_pocket input.pdb --ligand_center 10.0 20.0 30.0 -o pocket.pdb
```

**Python**:
```python
extract_pocket("input.pdb", "pocket.pdb",
              ligand_center=[10.0, 20.0, 30.0],
              radius=10.0)  # Default radius
```

---

### ⚗️ Extracting Ligands

#### Case 1: Specific Ligand by Name
**CLI**:
```bash
extract_ligand input.pdb -o nad.pdb --ligand_names NAD
```

**Python**:
```python
extract_ligand("input.pdb", "nad.pdb",
              ligand_names=["NAD"],
              model_id=0,  # First model
              chain_id="A")
```

#### Case 2: Multiple Ligands Separately
**CLI**:
```bash
extract_ligand input.pdb -o ligands/ --ligand_names ATP NAD --multi_ligand
```
_Outputs_: `ligands/ligand_1.pdb`, `ligands/ligand_2.pdb`

**Python**:
```python
extract_ligand("input.pdb", "output_dir/",
              ligand_names=["ATP", "NAD"],
              multi_ligand=True)
```

#### Case 3: All Non-Solvent HETATM Residues
**CLI**:
```bash
extract_ligand input.pdb -o all_ligands.pdb
```

**Python**:
```python
extract_ligand("input.pdb", "all_ligands.pdb")
```

### Extracting Ligands and Pockets Together

Efficiently extract **ligands and their binding pockets** in one step using unified workflows. Choose from three modes depending on your use case:

---

#### Case 1: Merged Multi-Residue Ligands with Unified Pocket 
*Combine fragmented ligand residues into a single structure and extract their shared binding pocket.*  
  
**Use Cases**:  
• **Large/complex ligands** (e.g., ATP, NADH) split across multiple residues in PDB files  
• **Metal-cofactor systems** where ligands consist of multiple coordinated residues  
• **Cryo-EM/X-ray structures** with discontinuous ligand density assignments  

**CLI**:
```bash
extract_ligand_and_pocket input.pdb \
  -l ligand.pdb \
  -p pocket.pdb \
  --ligand_names HIS ARG \
  --model_id 0 \
  --chain_id E \
  --radius 12.0
```

**Python**:
```python
from pocket_extraction import extract_ligand_and_pocket

extract_ligand_and_pocket(
    pdb_path="input.pdb",
    ligand_output="ligand.pdb",
    pocket_output="pocket.pdb",
    ligand_names=["ATP", "NAD"],
    model_id=0,
    chain_id="E",
    radius=12.0
)
```

---

#### Case 2: Individual Pockets for Each Ligand  
*Extract ligands and pockets into separate files.*  
**Use Case**: Compare binding environments of distinct ligands.

**CLI**:
```bash
extract_ligand_and_pocket input.pdb \
  -l ligands/ \
  -p pockets/ \
  --ligand_names ATP NAD \
  --multi_ligand \
  --radius 10.0
```


**Python**:
```python
extract_ligand_and_pocket(
    pdb_path="input.pdb",
    ligand_output="ligands/",
    pocket_output="pockets/",
    ligand_names=["ATP", "NAD"],
    multi_ligand=True,
    radius=10.0
)
```

---

#### Case 3: Extract All Ligands & Pockets  
*Automatically process all non-solvent HETATM residues.*  
**Use Case**: High-throughput screening of unknown ligands.

**CLI**:
```bash
extract_ligand_and_pocket input.pdb \
  -l auto_ligands/ \
  -p auto_pockets/ \
  --multi_ligand \
  --radius 10.0
```

**Python**:
```python
extract_ligand_and_pocket(
    pdb_path="input.pdb",
    ligand_output="auto_ligands/",
    pocket_output="auto_pockets/",
    multi_ligand=True,
    radius=10.0
)
```

---

## License
MIT License. See [LICENSE](LICENSE) for details.

---

## Author
**Hanker Wu**  
📧 GitHub: [HankerWu](https://github.com/HankerWu/pocket_extraction)  
💬 *For bug reports or feature requests, please open a GitHub issue.*