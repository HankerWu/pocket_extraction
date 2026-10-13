# Pocket Extraction 

**Pocket Extraction** is a Python package built on **Biopython** for extracting ligands and binding pockets from structural biology files (PDB/mmCIF). It supports high-throughput screening as well as detailed structural analyses.

---

## Key Features ✨

- **Binding Pocket Extraction**  
  Extract pockets around ligands using either:
  - An existing ligand file, or  
  - Manually specified coordinates  
  *(Adjust search radius with `--radius`)*

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

### Global Arguments (All Commands)
| Argument | Description | Required | Default | Type |
|----------|-------------|----------|---------|------|
| `pdb_file` | Input PDB/mmCIF file or PDB ID | Yes | - | Path/Str |
| `-o/--output` | Output path (file/directory) | No | `output.pdb` | Path |
| `--ext` | Override output format | No | Auto-detect | {pdb, cif} |
| `-q/--quiet` | Suppress informational output | No | `False` | Flag |
| `--debug` | Enable debug logging | No | `False` | Flag |
| `--logfile` | Path to save log file | No | None | Path |
| `--ignore_duplicates` | Whether to ignore duplicate ligands | No | `True` | Str2Bool |
---

### `extract_ligand` (Ligand Extraction)
| Argument | Description | Required | Default | Type |
|----------|-------------|----------|---------|------|
| `--ligand_names` | Ligand residue names to extract (e.g., `RLZ HEM`) | No | All non-solvent | 1+ Strings |
| `--exclude` | Residues to exclude (e.g., `HOH`) | No | None | 1+ Strings |
| `-m/--model_ids` | Model IDs (0-based) to process | No | All models | 1+ Integers |
| `-c/--chain_ids` | Chain IDs to process | No | All chains | 1+ Strings |
| `--multi` | Save separate files per ligand | No | `False` | Flag |

---

### `extract_pocket` (Pocket Extraction)
| Argument | Description | Required | Default | Type |
|----------|-------------|----------|---------|------|
| `--ligand_file` | Reference ligand structure file | Mutually Exclusive | - | Path |
| `--ligand_center` | Manual center coordinates (X Y Z) | Mutually Exclusive | - | 3 Floats |
| `-r/--radius` | Pocket radius in Angstroms | No | 10.0 | Float |

---

### `extract_ligand_and_pocket` (Combined Extraction)
| Argument | Description | Required | Default | Type |
|----------|-------------|----------|---------|------|
| `-ol/--output_ligand` | Output path for ligands | No | `ligand.pdb` | Path |
| `-op/--output_pocket` | Output path for pockets | No | `pocket.pdb` | Path |
| `-r/--radius` | Pocket radius around ligands | No | 10.0 | Float |
| *Inherits all arguments from `extract_ligand`* | | | | |

## 💡 Usage Examples

`pocket_extraction` supports both **command-line** and **Python API** usage. Below are examples demonstrating typical workflows using the structure `7AHN`.

---

### 🧪 1. **Extracting Binding Pockets**

#### 🔧 Command-Line

```bash
# Extract a binding pocket around an existing ligand file
extract_pocket 7AHN.pdb -o 7AHN_pocket.cif --ligand_file 7AHN_ligand.pdb --radius 12.5 --quiet

# Extract a pocket using manually provided coordinates
extract_pocket 7AHN -o 7AHN_pocket.pdb --ligand_center 117.21642 132.5165 129.84128 --radius 10.0
```

#### 🐍 Python

```python
from pocket_extraction import extract_pocket, get_ligand_coords

# Option 1: Use coordinates extracted from a ligand file
ligand_coords = get_ligand_coords("7AHN_ligand.pdb")
extract_pocket("7AHN.pdb", "7AHN_pocket.pdb", ligand_coords=ligand_coords, radius=12.5, quiet=True)

# Option 2: Use manually specified coordinates
extract_pocket("7AHN", "7AHN_pocket.cif", ligand_center=[117.21642, 132.5165, 129.84128], radius=10.0, quiet=True)
```

---

### 💊 2. **Extracting Ligands**

#### 🔧 Command-Line

```bash
# Extract a specific ligand by name and chain
extract_ligand 7AHN -o 7AHN_RLZ_B_ligand.pdb --ligand_names RLZ --chain_ids B --quiet

# Extract multiple ligands, excluding HIC, and save each one separately
extract_ligand 7AHN -o 7AHN_ligands/ --multi --exclude HIC
```

#### 🐍 Python

```python
from pocket_extraction import extract_ligand

# Extract a specific ligand with optional model and chain filtering
extract_ligand("7AHN.pdb", "7AHN_RLZ_B_ligand.pdb", ligand_names=["RLZ"], model_ids=[0], chain_ids=["B"], quiet=True)

# Extract multiple ligands, each saved individually
extract_ligand("7AHN", "7AHN_ligands/", ligand_names=["RLZ", "HIC"], multi=True)
```

---

### 🔗 3. **Extracting Both Ligands and Binding Pockets**

Each pocket uses ligand coordinates from the same structure model. Model filters
and duplicate-ligand selection also restrict which models contribute pocket
residues. Combined output keeps the selected models together; `--multi` writes
each ligand's own-model pocket. Pocket selection may include nearby residues
from other chains in that model, preserving protein interface pockets.

#### 🔧 Command-Line

```bash
# Example 1: Extract ligand (RLZ) and its pocket from model 0, chain B
extract_ligand_and_pocket 7AHN \
  -ol 7AHN_RLZ.pdb \
  -op 7AHN_pocket.pdb \
  --ligand_names RLZ \
  -m 0 \
  -c B \
  -r 12.0 \
  -q

# Example 2: Extract RLZ and HIC ligands and their pockets, save separately
extract_ligand_and_pocket 7AHN.pdb \
  -ol 7AHN_ligands/ \
  -op 7AHN_pockets/ \
  --ligand_names RLZ HIC \
  --multi \
  -r 10.0 

# Example 3: Automatically extract all non-solvent ligands and pockets
extract_ligand_and_pocket 7AHN.pdb \
  -ol auto_ligands/ \
  -op auto_pockets/ \
  --multi \
  -r 10.0 \
  -q
```

#### 🐍 Python

```python
from pocket_extraction import extract_ligand_and_pocket

extract_ligand_and_pocket(
    pdb_file="7AHN",             # or your structural file path
    output_ligand="ligand.pdb",  # or a directory (e.g., "ligands/")
    output_pocket="pocket.pdb",  # or a directory (e.g., "pockets/")
    ligand_names=["RLZ"],        # omit to auto detect
    model_ids=[0],               # optional filtering
    chain_ids=["B"],             # optional filtering
    multi=True,                  # set to True for separate files
    radius=12.0,                 # adjust the search radius
    quiet=True,                  # suppress output messages
    exclude=["HIC"]              # exclude specific ligands or chains
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


## Ligand output semantics

`extract_ligand` and `extract_ligand_and_pocket` share filtered ligand collection.
Without `--multi`, the single ligand output contains **all selected residues**,
retaining their original model, chain, residue and atom identities; the paired
pocket is selected around that same collection. `--ignore_duplicates true` keeps
the first residue of each ligand name after model/chain/name filters. With
`--multi`, existing per-ligand output naming and pocket selection are retained.
The returned count reports the number of selected ligand residues, not files.
Historical extracted files are not regenerated by updating the implementation.

The exported default exclusion collection `REMOVE_LIGANDS` is immutable.
Additional `exclude` names are local to each Python extraction call and override
explicit `ligand_names`; the two ligand CLIs pass `--exclude` through this same
selector without changing the package's default blocklist. Default automatic
extraction still applies that blocklist; explicit names still permit a requested
residue outside the automatic selection.

Both pocket APIs require a finite positive radius. Reference coordinates must
be a nonempty finite Nx3 array, or a finite three-coordinate center. These
geometry checks run before output publication. Ligand and pocket destinations
must resolve to distinct file paths; colliding destinations are rejected before
the first file is written, including multi mode. A protein residue is included
when any of its atoms is within the requested radius; approximate average amino
acid sizes cannot exclude a qualifying side chain.

`save_structure` writes a same-directory temporary PDB/mmCIF and replaces the
requested file only after successful serialization, retaining an existing file's
permission mode and the normal creation umask for new files. A failed writer
removes its temporary file and preserves the original output and exception.
Each file is published separately: combined ligand/pocket and multi-file outputs
are not one atomic transaction. A later filesystem failure can leave an already
completed file alongside the previous version of another file. No frozen
structure or downstream scientific result is regenerated by these fixes.

The three CLI parsers build only their supported arguments. Standalone pocket
extraction requires exactly one reference file or center; combined extraction
uses the selected ligand coordinates and has no reference or generic output
argument. The existing `requests` download client is declared as a dependency;
offline local-file extraction needs no network request.
