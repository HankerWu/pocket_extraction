from Bio.PDB import Select
from Bio import PDB
import numpy as np
from typing import Optional, List
import os
from os.path import dirname
from pathlib import Path

# Mapping of three-letter residue codes to one-letter codes
residue_name_to_one_letter_code = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M", "UNK": "X"  # Selenomethionine and unknown residue
}

# Bounding sphere radii for amino acids (in Angstroms)
amino_acid_bounding_radii = {
    "A": 1.5, "C": 1.7, "D": 2.0, "E": 2.2, "F": 2.8, "G": 1.0,
    "H": 2.5, "I": 2.2, "K": 2.8, "L": 2.2, "M": 2.3, "N": 2.0,
    "P": 1.9, "Q": 2.2, "R": 3.0, "S": 1.6, "T": 1.9, "V": 2.0,
    "W": 2.8, "Y": 2.5  # Based on average side chain dimensions
}

def get_remove_ligands():
    # Simplified quality assessment for small-molecule ligands in the Protein Data Bank
    # https://doi.org/10.1016/j.str.2021.10.003
    # https://github.com/rcsb/PDB_ligand_quality_composite_score
    remove_ligands = []
    with open(dirname(__file__) + "/non-LOI-blocklist.tsv", "r") as f:
        for line in f:
            remove_ligands.append(line.strip().split("\t")[0])
    remove_ligands += ["DMS", "ZN", "SO4", "GOL", "BTB"]
    remove_ligands = set(remove_ligands)
    return remove_ligands

def process_output_path(output_path: str, base_name: str, ext: Optional[str] = None, index: Optional[int] = None) -> str:
    """Process output path with proper extension handling.
    
    Args:
        output_path: User-specified output path (could be file or directory)
        base_name: Default base name if path is directory
        ext: User-specified extension override
        index: Counter for multi-file mode
    
    Returns:
        Full output path with proper extension
    """
    output_dir, filename = os.path.split(output_path)
    
    # Handle directory path case
    if not filename:
        filename = f"{base_name}_{index}" if index else base_name
    
    # Split name and original extension
    name_part, orig_ext = os.path.splitext(filename)
    orig_ext = orig_ext.lstrip('.').lower()  # Normalize extension
    
    # Determine final extension priority: user specified > filename extension > default pdb
    final_ext = (ext or orig_ext or 'pdb').lower()
    
    # Validate supported formats
    if final_ext not in ('pdb', 'cif'):
        raise ValueError(f"Unsupported output format: {final_ext}. Use 'pdb' or 'cif'.")
    
    # Construct final filename
    final_name = f"{name_part}.{final_ext}"
    
    # Create directory if needed
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    
    return os.path.join(output_dir, final_name)

def load_structure(pdb_file):
    # Load the structure
    input_ext = pdb_file.split(".")[-1]
    if input_ext == "pdb":
        structure = PDB.PDBParser(QUIET=True).get_structure("pdb_file", pdb_file)
    elif input_ext == "cif":
        structure = PDB.MMCIFParser(QUIET=True).get_structure("pdb_file", pdb_file)
    else:
        raise ValueError(f"Unsupported file extension: {input_ext}")
    return structure
            
def save_structure(filename, structure, select=PDB.Select()):
    """Save structural data to file with robust error handling.
    
    Args:
        filename: Output path (str or Path)
        structure: Biopython Structure object
        select: Residue selection criteria (default: all)
        
    Raises:
        ValueError: Unsupported file format
        RuntimeError: Chain ID issues or cleanup failures
        IOError: File system errors during save
    """
    path = Path(filename).resolve()
    suffix = path.suffix.lower()
    
    # Validate file format
    if suffix not in (".pdb", ".cif"):
        raise ValueError(
            f"Unsupported format '{suffix}'. Use .pdb or .cif"
        )

    # Initialize appropriate IO object
    io = PDB.MMCIFIO() if suffix == ".cif" else PDB.PDBIO()
    io.set_structure(structure)

    try:
        # Ensure parent directory exists
        path.parent.mkdir(parents=True, exist_ok=True)
        io.save(str(path), select=select)
        
    except Exception as save_error:
        # Clean up partial files on any error
        if path.exists():
            try:
                path.unlink()
            except Exception as cleanup_error:
                raise RuntimeError(
                    f"Failed to clean up corrupted file: {cleanup_error}"
                ) from save_error

        # Handle specific PDB format constraints
        if (
            suffix == ".pdb" 
            and isinstance(save_error, TypeError)
            and "requires int or char" in str(save_error)
        ):
            raise RuntimeError(
                "Invalid PDB chain ID detected. Must be single character."
            ) from save_error

        # Re-raise original error with additional context
        raise type(save_error)(
            f"Failed to save structure to {path}: {save_error}"
        ) from save_error
        
