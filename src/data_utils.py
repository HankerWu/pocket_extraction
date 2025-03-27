from Bio.PDB import Select
from Bio import PDB
import numpy as np
from typing import Optional, List
from os.path import dirname

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
            
def save_structure(filename, structure, select, ext):
    if ext == ".cif":
        io = PDB.MMCIFIO()
    else:
        io = PDB.PDBIO()
    io.set_structure(structure)
    try:
        io.save(str(filename), select=select)
    except TypeError as e:
        # Clean bad file.
        if filename.exists():
            filename.unlink()
        if str(e) == "%c requires int or char":
            raise RuntimeError("Chain id is not a single char") from e
        else:
            raise
        
        
