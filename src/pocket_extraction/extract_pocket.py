import argparse
import os
import logging
from Bio import PDB
from rdkit import Chem
import numpy as np
from .data_utils import load_structure, save_structure, process_output_path
from .selection import PocketSelect

logger = logging.getLogger(__name__)

def get_ligand_coords(ligand_file):
    """Extract ligand coordinates from various chemical file formats."""
    ext = os.path.splitext(ligand_file)[1].lower().lstrip('.')
    
    if ext == "sdf":
        return _get_coords_from_sdf(ligand_file)
    if ext == "mol2":
        return _get_coords_from_mol2(ligand_file)
    if ext in ["pdb", "ent"]:
        return _get_coords_from_pdb(ligand_file)
    if ext in ["cif", "mmcif"]:
        return _get_coords_from_cif(ligand_file)
    
    raise ValueError(f"Unsupported ligand file format: {ext}")

def _get_coords_from_sdf(sdf_file):
    """Extract coordinates from SDF file using RDKit."""
    suppl = Chem.SDMolSupplier(sdf_file)
    for mol in suppl:
        if mol:
            return np.array([mol.GetConformer().GetAtomPosition(i) for i in range(mol.GetNumAtoms())])
    raise ValueError(f"No valid molecules in SDF file: {sdf_file}")

def _get_coords_from_mol2(mol2_file):
    """Extract coordinates from MOL2 file using RDKit."""
    mol = Chem.MolFromMol2File(mol2_file)
    if not mol:
        raise ValueError(f"Failed to read MOL2 file: {mol2_file}")
    return np.array([mol.GetConformer().GetAtomPosition(i) for i in range(mol.GetNumAtoms())])

def _get_coords_from_pdb(pdb_file):
    """Extract coordinates from PDB file using Biopython."""
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("ligand", pdb_file)
    return np.array([atom.coord for atom in structure.get_atoms()])

def _get_coords_from_cif(cif_file):
    """Extract coordinates from CIF file using Biopython."""
    parser = PDB.MMCIFParser(QUIET=True)
    structure = parser.get_structure("ligand", cif_file)
    return np.array([atom.coord for atom in structure.get_atoms()])


def extract_pocket(pdb_file, output_path, ligand_coords=None, ligand_center=None, radius=10.0, ext=None):
    """Main pocket extraction logic with improved file handling."""
    # Validate input coordinates
    if ligand_coords is None and ligand_center is None:
        raise ValueError("Must provide either ligand coordinates or center")
    
    # Process output path
    output_path = process_output_path(output_path, "pocket", ext)
    
    # Load structure and create selector
    structure = load_structure(pdb_file)
    selector = PocketSelect(
        ligand_coords=ligand_coords,
        ligand_center=ligand_center,
        radius=radius
    )
    
    # Save the pocket structure
    save_structure(output_path, structure, selector)
    return output_path

def main():
    parser = argparse.ArgumentParser(description="Extract protein pocket based on ligand proximity")
    parser.add_argument("pdb_file", help="Input protein structure (PDB/mmCIF)")
    parser.add_argument("-o", "--output", default="pocket.pdb",
                      help="Output path for pocket structure (directory or filename)")
    parser.add_argument("--ligand_file", help="Ligand structure file (SDF/MOL2/PDB/CIF)")
    parser.add_argument("--ligand_center", nargs=3, type=float,
                      help="Manual ligand center coordinates (X Y Z)")
    parser.add_argument("-r", "--radius", type=float, default=10.0,
                      help="Pocket radius in Angstroms (default: 10.0)")
    parser.add_argument("--ext", choices=["pdb", "cif"], 
                      help="Force output format (default: auto from filename)")
    
    args = parser.parse_args()
    
    # Validate input
    if not (args.ligand_file or args.ligand_center):
        parser.error("Must provide either --ligand_file or --ligand_center")
    
    # Get ligand coordinates
    ligand_coords = None
    if args.ligand_file:
        try:
            ligand_coords = get_ligand_coords(args.ligand_file)
        except Exception as e:
            logger.error(f"Failed to process ligand file: {str(e)}")
            return
    
    # Extract and save pocket
    try:
        output_path = extract_pocket(
            pdb_file=args.pdb_file,
            output_path=args.output,
            ligand_coords=ligand_coords,
            ligand_center=args.ligand_center,
            radius=args.radius,
            ext=args.ext
        )
        print(f"Successfully saved pocket to: {output_path if output_path else '.'}")
    except Exception as e:
        logger.error(f"Pocket extraction failed: {str(e)}")

if __name__ == "__main__":
    main()
    