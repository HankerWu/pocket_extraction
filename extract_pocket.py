import argparse
from Bio import PDB
from data_utils import save_structure
from selection import PocketSelect
from rdkit import Chem
import numpy as np

def check_ligand_file(ligand_file):
    """Check if the ligand file is valid and contains a ligand."""
    ext = ligand_file.split(".")[-1]
    if ext not in ["sdf", "mol2", "pdb"]:
        raise ValueError(f"Unsupported ligand file format: {ext}")
    return ext

def get_ligand_coords(ligand_file):
    """Extract ligand coordinates from the ligand file."""
    ext = check_ligand_file(ligand_file)
    if ext == "sdf":
        return _get_coords_from_sdf(ligand_file)
    elif ext == "mol2":
        return _get_coords_from_mol2(ligand_file)
    elif ext == "pdb":
        return _get_coords_from_pdb(ligand_file)
    else:
        raise ValueError(f"Unsupported ligand file format: {ext}")

def _get_coords_from_sdf(ligand_file):
    """Extract coordinates from SDF file."""
    suppl = Chem.SDMolSupplier(ligand_file)
    coords = []
    for mol in suppl:
        if mol is not None:
            conf = mol.GetConformer()
            for atom in mol.GetAtoms():
                coords.append(conf.GetAtomPosition(atom.GetIdx()))
    return np.array(coords)

def _get_coords_from_mol2(ligand_file):
    """Extract coordinates from MOL2 file."""
    suppl = Chem.MolFromMol2File(ligand_file)
    coords = []
    if suppl is not None:
        conf = suppl.GetConformer()
        for atom in suppl.GetAtoms():
            coords.append(conf.GetAtomPosition(atom.GetIdx()))
    return np.array(coords)

def _get_coords_from_pdb(ligand_file):
    """Extract coordinates from PDB file."""
    ligand_structure = PDB.PDBParser(QUIET=True).get_structure("LIGAND", ligand_file)
    return np.array([atom.coord for atom in ligand_structure.get_atoms()])

def main():
    parser = argparse.ArgumentParser(description="Extract pocket from PDB file")
    parser.add_argument("pdb_file", type=str, help="Input PDB file")
    parser.add_argument("-o", "--pocket_file", type=str, default="pocket.pdb", help="Output pocket file")
    parser.add_argument("--ligand_file", type=str, default=None, help="Input ligand file (optional)")
    parser.add_argument("--ligand_center", type=float, nargs=3, default=None, help="Ligand center coordinates (optional)")
    parser.add_argument("--radius", type=float, default=10.0, help="Pocket search radius (default: 10.0)")
    parser.add_argument("--ext", type=str, default=".pdb", choices=[".pdb", ".cif"], help="File extension (default: .pdb)")
    args = parser.parse_args()
    
    # Load the structure
    structure = PDB.PDBParser(QUIET=True).get_structure("POCKET", args.pdb_file)
    # Create a selection object for the pocket
    if args.ligand_file:
        ligand_coords = get_ligand_coords(args.ligand_file)
        ligand_center = None
    else:
        ligand_coords = None
        assert args.ligand_center is not None, "Must provide either ligand_file or ligand_center"
        ligand_center = args.ligand_center
        
    select = PocketSelect(
        ligand_coords=ligand_coords,
        ligand_center=ligand_center,
        radius=args.radius
    )
    # Save the pocket structure
    save_structure(args.pocket_file, structure, select, args.ext)
    print(f"Pocket saved to {args.pocket_file}")
    
if __name__ == "__main__":
    main()