import argparse
from Bio import PDB
from data_utils import save_structure
from selection import LigandSelect, PocketSelect
import numpy as np

def extract_ligand_and_pocket(pdb_file, ligand_file, pocket_file, ligand_names, multi_ligand, radius, ext):
    """Extract ligands and corresponding binding pockets from a PDB file."""

    # Load the PDB structure
    structure = PDB.PDBParser(QUIET=True).get_structure("STRUCTURE", pdb_file)

    # Extract ligands
    ligand_select = LigandSelect(ligand_names=ligand_names)
    ligand_structures = []
    for model in structure.get_list():
        for chain in model.get_list():
            for residue in chain.get_unpacked_list():
                if ligand_select.accept_residue(residue):
                    ligand_structures.append(residue)

    if not ligand_structures:
        raise ValueError("No ligands found matching the given criteria.")

    # Multi-ligand extraction: separate files for each ligand
    if multi_ligand:
        for i, ligand_residue in enumerate(ligand_structures):
            ligand_filename = ligand_file.replace(ext, f"_{i+1}{ext}")
            save_structure(ligand_filename, ligand_residue, ligand_select, ext)
            print(f"Ligand saved to {ligand_filename}")

            # Compute ligand center for each individual ligand
            ligand_coords = np.array([atom.coord for atom in ligand_residue.get_atoms()])

            # Extract and save corresponding pocket
            pocket_filename = pocket_file.replace(ext, f"_{i+1}{ext}")
            pocket_select = PocketSelect(radius=radius, ligand_coords=ligand_coords)
            save_structure(pocket_filename, structure, pocket_select, ext)
            print(f"Pocket saved to {pocket_filename}")

    # Single-ligand mode: merge all ligands and extract one pocket
    else:
        save_structure(ligand_file, structure, ligand_select, ext)
        print(f"Ligand saved to {ligand_file}")

        # Compute ligand center from all ligands combined
        ligand_coords = np.array([atom.coord for ligand in ligand_structures for atom in ligand.get_atoms()])

        # Extract and save pocket
        pocket_select = PocketSelect(radius=radius, ligand_coords=ligand_coords)
        save_structure(pocket_file, structure, pocket_select, ext)
        print(f"Pocket saved to {pocket_file}")

def main():
    parser = argparse.ArgumentParser(description="Extract ligands and their binding pockets from a PDB file.")
    parser.add_argument("pdb_file", type=str, help="Input PDB file")
    parser.add_argument("-l", "--ligand_file", type=str, default="ligand.pdb", help="Output ligand file")
    parser.add_argument("-p", "--pocket_file", type=str, default="pocket.pdb", help="Output pocket file")
    parser.add_argument("--ligand_names", type=str, nargs="+", default=None, help="Ligand names in PDB file")
    parser.add_argument("--multi_ligand", action="store_true", help="Extract multiple ligands separately")
    parser.add_argument("--radius", type=float, default=10.0, help="Pocket search radius")
    parser.add_argument("--ext", type=str, choices=[".pdb", ".cif"], default=".pdb", help="Output file format")

    args = parser.parse_args()
    extract_ligand_and_pocket(
        args.pdb_file, args.ligand_file, args.pocket_file,
        args.ligand_names, args.multi_ligand, args.radius, args.ext
    )

if __name__ == "__main__":
    main()
