import argparse
from Bio import PDB
from .data_utils import save_structure
from .selection import LigandSelect, PocketSelect
import numpy as np
import os

def extract_ligand_and_pocket(pdb_file, ligand_file, pocket_file, ligand_names=None, 
                             model_id=None, chain_id=None, multi_ligand=False, 
                             radius=10.0, ext=".pdb"):
    """Extract ligands and corresponding binding pockets from a PDB file."""
    
    # Load the structure
    structure = PDB.PDBParser(QUIET=True).get_structure("STRUCTURE", pdb_file)

    # Extract ligands
    ligand_select = LigandSelect(ligand_names=ligand_names, model_id=model_id, chain_id=chain_id)
    ligand_structures = []
    for model in structure:
        if not ligand_select.accept_model(model):
            continue
        for chain in model:
            if not ligand_select.accept_chain(chain):
                continue
            for residue in chain.get_unpacked_list():
                if ligand_select.accept_residue(residue):
                    ligand_structures.append(residue)

    if not ligand_structures:
        raise ValueError("No ligands found matching the given criteria.")

    if multi_ligand:
        # Parse output paths
        lig_dir, lig_filename = os.path.split(ligand_file)
        pocket_dir, pocket_filename = os.path.split(pocket_file)
        
        # Get base names (use default if directory)
        lig_base = os.path.splitext(lig_filename)[0] if lig_filename else "ligand"
        pocket_base = os.path.splitext(pocket_filename)[0] if pocket_filename else f"pocket"

        # Create output directories
        os.makedirs(lig_dir, exist_ok=True) if lig_dir else None
        os.makedirs(pocket_dir, exist_ok=True) if pocket_dir else None

        # Process each ligand
        for i, lig_res in enumerate(ligand_structures, 1):
            # Generate ligand filename
            ligand_name = lig_res.get_resname().strip()
            lig_output = os.path.join(lig_dir, f"{ligand_name}_{lig_base}_{i}{ext}")
            save_structure(lig_output, lig_res, ligand_select, ext)

            # Generate pocket filename
            pocket_output = os.path.join(pocket_dir, f"{ligand_name}_{pocket_base}_{i}{ext}")
            
            # Extract pocket
            lig_coords = np.array([atom.coord for atom in lig_res.get_atoms()])
            pocket_select = PocketSelect(radius=radius, ligand_coords=lig_coords)
            save_structure(pocket_output, structure, pocket_select, ext)
        return len(ligand_structures)
    else:  
        # Save merged ligands
        save_structure(ligand_file, structure, ligand_select, ext)

        # Extract combined pocket
        all_coords = np.array([atom.coord for lig in ligand_structures for atom in lig.get_atoms()])
        pocket_select = PocketSelect(radius=radius, ligand_coords=all_coords)
        save_structure(pocket_file, structure, pocket_select, ext)

        return 1

def main():
    parser = argparse.ArgumentParser(description="Extract ligands and their binding pockets from a PDB file.")
    parser.add_argument("pdb_file", type=str, help="Input PDB file")
    parser.add_argument("-l", "--ligand_file", type=str, default="ligand.pdb", help="Output ligand file")
    parser.add_argument("-p", "--pocket_file", type=str, default="pocket.pdb", help="Output pocket file")
    parser.add_argument("--ligand_names", type=str, nargs="+", default=None, help="Ligand names in PDB file")
    parser.add_argument("--model_id", type=int, default=None, help="Model ID (default: None)")
    parser.add_argument("--chain_id", type=str, default=None, help="Chain ID (default: None)")
    parser.add_argument("--multi_ligand", action="store_true", help="Extract multiple ligands separately")
    parser.add_argument("--radius", type=float, default=10.0, help="Pocket search radius")
    parser.add_argument("--ext", type=str, choices=[".pdb", ".cif"], default=".pdb", help="Output file format")

    args = parser.parse_args()
    num_ligands = extract_ligand_and_pocket(
        args.pdb_file, args.ligand_file, args.pocket_file,
        args.ligand_names, args.model_id, args.chain_id, 
        args.multi_ligand, args.radius, args.ext
    )
    print(f"Extraction complete. Ligands saved to {args.ligand_file} and pockets saved to {args.pocket_file}")
    if args.multi_ligand:
        print(f"Number of pocket-ligand pairs extracted: {num_ligands}")

if __name__ == "__main__":
    main()
