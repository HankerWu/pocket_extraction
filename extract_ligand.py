import argparse
from Bio import PDB
from data_utils import save_structure
from selection import LigandSelect
import os

def main():
    parser = argparse.ArgumentParser(description="Extract ligand from PDB file")
    parser.add_argument("pdb_file", type=str, help="Input PDB file")
    parser.add_argument("-o", "--ligand_file", type=str, default="ligand.pdb", help="Output ligand file")
    parser.add_argument("--multi_ligand", action="store_true", help="Extract multiple ligands")
    parser.add_argument("--model_id", type=int, default=None, help="Model ID (default: None)")
    parser.add_argument("--chain_id", type=str, default=None, help="Chain ID (default: None)")
    parser.add_argument("--ligand_names", type=str, nargs="+", default=None, help="Ligand names (default: None)")
    parser.add_argument("--ext", type=str, default=".pdb", choices=[".pdb", ".cif"], help="File extension (default: .pdb)")
    args = parser.parse_args()

    # Load the structure
    structure = PDB.PDBParser(QUIET=True).get_structure("LIGAND", args.pdb_file)

    # Create a selection object for the ligand
    select = LigandSelect(
        ligand_names=args.ligand_names,
        model_id=args.model_id,
        chain_id=args.chain_id
    )
    if not args.multi_ligand:
        # Save the ligand structure
        save_structure(args.ligand_file, structure, select, args.ext)
        print(f"Ligand saved to {args.ligand_file}")
    else:
        ligand_file = args.ligand_file
        ligand_file_dir = os.path.dirname(ligand_file)
        ligand_file_basename = os.path.basename(ligand_file)
        
        ligand_structures = []
        if len(ligand_file_dir) > 0:
            os.makedirs(ligand_file_dir, exist_ok=True)
        for model in structure.get_list():
            if not select.accept_model(model):
                continue
            for chain in model.get_list():
                if not select.accept_chain(chain):
                    continue
                for residue in chain.get_unpacked_list():
                    if not select.accept_residue(residue):
                        continue
                    ligand_structures.append(residue)
        for i, ligand_structure in enumerate(ligand_structures):
            ligand_filename = os.path.join(ligand_file_dir, f"{ligand_file_basename}_{i+1}{args.ext}")
            save_structure(ligand_filename, ligand_structure, select, args.ext)
    
if __name__ == "__main__":
    main()
    