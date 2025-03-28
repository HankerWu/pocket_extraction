import argparse
from Bio import PDB
from .data_utils import save_structure
from .selection import LigandSelect
import os

def extract_ligand(pdb_file, ligand_file, ligand_names=None, multi_ligand=False, model_id=None, chain_id=None, ext=".pdb"):
    """Extract ligands from a PDB file and save them to a new file."""
    
    structure = PDB.PDBParser(QUIET=True).get_structure("LIGAND", pdb_file)
    select = LigandSelect(ligand_names=ligand_names, model_id=model_id, chain_id=chain_id)

    if not multi_ligand:
        save_structure(ligand_file, structure, select, ext)
        return 1
    else:
        # Handle multi-ligand output path
        output_dir, filename = os.path.split(ligand_file)
        base_name = "ligand"  # Default base name
        
        # Extract base name if user provided a file path
        if filename:
            base_name = os.path.splitext(filename)[0]
        
        # Create output directory if needed
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)

        # Collect ligand residues
        ligand_structures = []
        for model in structure:
            if not select.accept_model(model): continue
            for chain in model:
                if not select.accept_chain(chain): continue
                for residue in chain.get_unpacked_list():
                    if select.accept_residue(residue):
                        ligand_structures.append(residue)

        # Save individual ligand files
        for i, lig in enumerate(ligand_structures, 1):
            output_path = os.path.join(output_dir, f"{base_name}_{i}{ext}")
            save_structure(output_path, lig, select, ext)
            
        return len(ligand_structures)

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
    
    # Extract ligands
    number_of_ligands = extract_ligand(
        pdb_file=args.pdb_file,
        ligand_file=args.ligand_file,
        ligand_names=args.ligand_names,
        multi_ligand=args.multi_ligand,
        model_id=args.model_id,
        chain_id=args.chain_id,
        ext=args.ext
    )
    print(f"Ligands extracted from {args.pdb_file} and saved to {args.ligand_file}")
    if args.multi_ligand:
        print(f"Number of ligands extracted: {number_of_ligands}")
    else:
        print("Single ligand extracted.")
    
if __name__ == "__main__":
    main()
    