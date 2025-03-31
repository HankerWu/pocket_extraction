import argparse
import logging
import os
import numpy as np
from typing import Tuple, List, Optional
from pathlib import Path
from Bio import PDB
from .data_utils import load_structure, save_structure, process_output_path
from .selection import LigandSelect, PocketSelect

logger = logging.getLogger(__name__)

def extract_ligand_and_pocket(
    pdb_file: str,
    ligand_path: str,
    pocket_path: str,
    ligand_names: Optional[List] = None,
    model_id: Optional[int] = None,
    chain_id: Optional[str] = None,
    multi_mode: bool = False,
    radius: float = 10.0,
    ext: Optional[str] = None
) -> Tuple[int, str|None, str|None]:
    """Main extraction logic with unified error handling."""
    try:
        structure = load_structure(pdb_file)
        ligand_selector = LigandSelect(ligand_names, model_id, chain_id)
        
        # Find matching ligands
        ligands = []
        for model in structure:
            if not ligand_selector.accept_model(model):
                continue
            for chain in model:
                if not ligand_selector.accept_chain(chain):
                    continue
                ligands.extend(res for res in chain.get_unpacked_list() 
                             if ligand_selector.accept_residue(res))
        
        if not ligands:
            logger.error("No ligands found matching criteria")
            return 0, None, None
        
        # Handle output modes
        if not multi_mode:
            # Single file mode
            lig_file = process_output_path(ligand_path, "ligand", ext)
            save_structure(lig_file, structure, ligand_selector)
            
            # Extract combined pocket
            all_coords = np.array([atom.coord for lig in ligands for atom in lig.get_atoms()])
            pocket_selector = PocketSelect(radius=radius, ligand_coords=all_coords)
            pocket_file = process_output_path(pocket_path, "pocket", ext)
            save_structure(pocket_file, structure, pocket_selector)
            return 1, lig_file, pocket_file
        
        else:
            # Multi-file mode
            count = 0
            for idx, lig in enumerate(ligands, 1):
                # Save ligand
                lig_name = lig.get_resname().strip()
                lig_file = process_output_path(
                    ligand_path, lig_name, ext, idx
                )
                save_structure(lig_file, lig)
                
                # Save corresponding pocket
                pocket_coords = np.array([atom.coord for atom in lig.get_atoms()])
                pocket_selector = PocketSelect(radius, pocket_coords)
                pocket_file = process_output_path(
                    pocket_path, f"{lig_name}_pocket", ext, idx
                )
                save_structure(pocket_file, structure, pocket_selector)
                count += 1
            return count, os.path.dirname(ligand_path), os.path.dirname(pocket_path)
        
    except Exception as e:
        logger.error(f"Extraction failed: {str(e)}")
        raise

def main():
    """CLI interface with unified argument handling."""
    parser = argparse.ArgumentParser(
        description="Extract ligands and corresponding binding pockets"
    )
    parser.add_argument("pdb_file", help="Input structure (PDB/mmCIF)")
    parser.add_argument("-l", "--ligand", default="ligand.pdb",
                      help="Output ligand path (file/directory)")
    parser.add_argument("-p", "--pocket", default="pocket.pdb",
                      help="Output pocket path (file/directory)")
    parser.add_argument("--ligands", nargs="+",
                      help="Target ligand names (e.g. ATP PO4)")
    parser.add_argument("--model", type=int,
                      help="Specific model ID")
    parser.add_argument("--chain", type=str,
                      help="Specific chain ID")
    parser.add_argument("--multi", action="store_true",
                      help="Separate files per ligand-pocket pair")
    parser.add_argument("-r", "--radius", type=float, default=10.0,
                      help="Pocket radius in Å (default: 10.0)")
    parser.add_argument("--ext", choices=["pdb", "cif"],
                      help="Output format override")
    
    args = parser.parse_args()
    
    try:
        count, lig_path, pocket_path = extract_ligand_and_pocket(
            args.pdb_file,
            args.ligand,
            args.pocket,
            args.ligands,
            args.model,
            args.chain,
            args.multi,
            args.radius,
            args.ext
        )
        
        if args.multi:
            print(f"Extracted {count} ligand-pocket pairs to directories:")
            print(f"Ligands: {lig_path if lig_path else '.'}")
            print(f"Pockets: {pocket_path if pocket_path else '.'}")
        else:
            print(f"Successfully extracted ligand-pocket pair to:")
            print(f"Ligand: {lig_path if lig_path else './'}")
            print(f"Pocket: {pocket_path if pocket_path else './'}")
            
    except Exception as e:
        logger.error(f"Fatal error: {str(e)}")
        exit(1)

if __name__ == "__main__":
    main()
    