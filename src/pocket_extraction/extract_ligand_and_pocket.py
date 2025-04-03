import argparse
import numpy as np
from typing import Optional, List
from pathlib import Path
from .data_utils import load_structure, save_structure, process_output_path
from .selection import LigandSelect, PocketSelect, REMOVE_LIGANDS
from .logger import logger, setup_logger

def extract_ligand_and_pocket(
    pdb_file: str,
    ligand_path: str,
    pocket_path: str,
    ligand_names: Optional[List[str]] = None,
    model_id: Optional[int] = None,
    chain_id: Optional[str] = None,
    multi_mode: bool = False,
    radius: float = 10.0,
    ext: Optional[str] = None,
    quiet: bool = False
) -> int:
    """Simultaneous extraction of ligands and binding pockets."""
    try:
        structure = load_structure(pdb_file, quiet)
        ligand_selector = LigandSelect(
            ligand_names=ligand_names,
            model_id=model_id,
            chain_id=chain_id,
            quiet=quiet
        )
        
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
            logger.warning("No matching ligands found")
            return 0
        
        # Process based on mode
        count = 0
        if not multi_mode:
            # Combined output
            lig_file = process_output_path(ligand_path, "ligand", ext)
            save_structure(lig_file, structure, ligand_selector, quiet)
            
            # Combined pocket
            all_coords = np.concatenate([np.array([atom.coord for atom in lig.get_atoms()]) 
                                       for lig in ligands])
            pocket_selector = PocketSelect(
                radius=radius,
                ligand_coords=all_coords,
                quiet=quiet
            )
            pocket_file = process_output_path(pocket_path, "pocket", ext)
            save_structure(pocket_file, structure, pocket_selector, quiet)
            
            count = len(ligands)
        else:
            # Per-ligand output
            for idx, lig in enumerate(ligands, 1):
                # Save ligand
                lig_name = lig.get_resname().strip()
                lig_file = process_output_path(ligand_path, lig_name, ext, idx)
                save_structure(lig_file, lig, quiet=quiet)
                
                # Save pocket
                coords = np.array([atom.coord for atom in lig.get_atoms()])
                pocket_selector = PocketSelect(
                    radius=radius,
                    ligand_coords=coords,
                    quiet=quiet
                )
                pocket_file = process_output_path(
                    pocket_path,
                    f"{lig_name}_pocket",
                    ext,
                    idx
                )
                save_structure(pocket_file, structure, pocket_selector, quiet)
                
                count += 1
        
        logger.info(f"Processed {count} ligand-pocket pairs")
        return count
        
    except Exception as e:
        logger.exception("Combined extraction failed")
        raise

def main():
    """CLI for combined ligand and pocket extraction."""
    parser = argparse.ArgumentParser(
        description="Extract ligands and binding pockets simultaneously",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    # Input/output
    parser.add_argument("pdb_file", help="Input structure file")
    parser.add_argument("-l", "--ligand", default="ligand.pdb",
                      help="Output ligand path (file/directory)")
    parser.add_argument("-p", "--pocket", default="pocket.pdb",
                      help="Output pocket path (file/directory)")
    
    # Selection criteria
    parser.add_argument("--ligands", nargs="+",
                      help="Specific ligand names to extract")
    parser.add_argument("--exclude", nargs="+",
                      help="Additional ligands to exclude")
    parser.add_argument("--model", type=int,
                      help="Model ID to extract from")
    parser.add_argument("--chain",
                      help="Chain ID to extract from")
    
    # Modes/parameters
    parser.add_argument("--multi", action="store_true",
                      help="Save separate files per ligand-pocket pair")
    parser.add_argument("-r", "--radius", type=float, default=10.0,
                      help="Pocket radius in Angstroms")
    parser.add_argument("--ext", choices=["pdb", "cif"],
                      help="Output format override")
    
    # Logging
    parser.add_argument("-q", "--quiet", action="store_true",
                      help="Suppress informational output")
    parser.add_argument("--debug", action="store_true",
                      help="Enable debug logging")
    parser.add_argument("--logfile",
                      help="Path to log file")
    
    args = parser.parse_args()
    
    try:
        # Configure logging
        setup_logger(
            quiet=args.quiet,
            debug=args.debug,
            logfile=args.logfile
        )
        
        if args.exclude:
            logger.info(f"Excluding ligands: {args.exclude}")
            for ligand in args.exclude:
                REMOVE_LIGANDS.add(ligand)
        
        # Perform extraction
        count = extract_ligand_and_pocket(
            pdb_file=args.pdb_file,
            ligand_path=args.ligand,
            pocket_path=args.pocket,
            ligand_names=args.ligands,
            model_id=args.model,
            chain_id=args.chain,
            multi_mode=args.multi,
            radius=args.radius,
            ext=args.ext,
            quiet=args.quiet
        )
        
        if not args.quiet:
            logger.info(f"Successfully processed {count} ligand-pocket pairs")
            
    except Exception as e:
        logger.exception("Fatal error")
        exit(1)

if __name__ == "__main__":
    main()
    