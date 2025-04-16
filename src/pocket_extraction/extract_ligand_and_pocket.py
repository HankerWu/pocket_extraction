import numpy as np
from typing import Optional, List
from pathlib import Path
from .data_utils import load_structure, save_structure, process_output_path
from .selection import LigandSelect, PocketSelect, REMOVE_LIGANDS
from .logger import logger, setup_logger
from .arguments import get_ligand_and_pocket_parser

def extract_ligand_and_pocket(
    pdb_file: str,
    output_ligand: str,
    output_pocket: str,
    ligand_names: Optional[List[str]] = None,
    model_ids: Optional[List[int]] = None,
    chain_ids: Optional[List[str]] = None,
    multi: bool = False,
    radius: float = 10.0,
    ext: Optional[str] = None,
    quiet: bool = False,
    ignore_duplicates: bool = True
) -> int:
    """Simultaneous extraction of ligands and binding pockets."""
    try:
        structure = load_structure(pdb_file, quiet)
        ligand_selector = LigandSelect(
            ligand_names=ligand_names,
            model_ids=model_ids,
            chain_ids=chain_ids,
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
        
        if ignore_duplicates:
            seen = set()
            unique_ligands = []
            for lig in ligands:
                resname = lig.get_resname().strip()
                if resname not in seen:
                    seen.add(resname)
                    unique_ligands.append(lig)
            ligands = unique_ligands
        
        if not ligands:
            if not quiet:
                logger.warning(f"No matching ligands found for {pdb_file}.")
            return 0
        
        # Process based on mode
        count = 0
        if not multi or len(ligands) == 1:
            # Combined output
            lig_file = process_output_path(output_ligand, "ligand", ext)
            save_structure(lig_file, structure, ligand_selector, quiet)
            
            # Combined pocket
            all_coords = np.concatenate([np.array([atom.coord for atom in lig.get_atoms()]) 
                                       for lig in ligands])
            pocket_selector = PocketSelect(
                radius=radius,
                ligand_coords=all_coords,
                quiet=quiet
            )
            pocket_file = process_output_path(output_pocket, "pocket", ext)
            save_structure(pocket_file, structure, pocket_selector, quiet)
            
            count = len(ligands)
        else:
            # Per-ligand output
            for idx, lig in enumerate(ligands):
                # Save ligand
                lig_name = f"{lig.get_resname()}_{lig.get_parent().id}"
                lig_file = process_output_path(output_ligand, lig_name, ext, idx)
                save_structure(lig_file, lig, quiet=quiet)
                
                # Save pocket
                coords = np.array([atom.coord for atom in lig.get_atoms()])
                pocket_selector = PocketSelect(
                    radius=radius,
                    ligand_coords=coords,
                    quiet=quiet
                )
                pocket_file = process_output_path(
                    output_pocket,
                    f"{lig_name}_pocket",
                    ext,
                    idx
                )
                save_structure(pocket_file, structure, pocket_selector, quiet)
                
                count += 1
        
        if not quiet:
            logger.info(f"Processed {count} ligand-pocket pairs from {pdb_file}.")
        return count
        
    except Exception as e:
        logger.exception(f"Extraction failed for {pdb_file}: {e}")
        raise

def main():
    parser = get_ligand_and_pocket_parser()
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
            output_ligand=args.output_ligand,
            output_pocket=args.output_pocket,
            ligand_names=args.ligand_names,
            model_ids=args.model_ids,
            chain_ids=args.chain_ids,
            multi=args.multi,
            radius=args.radius,
            ext=args.ext,
            quiet=args.quiet,
            ignore_duplicates=args.ignore_duplicates
        )
        
        if not args.quiet:
            logger.info(f"Successfully processed {count} ligand-pocket pairs to: {args.output_ligand} and {args.output_pocket}")
            
    except Exception as e:
        logger.exception("Fatal error")
        exit(1)

if __name__ == "__main__":
    main()
    