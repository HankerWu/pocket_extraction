from typing import Optional, List
from pathlib import Path
from .data_utils import load_structure, save_structure, process_output_path
from .selection import LigandSelect, REMOVE_LIGANDS
from .logger import logger, setup_logger
from .arguments import get_ligand_parser

def extract_ligand(
    pdb_file: str,
    output_path: str,
    ligand_names: Optional[List[str]] = None,
    multi: bool = False,
    model_ids: Optional[List[int]] = None,
    chain_ids: Optional[List[str]] = None,
    ext: Optional[str] = None,
    quiet: bool = False,
    ignore_duplicates: bool = True
) -> int:
    """Extract ligands from structure with flexible filtering."""
    output_is_dir = output_path.endswith(("/", "\\")) or len(Path(output_path).suffix) <= 1
    
    try:
        structure = load_structure(pdb_file, quiet)
        selector = LigandSelect(
            ligand_names=ligand_names,
            model_ids=model_ids,
            chain_ids=chain_ids,
            quiet=quiet
        )
        
        # Collect matching ligands
        ligands = []
        for model in structure:
            if not selector.accept_model(model):
                continue
            for chain in model:
                if not selector.accept_chain(chain):
                    continue
                for res in chain.get_unpacked_list():
                    if selector.accept_residue(res):
                        ligands.append(res)
        
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
        
        # Save results
        count = 0
        if not multi or len(ligands) == 1:
            out_path = process_output_path(
                output_path,
                "ligand" if output_is_dir else Path(output_path).stem,
                ext
            )
            save_structure(out_path, ligands[0], quiet=quiet)
            count = len(ligands)
        else:
            for i, lig in enumerate(ligands):
                out_path = process_output_path(
                    output_path,
                    f"{lig.get_resname()}_{lig.get_parent().id}",
                    ext,
                    i
                )
                logger.debug(f"Saving ligand {i}: {out_path}, {lig.get_resname()}, {lig}")
                save_structure(out_path, lig, quiet=quiet)
                count += 1
        
        if not quiet:
            logger.info(f"Extracted {count} ligand(s) from {pdb_file} to: {output_path}")
        return count
        
    except Exception as e:
        logger.exception(f"Ligand extraction failed for {pdb_file}: {e}")
        raise

def main():
    parser = get_ligand_parser()
    
    args = parser.parse_args()
    
    try:
        # Configure logging
        setup_logger(
            quiet=args.quiet,
            debug=args.debug,
            logfile=args.logfile
        )
        
        if args.exclude:
            logger.warning(f"Excluding ligands: {args.exclude}")
            for ligand in args.exclude:
                REMOVE_LIGANDS.add(ligand)
        
        # Perform extraction
        count = extract_ligand(
            pdb_file=args.pdb_file,
            output_path=args.output,
            ligand_names=args.ligand_names,
            multi=args.multi,
            model_ids=args.model_ids,
            chain_ids=args.chain_ids,
            ext=args.ext,
            quiet=args.quiet,
            ignore_duplicates=args.ignore_duplicates
        )
        
        if not args.quiet:
            logger.info(f"Extracted {count} ligand(s) to: {args.output}")
            
    except Exception as e:
        logger.exception("Fatal error")
        exit(1)

if __name__ == "__main__":
    main()
    