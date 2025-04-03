import argparse
from typing import Optional, List
from .data_utils import load_structure, save_structure, process_output_path
from .selection import LigandSelect, REMOVE_LIGANDS
from .logger import logger, setup_logger
from pathlib import Path

def extract_ligand(
    pdb_file: str,
    output_path: str,
    ligand_names: Optional[List[str]] = None,
    multi_mode: bool = False,
    model_id: Optional[int] = None,
    chain_id: Optional[str] = None,
    ext: Optional[str] = None,
    quiet: bool = False
) -> int:
    """Extract ligands from structure with flexible filtering."""
    output_is_dir = output_path.endswith(("/", "\\")) or len(Path(output_path).suffix) <= 1
    
    try:
        structure = load_structure(pdb_file, quiet)
        selector = LigandSelect(
            ligand_names=ligand_names,
            model_id=model_id,
            chain_id=chain_id,
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
        
        if not ligands:
            logger.warning("No matching ligands found")
            return 0
        
        # Save results
        count = 0
        if not multi_mode:
            out_path = process_output_path(
                output_path,
                "ligand" if output_is_dir else Path(output_path).stem,
                ext
            )
            save_structure(out_path, structure, selector, quiet)
            count = len(ligands)
        else:
            for i, lig in enumerate(ligands, 1):
                out_path = process_output_path(
                    output_path,
                    lig.get_resname().strip(),
                    ext,
                    i
                )
                logger.debug(f"Saving ligand {i}: {out_path}, {lig.get_resname()}, {lig}")
                save_structure(out_path, lig, quiet=quiet)
                count += 1
        
        logger.info(f"Extracted {count} ligand(s)")
        return count
        
    except Exception as e:
        logger.exception("Ligand extraction failed")
        raise

def main():
    """CLI for ligand extraction."""
    parser = argparse.ArgumentParser(
        description="Extract ligands from structure files",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    # Input/output
    parser.add_argument("pdb_file", help="Input structure file")
    parser.add_argument("-o", "--output", default="ligand.pdb",
                      help="Output path (file/directory)")
    
    # Selection criteria
    parser.add_argument("--ligands", nargs="+",
                      help="Specific ligand names to extract")
    parser.add_argument("--exclude", nargs="+",
                      help="Additional ligands to exclude")
    parser.add_argument("--model", type=int,
                      help="Model ID to extract from")
    parser.add_argument("--chain",
                      help="Chain ID to extract from")
    
    # Modes
    parser.add_argument("--multi", action="store_true",
                      help="Save each ligand separately")
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
            logger.warning(f"Excluding ligands: {args.exclude}")
            for ligand in args.exclude:
                REMOVE_LIGANDS.add(ligand)
        
        # Perform extraction
        count = extract_ligand(
            pdb_file=args.pdb_file,
            output_path=args.output,
            ligand_names=args.ligands,
            multi_mode=args.multi,
            model_id=args.model,
            chain_id=args.chain,
            ext=args.ext,
            quiet=args.quiet
        )
        
        if not args.quiet:
            logger.info(f"Extracted {count} ligand(s) to: {args.output}")
            
    except Exception as e:
        logger.exception("Fatal error")
        exit(1)

if __name__ == "__main__":
    main()
    