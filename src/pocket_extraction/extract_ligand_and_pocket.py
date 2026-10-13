from typing import Optional, List
from pathlib import Path
from .data_utils import load_structure, save_structure, process_output_path
from .selection import LigandSelect, LigandPocketSelect, collect_ligands, validate_radius
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
    ignore_duplicates: bool = True,
    exclude: Optional[List[str]] = None
) -> int:
    """Simultaneous extraction of ligands and binding pockets."""
    radius = validate_radius(radius)
    try:
        structure = load_structure(pdb_file, quiet)
        ligand_selector = LigandSelect(
            ligand_names=ligand_names,
            model_ids=model_ids,
            chain_ids=chain_ids,
            quiet=quiet,
            exclude=exclude
        )
        
        ligands = collect_ligands(structure, ligand_selector, ignore_duplicates)

        if not ligands:
            if not quiet:
                logger.warning(f"No matching ligands found for {pdb_file}.")
            return 0
        
        # Validate every geometry and output path before publishing the first file.
        combined = not multi or len(ligands) == 1
        groups = [ligands] if combined else [[ligand] for ligand in ligands]
        outputs = []
        for index, group in enumerate(groups):
            pocket_selector = LigandPocketSelect(group, radius=radius, quiet=quiet)
            if combined:
                ligand_file = process_output_path(output_ligand, "ligand", ext)
                pocket_file = process_output_path(output_pocket, "pocket", ext)
                ligand_source = structure
                selection = LigandSelect(selected_residues=group, quiet=quiet)
            else:
                ligand = group[0]
                name = f"{ligand.get_resname()}_{ligand.get_parent().id}"
                ligand_file = process_output_path(output_ligand, name, ext, index)
                pocket_file = process_output_path(output_pocket, f"{name}_pocket", ext, index)
                ligand_source, selection = ligand, None
            outputs.append((ligand_file, pocket_file, ligand_source, selection, pocket_selector))
        paths = [Path(path).resolve() for output in outputs for path in output[:2]]
        if len(set(paths)) != len(paths):
            raise ValueError("Ligand and pocket outputs must use distinct file paths")
        for ligand_file, pocket_file, ligand_source, selection, pocket_selector in outputs:
            if selection is None:
                save_structure(ligand_file, ligand_source, quiet=quiet)
            else:
                save_structure(ligand_file, ligand_source, selection, quiet)
            save_structure(pocket_file, structure, pocket_selector, quiet)
        count = len(ligands)

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
            ignore_duplicates=args.ignore_duplicates,
            exclude=args.exclude
        )
        
        if not args.quiet:
            logger.info(f"Successfully processed {count} ligand-pocket pairs to: {args.output_ligand} and {args.output_pocket}")
            
    except Exception as e:
        logger.exception("Fatal error")
        exit(1)

if __name__ == "__main__":
    main()
