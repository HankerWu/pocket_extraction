import argparse
def str_to_bool(value):
    if value.lower() in {"true", "t", "yes", "y", "1", "True", "T", "Yes", "Y", "1"}:
        return True
    elif value.lower() in {"false", "f", "no", "n", "0", "False", "F", "No", "N", "0"}:
        return False
    else:
        raise argparse.ArgumentTypeError("Boolean value expected.")

def add_common_arguments(parser: argparse.ArgumentParser, include_output: bool = True) -> None:
    """Add arguments common to all scripts."""
    # Input/output
    parser.add_argument("pdb_file", help="Input structure file, or a PDB ID")
    if include_output:
        parser.add_argument("-o", "--output", default="output.pdb",
                            help="Output path (file/directory)")
    parser.add_argument("--ext", choices=["pdb", "cif"],
                      help="Output format override")
    
    # Logging
    parser.add_argument("-q", "--quiet", action="store_true",
                      help="Suppress informational output")
    parser.add_argument("--debug", action="store_true",
                      help="Enable debug logging")
    parser.add_argument("--logfile", help="Path to log file")

def add_ligand_selection_arguments(parser: argparse.ArgumentParser) -> None:
    """Add arguments for ligand selection."""
    parser.add_argument("--ligand_names", nargs="+",
                      help="Specific ligand names to extract")
    parser.add_argument("--exclude", nargs="+",
                      help="Additional ligands to exclude")
    parser.add_argument("-m", "--model_ids", type=int, nargs="+",
                      help="Model IDs to extract from")
    parser.add_argument("-c", "--chain_ids", nargs="+",
                      help="Chain IDs to extract from")
    parser.add_argument("--multi", action="store_true",
                      help="Save separate files per ligand")
    parser.add_argument("--ignore_duplicates", type=str_to_bool, default=True,
                      help="Whether to ignore duplicate ligands in the same file")

def add_pocket_arguments(parser: argparse.ArgumentParser, include_reference: bool = True) -> None:
    """Add arguments for pocket extraction."""
    if include_reference:
        group = parser.add_mutually_exclusive_group(required=True)
        group.add_argument("--ligand_file", help="Ligand structure file")
        group.add_argument("--ligand_center", nargs=3, type=float,
                           help="Manual ligand center coordinates (X Y Z)")
    parser.add_argument("-r", "--radius", type=float, default=10.0,
                      help="Pocket radius in Angstroms")

def get_ligand_parser() -> argparse.ArgumentParser:
    """Create parser for ligand extraction."""
    parser = argparse.ArgumentParser(
        description="Extract ligands from structure files",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    add_common_arguments(parser)
    add_ligand_selection_arguments(parser)
    return parser

def get_pocket_parser() -> argparse.ArgumentParser:
    """Create parser for pocket extraction."""
    parser = argparse.ArgumentParser(
        description="Extract binding pockets from structures",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    add_common_arguments(parser)
    add_pocket_arguments(parser)
    return parser

def get_ligand_and_pocket_parser() -> argparse.ArgumentParser:
    """Create parser for ligand and pocket extraction."""
    parser = argparse.ArgumentParser(
        description="Extract ligands and binding pockets simultaneously",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    # Override output argument for combined extraction
    parser.add_argument("-ol", "--output_ligand", default="ligand.pdb",
                      help="Output ligand path (file/directory)")
    parser.add_argument("-op", "--output_pocket", default="pocket.pdb",
                      help="Output pocket path (file/directory)")
    add_common_arguments(parser, include_output=False)
    add_ligand_selection_arguments(parser)
    add_pocket_arguments(parser, include_reference=False)
        
    return parser
