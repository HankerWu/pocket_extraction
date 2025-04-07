import argparse
from typing import Union, List, Dict

def remove_arguments(parser: argparse.ArgumentParser, 
                    target_opts: Union[str, List[str]]) -> None:
    """
    Remove arguments while maintaining integrity of multiple groups.
    
    Features:
    - Handles multiple independent groups
    - Preserves group-specific requirements
    - Cross-group removal safety
    
    Args:
        parser: Configured ArgumentParser with potential multiple groups
        target_opts: Option string(s) to remove (supports multi-group targets)
    """
    target_opts = [target_opts] if isinstance(target_opts, str) else target_opts

    # Phase 1: Map each group to its original state
    group_registry: Dict[argparse._MutuallyExclusiveGroup, Dict] = {}
    for group in parser._mutually_exclusive_groups:
        group_registry[group] = {
            'required': group.required,
            'original_actions': set(group._group_actions),
            'remaining_actions': set(group._group_actions)
        }

    # Phase 2: Remove targets across all groups
    target_actions = []
    for action in parser._actions:
        if any(opt in action.option_strings for opt in target_opts):
            target_actions.append(action)

    for action in target_actions:
        # Track parent groups
        parent_groups = set()
        
        # Remove from all containers
        for group in parser._action_groups + parser._mutually_exclusive_groups:
            if action in group._group_actions:
                group._group_actions.remove(action)
                if group in group_registry:
                    group_registry[group]['remaining_actions'].discard(action)
                    parent_groups.add(group)

        # Clean option strings
        for opt in action.option_strings:
            if opt in parser._option_string_actions:
                del parser._option_string_actions[opt]
        parser._actions.remove(action)

    # Phase 3: Reconcile each group independently
    for group in list(parser._mutually_exclusive_groups):
        state = group_registry.get(group, {})
        current_actions = group._group_actions
        
        # Case 1: Empty group - remove completely
        if not current_actions:
            parser._mutually_exclusive_groups.remove(group)
            continue
            
        # Case 2: Single argument - convert to standalone
        if len(current_actions) == 1:
            arg = current_actions[0]
            arg.required = state.get('required', False)
            
            # Migrate to main group
            main_group = parser._action_groups[0]
            if arg not in main_group._group_actions:
                main_group._group_actions.append(arg)
            
            # Remove from mutex groups
            parser._mutually_exclusive_groups.remove(group)
        
        # Case 3: Multiple remain - preserve original group behavior
        else:
            original_size = len(state.get('original_actions', []))
            current_size = len(current_actions)
            
            # Only update requirement if group composition changed
            if current_size < original_size:
                group.required = state.get('required', False)

    # Force help system reset
    parser._get_formatter = lambda: parser.formatter_class(prog=parser.prog)

def add_common_arguments(parser: argparse.ArgumentParser) -> None:
    """Add arguments common to all scripts."""
    # Input/output
    parser.add_argument("pdb_file", help="Input structure file, or a PDB ID")
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

def add_pocket_arguments(parser: argparse.ArgumentParser) -> None:
    """Add arguments for pocket extraction."""
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
    add_common_arguments(parser)
    add_ligand_selection_arguments(parser)
    add_pocket_arguments(parser)
    
    # Remove arguments for combined extraction
    remove_arguments(parser, ["--output", "--ligand_file", "--ligand_center"])
        
    return parser
