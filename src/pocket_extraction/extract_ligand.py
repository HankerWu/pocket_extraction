import argparse
import os
from typing import Tuple, List, Optional
import logging
from Bio import PDB
from .data_utils import load_structure, save_structure, process_output_path
from .selection import LigandSelect

logger = logging.getLogger(__name__)

def extract_ligand(
    pdb_file: str,
    output_path: str,
    ligand_names: Optional[List] = None,
    multi_mode: bool = False,
    model_id: Optional[int] = None,
    chain_id: Optional[str] = None,
    ext: Optional[str] = None
) -> Tuple[int, str|None]:
    """Main ligand extraction logic with exclusion support."""
    try:
        structure = load_structure(pdb_file)
        selector = LigandSelect(
            ligand_names=ligand_names,
            model_id=model_id,
            chain_id=chain_id
        )
        
        # Collect matching ligands
        ligands = []
        for model in structure:
            if not selector.accept_model(model):
                continue
            for chain in model:
                if not selector.accept_chain(chain):
                    continue
                for residue in chain.get_unpacked_list():
                    if selector.accept_residue(residue):
                        ligands.append(residue)
        
        if not ligands:
            logger.warning("No ligands found matching criteria")
            return 0, None
        
        # Handle output
        if not multi_mode:
            output_file = process_output_path(output_path, "ligand", ext)
            save_structure(output_file, structure, selector)
            return 1, output_file
        else:
            for i, lig in enumerate(ligands, 1):
                output_file = process_output_path(output_path, "ligand", ext, i)
                save_structure(output_file, lig)
            return len(ligands), os.path.dirname(output_path)
        
    except Exception as e:
        logger.error(f"Extraction failed: {str(e)}")
        raise

def main():
    """CLI entry point with improved error handling."""
    parser = argparse.ArgumentParser(description="Extract ligands from structure files")
    parser.add_argument("pdb_file", help="Input complex structure (PDB/mmCIF)")
    parser.add_argument("-o", "--output", default="ligand.pdb",
                      help="Output path (file/directory)")
    parser.add_argument("--ligands", nargs="+", 
                      help="Include specific ligand names")
    parser.add_argument("--multi", action="store_true",
                      help="Save each ligand in separate files")
    parser.add_argument("--model", type=int,
                      help="Select specific model")
    parser.add_argument("--chain", type=str,
                      help="Select specific chain")
    parser.add_argument("--ext", choices=["pdb", "cif"], 
                      help="Force output format (default: auto from filename)")
    
    args = parser.parse_args()
    
    try:
        count, output = extract_ligand(
            pdb_file=args.pdb_file,
            output_path=args.output,
            ligand_names=args.ligands,
            multi_mode=args.multi,
            model_id=args.model,
            chain_id=args.chain,
            ext=args.ext
        )
        
        print(f"Successfully extracted {count} ligands and saved to {output if output else '.'}")
            
    except Exception as e:
        logger.error(f"Fatal error: {str(e)}")
        exit(1)

if __name__ == "__main__":
    main()
    