from Bio.PDB import Select
import numpy as np
from typing import Optional, List, Dict
from .logger import logger
from .data_utils import RESIDUE_MAP, get_remove_ligands

# Immutable default exclusions; additional exclusions belong to each selector.
REMOVE_LIGANDS = frozenset(get_remove_ligands())

def validate_radius(radius):
    """Use the same finite, positive distance contract for both pocket APIs."""
    radius = float(radius)
    if not np.isfinite(radius) or radius <= 0:
        raise ValueError("Radius must be finite and positive")
    return radius

class PocketSelect(Select):
    """Advanced binding pocket selection based on spatial criteria."""
    
    def __init__(
        self,
        radius: float,
        ligand_coords: Optional[np.ndarray] = None,
        ligand_center: Optional[np.ndarray] = None,
        quiet: bool = False
    ):
        """Initialize pocket selector.
        
        Args:
            radius: Search radius in Angstroms
            ligand_coords: Nx3 array of ligand atom coordinates
            ligand_center: Precomputed ligand centroid [x,y,z]
            quiet: Suppress info messages
        """
        if (ligand_coords is None) == (ligand_center is None):
            raise ValueError("Provide exactly one of ligand_coords or ligand_center")
            
        self.radius = validate_radius(radius)
        self.quiet = quiet
        
        if ligand_coords is not None:
            self.ligand_coords = np.asarray(ligand_coords, dtype=np.float32)
            if (self.ligand_coords.ndim != 2 or self.ligand_coords.shape[1] != 3
                    or not len(self.ligand_coords) or not np.isfinite(self.ligand_coords).all()):
                raise ValueError("Ligand coordinates must be a nonempty finite Nx3 array")
        else:
            self.ligand_center = np.asarray(ligand_center, dtype=np.float32)
            if self.ligand_center.shape != (3,) or not np.isfinite(self.ligand_center).all():
                raise ValueError("Ligand center must contain three finite coordinates")

    def accept_residue(self, residue):
        """Select residues by their actual atom distances, including long side chains."""
        resname = residue.get_resname().strip()
        one_letter = RESIDUE_MAP.get(resname)
        if one_letter is None:
            return False
            
        return (self._check_ligand_coords(residue) if hasattr(self, 'ligand_coords') 
                else self._check_ligand_center(residue))

    def _check_ligand_center(self, residue):
        """Check against ligand centroid."""
        for atom in residue:
            if np.linalg.norm(atom.coord - self.ligand_center) <= self.radius:
                return True
        return False

    def _check_ligand_coords(self, residue):
        """Check against all ligand atoms."""
        for atom in residue:
            if np.any(np.linalg.norm(
                self.ligand_coords - atom.coord,
                axis=1
            ) <= self.radius):
                return True
        return False

class LigandSelect(Select):
    """Precise ligand selection with filtering options."""
    
    def __init__(
        self,
        ligand_names: Optional[List[str]] = None,
        model_ids: Optional[List[int]] = None,
        chain_ids: Optional[List[str]] = None,
        quiet: bool = False,
        selected_residues=None,
        exclude: Optional[List[str]] = None
    ):
        """Initialize ligand selector.
        
        Args:
            ligand_names: Specific ligands to extract (None for all)
            model_id: Model ID filter
            chain_id: Chain ID filter
            quiet: Suppress warnings
        """
        self.ligand_names = None if ligand_names is None else set(ligand_names)
        self.model_ids = model_ids
        self.chain_ids = chain_ids
        self.quiet = quiet
        self.exclude = set(exclude or ())
        self.selected_residues = (None if selected_residues is None
                                  else {id(residue) for residue in selected_residues})
        
        if self.ligand_names:
            invalid = self.ligand_names & REMOVE_LIGANDS
            if invalid and not quiet:
                logger.warning(f"Excluded ligands in request: {invalid}")

    def accept_model(self, model):
        if self.model_ids is None:
            return True
        return True if model.id in self.model_ids else False

    def accept_chain(self, chain):
        if self.chain_ids is None:
            return True
        return True if chain.id in self.chain_ids else False
    
    def accept_residue(self, residue):
        if self.selected_residues is not None:
            return id(residue) in self.selected_residues
        resname = residue.get_resname().strip()
        if resname in self.exclude:
            return False
        
        if self.ligand_names is not None:
            return resname in self.ligand_names
            
        if not residue.id[0].startswith("H_"):  # Skip standard residues
            return False
            
        return resname not in REMOVE_LIGANDS

def collect_ligands(structure, selector, ignore_duplicates=True):
    """Collect filtered residues once; optionally keep the first of each name."""
    ligands = []
    seen = set()
    for model in structure:
        if not selector.accept_model(model):
            continue
        for chain in model:
            if not selector.accept_chain(chain):
                continue
            for residue in chain.get_unpacked_list():
                if not selector.accept_residue(residue):
                    continue
                name = residue.get_resname().strip()
                if ignore_duplicates and name in seen:
                    continue
                seen.add(name)
                ligands.append(residue)
    return ligands
