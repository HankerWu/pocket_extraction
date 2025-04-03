from Bio.PDB import Select
import numpy as np
from typing import Optional, List, Dict
from .logger import logger
from .data_utils import RESIDUE_MAP, AA_RADII, get_remove_ligands

# Global exclusion set
REMOVE_LIGANDS = get_remove_ligands()

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
            
        self.radius = float(radius)
        self.quiet = quiet
        
        if ligand_coords is not None:
            self.ligand_coords = np.asarray(ligand_coords, dtype=np.float32)
            self.ligand_center = np.mean(self.ligand_coords, axis=0)
            self.ligand_radius = np.max(np.linalg.norm(
                self.ligand_coords - self.ligand_center, 
                axis=1
            ))
        else:
            self.ligand_center = np.asarray(ligand_center, dtype=np.float32)
            self.ligand_radius = 0.0
            
        self.extended_radius = self.radius + self.ligand_radius

    def accept_residue(self, residue):
        """Two-stage residue selection for optimal performance."""
        resname = residue.get_resname().strip()
        one_letter = RESIDUE_MAP.get(resname)
        if one_letter is None:
            return False
            
        # Stage 1: Fast bounding sphere check
        res_center = self._get_residue_center(residue)
        if res_center is None:
            return False
            
        res_radius = AA_RADII.get(one_letter, 2.0)
        dist = np.linalg.norm(res_center - self.ligand_center)
        
        if dist > 1.5 * (self.extended_radius + res_radius):
            return False
            
        # Stage 2: Precise atom-level check
        return (self._check_ligand_coords(residue) if hasattr(self, 'ligand_coords') 
                else self._check_ligand_center(residue))

    def _get_residue_center(self, residue):
        """Get approximate residue center coordinates."""
        try:
            return residue["CA"].coord
        except KeyError:
            try:
                return next(atom.coord for atom in residue.get_atoms())
            except StopIteration:
                logger.debug(f"No coordinates for {residue.get_resname()}")
                return None

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
        model_id: Optional[int] = None,
        chain_id: Optional[str] = None,
        quiet: bool = False
    ):
        """Initialize ligand selector.
        
        Args:
            ligand_names: Specific ligands to extract (None for all)
            model_id: Model ID filter
            chain_id: Chain ID filter
            quiet: Suppress warnings
        """
        self.ligand_names = None if ligand_names is None else set(ligand_names)
        self.model_id = model_id
        self.chain_id = chain_id
        self.quiet = quiet
        
        if self.ligand_names:
            invalid = self.ligand_names & REMOVE_LIGANDS
            if invalid and not quiet:
                logger.warning(f"Excluded ligands in request: {invalid}")

    def accept_model(self, model):
        if self.model_id is None:
            return True
        return True if model.id == self.model_id else False

    def accept_chain(self, chain):
        if self.chain_id is None:
            return True
        return True if chain.id == self.chain_id else False
    
    def accept_residue(self, residue):
        resname = residue.get_resname().strip()
        
        if self.ligand_names is not None:
            return resname in self.ligand_names
            
        if not residue.id[0].startswith("H_"):  # Skip standard residues
            return False
            
        return resname not in REMOVE_LIGANDS
    