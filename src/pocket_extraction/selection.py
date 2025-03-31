from Bio.PDB import Select
from Bio import PDB
import numpy as np
from typing import Optional, List
from .data_utils import (
    residue_name_to_one_letter_code,
    amino_acid_bounding_radii,
    get_remove_ligands,
)
remove_ligands = get_remove_ligands()

class PocketSelect(Select):
    """Efficiently selects binding pocket residues using bounding sphere optimization.
    
    Supports two ligand positioning modes:
    1. Ligand center point (fast)
    2. Full ligand atom coordinates (precise)
    """
    
    def __init__(
        self,
        radius: float,
        ligand_coords: Optional[np.ndarray] = None,
        ligand_center: Optional[np.ndarray] = None
    ):
        """
        Args:
            radius: Search radius in Angstroms
            ligand_coords: Ligand atom coordinates (Nx3 numpy array)
            ligand_center: Precomputed ligand center (3-element array)
            Must provide exactly one of ligand_coords or ligand_center
        """
        assert (ligand_coords is not None) ^ (ligand_center is not None), \
            "Exactly one of ligand_coords or ligand_center must be provided"
            
        self.radius = radius
        
        # Process ligand spatial parameters
        if ligand_coords is not None:
            self.ligand_coords = np.array(ligand_coords)
            self.ligand_center = np.mean(ligand_coords, axis=0)
            # Calculate ligand bounding sphere radius
            self.ligand_radius = np.max(
                np.linalg.norm(ligand_coords - self.ligand_center, axis=1)
            )
        else:
            self.ligand_center = np.array(ligand_center)
            self.ligand_radius = 0.0  # Not used in center mode
            
        # Precompute extended radius for fast rejection
        self.extended_radius = self.radius + self.ligand_radius

    def accept_residue(self, residue):
        """Main selection logic for PDB residues"""
        # Skip non-standard amino acids
        resname = residue.get_resname().strip()
        one_letter = residue_name_to_one_letter_code.get(resname)
        if one_letter is None:
            return False
            
        # Get residue spatial parameters
        residue_radius = amino_acid_bounding_radii.get(one_letter, 0.0)
        residue_center = self._get_residue_center(residue)
        if residue_center is None:
            return False
            
        # Fast sphere-sphere collision check
        center_distance = np.linalg.norm(residue_center - self.ligand_center)
        if center_distance > 1.5 * (self.extended_radius + residue_radius):
            return False
            
        # Detailed atom-level check
        return self._check_atoms(residue)

    def _get_residue_center(self, residue):
        """Get approximate residue center using CA atom or first available atom"""
        try:
            return residue["CA"].coord
        except KeyError:
            try:
                return next(residue.get_atoms()).coord
            except StopIteration:
                return None

    def _check_atoms(self, residue):
        """Dispatch to appropriate checking method based on ligand input type"""
        if hasattr(self, 'ligand_coords'):
            return self._check_ligand_coords(residue)
        return self._check_ligand_center(residue)

    def _check_ligand_center(self, residue):
        """Check against ligand center with combined radii"""
        for atom in residue:
            distance = np.linalg.norm(atom.coord - self.ligand_center)
            if distance <= self.radius:
                return True
        return False

    def _check_ligand_coords(self, residue):
        """Vectorized check against all ligand atoms using numpy"""
        for atom in residue:
            # Vectorized distance calculation to all ligand atoms
            distances = np.linalg.norm(
                self.ligand_coords - atom.coord,
                axis=1
            )
            if np.any(distances <= self.radius):
                return True
        return False
    
class LigandSelect(Select):
    """Selects only the ligand atoms from a structure"""
    def __init__(
            self, 
            ligand_names: Optional[List[str]] = None, 
            model_id: Optional[int] = None, 
            chain_id: Optional[str] = None
        ):
        self.ligand_names = ligand_names
        if self.ligand_names is not None:
            for ligand_name in self.ligand_names:
                assert ligand_name not in remove_ligands, \
                    f"Ligand name {ligand_name} is in the remove_ligands list, please select another ligand"
        self.model_id = model_id
        self.chain_id = chain_id
        
    def accept_model(self, model):
        if self.model_id is None:
            return True
        return True if model.id == self.model_id else False

    def accept_chain(self, chain):
        if self.chain_id is None:
            return True
        return True if chain.id == self.chain_id else False
    
    def accept_residue(self, residue):
        """Accept only ligand residues"""
        residue_name = residue.get_resname().strip()
        if self.ligand_names is not None:
            if residue_name in self.ligand_names:
                return True
            else:
                return False
        if not residue.get_id()[0].startswith("H_"):
            return False
        if residue_name in remove_ligands:
            print(f"Remove ligand {residue_name} from the structure")
            return False
        return True
