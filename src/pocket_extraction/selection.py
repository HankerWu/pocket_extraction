from Bio.PDB import Select
import numpy as np
from pathlib import Path
from typing import Optional, List
from .logger import logger
from .data_utils import RESIDUE_MAP, get_remove_ligands, save_structure

# Immutable default exclusions; additional exclusions belong to each selector.
REMOVE_LIGANDS = frozenset(get_remove_ligands())
_ATOM_PAIR_BLOCK_SIZE = 4096
_MIN_BOUNDING_BOX_LIGAND_ATOMS = 16

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
            coords = np.asarray(ligand_coords, dtype=np.float32)
            if (coords.ndim != 2 or coords.shape[1] != 3
                    or not len(coords) or not np.isfinite(coords).all()):
                raise ValueError("Ligand coordinates must be a nonempty finite Nx3 array")
            self._ligand_bounds = None
            if len(coords) >= _MIN_BOUNDING_BOX_LIGAND_ATOMS:
                coords = coords.copy()
                self._ligand_bounds = (coords.min(axis=0), coords.max(axis=0))
            self._ligand_blocks = tuple(
                coords[start:start + _ATOM_PAIR_BLOCK_SIZE]
                for start in range(0, len(coords), _ATOM_PAIR_BLOCK_SIZE)
            )
            self._atoms_per_block = _ATOM_PAIR_BLOCK_SIZE // len(self._ligand_blocks[0])
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
            
        return (self._check_ligand_coords(residue) if hasattr(self, '_ligand_blocks')
                else self._check_ligand_center(residue))

    def accept_chain(self, chain):
        # PDBIO checks chain IDs before applying accept_residue. Skip an
        # out-of-pocket mmCIF chain even when its ID cannot be written to PDB.
        return len(chain.id) <= 1 or any(map(self.accept_residue, chain.get_unpacked_list()))

    def _check_ligand_center(self, residue):
        """Check against ligand centroid."""
        for atom in residue.get_unpacked_list():
            if np.linalg.norm(atom.coord - self.ligand_center) <= self.radius:
                return True
        return False

    def _check_ligand_coords(self, residue):
        """Check against all ligand atoms."""
        atoms = residue.get_unpacked_list()
        atoms_per_block = self._atoms_per_block
        for start in range(0, len(atoms), atoms_per_block):
            points = np.asarray([atom.coord for atom in atoms[start:start + atoms_per_block]])
            if self._ligand_bounds is not None:
                lower, upper = self._ligand_bounds
                # This lower bound uses the same dtype and norm as the exact check.
                nearest = np.maximum(np.maximum(lower - points, points - upper), 0)
                if not np.any(np.linalg.norm(nearest, axis=-1) <= self.radius):
                    continue
            for ligands in self._ligand_blocks:
                deltas = ligands - points[:, None, :]
                if np.any(np.linalg.norm(deltas, axis=-1) <= self.radius):
                    return True
        return False

class LigandPocketSelect(Select):
    """Select pocket residues against ligands from their own structure model."""

    def __init__(self, ligands, radius, quiet=False):
        coordinates = {}
        for ligand in ligands:
            model_id = ligand.get_parent().get_parent().id
            coordinates.setdefault(model_id, []).extend(
                atom.coord for atom in ligand.get_unpacked_list()
            )
        self.selectors = {
            model_id: PocketSelect(radius, ligand_coords=np.asarray(points), quiet=quiet)
            for model_id, points in coordinates.items()
        }

    def accept_model(self, model):
        return model.id in self.selectors

    def accept_chain(self, chain):
        return self.selectors[chain.get_parent().id].accept_chain(chain)

    def accept_residue(self, residue):
        selector = self.selectors.get(residue.get_parent().get_parent().id)
        return selector is not None and selector.accept_residue(residue)


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
        self.selected_residues = None
        self.selected_chains = None
        self.selected_models = None
        if selected_residues is not None:
            self.selected_residues = set()
            self.selected_chains = set()
            self.selected_models = set()
            for residue in selected_residues:
                chain = residue.get_parent()
                self.selected_residues.add(id(residue))
                self.selected_chains.add(id(chain))
                self.selected_models.add(id(chain.get_parent()))
        
        if self.ligand_names:
            invalid = self.ligand_names & REMOVE_LIGANDS
            if invalid and not quiet:
                logger.warning(f"Excluded ligands in request: {invalid}")

    def accept_model(self, model):
        # PDBIO writes MODEL before looking at its chains. Only models with
        # selected residues should reach the writer.
        if self.selected_models is not None and id(model) not in self.selected_models:
            return False
        return self.model_ids is None or model.id in self.model_ids

    def accept_chain(self, chain):
        if self.selected_chains is not None and id(chain) not in self.selected_chains:
            return False
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
                if ignore_duplicates:
                    name = residue.get_resname().strip()
                    if name in seen:
                        continue
                    seen.add(name)
                ligands.append(residue)
    return ligands


def save_pocket(path, structure, selector, quiet):
    if len(structure) > 1 and Path(path).resolve().suffix.lower() == '.pdb':
        # PDBIO emits MODEL before checking its contents. Select the contacts
        # once so only nonempty models and chains reach the writer.
        contacts = (
            residue for model in structure if selector.accept_model(model)
            for chain in model for residue in chain.get_unpacked_list()
            if selector.accept_residue(residue)
        )
        selector = LigandSelect(selected_residues=contacts, quiet=quiet)
    save_structure(path, structure, selector, quiet)
