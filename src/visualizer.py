"""PoseAI Visualizer Module — Interactive 3D Consensus Rendering.

Interfaces with ``py3Dmol`` to render protein-ligand complexes. Supports
high-clarity consensus visualization by rendering specific clusters
identified via HDBSCAN.

v2.0 — Refactored with robust error handling, validation, and fail-safe
visualization.
"""

from __future__ import annotations

import logging
import os
from typing import Dict, List, Optional, Tuple

import py3Dmol
from rdkit import Chem

from docking import DockingResult, EngineType

logger = logging.getLogger("poseai.visualizer")


# ─────────────────────────────────────────────────────────────────
# Configuration
# ─────────────────────────────────────────────────────────────────
# Mapping from file extensions to py3Dmol format strings
_FORMAT_MAP: Dict[str, str] = {
    ".pdbqt": "pdbqt",
    ".sdf": "sdf",
    ".mol2": "mol2",
    ".pdb": "pdb",
    ".dok": "pdb",
}

# Engine-specific colors for pose visualization
_DEFAULT_COLORS: Dict[EngineType, str] = {
    EngineType.SMINA: "#2ca02c",   # Green
    EngineType.GNINA: "#1f77b4",   # Blue
    EngineType.LEDOCK: "#d62728",  # Red
}

# Consensus pose color (distinctive)
CONSENSUS_COLOR = "#f1c40f"  # Gold

# Default view dimensions
DEFAULT_VIEW_WIDTH = 800
DEFAULT_VIEW_HEIGHT = 600


# ─────────────────────────────────────────────────────────────────
# Visualization Configuration
# ─────────────────────────────────────────────────────────────────
class VisualizerConfig:
    """Configuration for 3D visualization."""
    
    # View dimensions (pixels)
    DEFAULT_WIDTH = DEFAULT_VIEW_WIDTH
    DEFAULT_HEIGHT = DEFAULT_VIEW_HEIGHT
    
    # Molecule styling
    PROTEIN_STYLE = {"cartoon": {"color": "spectrum"}}
    CONSENSUS_STICK_STYLE = {"stick": {"color": CONSENSUS_COLOR, "radius": 0.2}}
    
    # Validation
    MAX_ATOMS_FOR_RENDERING = 50000  # Prevent browser crashes
    MIN_ATOMS_FOR_LIGAND = 5


# ─────────────────────────────────────────────────────────────────
# Main Visualizer
# ─────────────────────────────────────────────────────────────────
class DockingVisualizer:
    """Interactive 3D renderer for PoseAI consensus results.
    
    Parameters
    ----------
    receptor_path : str
        Path to receptor structure file (PDB, PDBQT, etc.).
    work_dir : str, optional
        Working directory for outputs.
        
    Raises
    ------
    FileNotFoundError
        If receptor file doesn't exist.
    ValueError
        If receptor format is unrecognized.
    """

    def __init__(
        self,
        receptor_path: str,
        work_dir: str = "/content/fast_lane",
    ) -> None:
        if not os.path.exists(receptor_path):
            raise FileNotFoundError(f"Receptor file not found: {receptor_path}")
        
        self.receptor_path = os.path.abspath(receptor_path)
        self.work_dir = work_dir
        self.view: Optional[py3Dmol.view] = None
        
        # Infer format from extension
        ext = os.path.splitext(receptor_path)[1].lower()
        self.rec_format = _FORMAT_MAP.get(ext, ext.lstrip("."))
        
        if self.rec_format is None:
            raise ValueError(f"Unrecognized receptor format: {ext}")
        
        logger.info(f"Initialized visualizer for {os.path.basename(receptor_path)} ({self.rec_format})")

    # ────────────────────────────────────────────────────────────
    # View Initialization
    # ────────────────────────────────────────────────────────────
    def create_interactive_view(
        self,
        width: int = VisualizerConfig.DEFAULT_WIDTH,
        height: int = VisualizerConfig.DEFAULT_HEIGHT,
    ) -> py3Dmol.view:
        """Initialise a py3Dmol view with the receptor loaded.
        
        Parameters
        ----------
        width : int, optional
            View width in pixels.
        height : int, optional
            View height in pixels.
            
        Returns
        -------
        py3Dmol.view
            The initialized view object.
            
        Raises
        ------
        IOError
            If receptor file cannot be read.
        ValueError
            If receptor contains too many atoms (overflow risk).
        """
        try:
            with open(self.receptor_path, "r") as fh:
                receptor_content = fh.read()
        except IOError as e:
            logger.error(f"Failed to read receptor: {e}")
            raise

        # Validate content size
        lines = receptor_content.split("\n")
        atom_lines = [l for l in lines if l.startswith(("ATOM", "HETATM"))]
        n_atoms = len(atom_lines)
        
        if n_atoms > VisualizerConfig.MAX_ATOMS_FOR_RENDERING:
            raise ValueError(
                f"Receptor has too many atoms ({n_atoms}); "
                f"exceeds browser rendering limit ({VisualizerConfig.MAX_ATOMS_FOR_RENDERING})"
            )
        
        logger.info(f"Rendering receptor with {n_atoms} atoms")
        
        try:
            self.view = py3Dmol.view(width=width, height=height)
            self.view.addModel(receptor_content, self.rec_format)
            self.view.setStyle(VisualizerConfig.PROTEIN_STYLE)
            self.view.zoomTo()
        except Exception as e:
            logger.error(f"Failed to create py3Dmol view: {e}", exc_info=True)
            raise
        
        return self.view

    # ────────────────────────────────────────────────────────────
    # Pose Rendering
    # ────────────────────────────────────────────────────────────
    def add_consensus_cluster(
        self,
        all_poses: List[Chem.Mol],
        cluster_indices: List[int],
        label: str = "Consensus",
    ) -> py3Dmol.view:
        """Render only the poses belonging to the identified consensus cluster.
        
        Parameters
        ----------
        all_poses : List[Chem.Mol]
            The full list of RDKit molecules used in clustering.
        cluster_indices : List[int]
            The indices of the 'winning' cluster from HDBSCAN.
        label : str, optional
            Descriptive label for the rendered cluster.
            
        Returns
        -------
        py3Dmol.view
            The updated view object.
            
        Raises
        ------
        ValueError
            If cluster_indices is empty or contains invalid indices.
        RuntimeError
            If view hasn't been initialized.
        """
        if self.view is None:
            logger.info("View not initialized; creating it now")
            self.create_interactive_view()

        if not cluster_indices:
            raise ValueError("cluster_indices cannot be empty")

        if not all_poses:
            raise ValueError("all_poses list is empty")

        # Validate indices
        invalid_indices = [i for i in cluster_indices if i >= len(all_poses) or i < 0]
        if invalid_indices:
            raise ValueError(
                f"Invalid cluster indices {invalid_indices} "
                f"(only {len(all_poses)} poses available)"
            )

        logger.info(f"Rendering {label} cluster with {len(cluster_indices)} poses")
        
        added_count = 0
        for idx in cluster_indices:
            mol = all_poses[idx]
            
            # Validate molecule
            if mol is None:
                logger.warning(f"Pose {idx} is None; skipping")
                continue
            
            if mol.GetNumAtoms() == 0:
                logger.warning(f"Pose {idx} has 0 atoms; skipping")
                continue

            try:
                # Convert RDKit Mol to SDF block
                block = Chem.MolToMolBlock(mol)
                
                if not block:
                    logger.warning(f"Failed to convert pose {idx} to MolBlock")
                    continue
                
                # Add to view
                self.view.addModel(block, "sdf")
                self.view.setStyle(
                    {"model": -1},
                    VisualizerConfig.CONSENSUS_STICK_STYLE,
                )
                added_count += 1
                
            except Exception as e:
                logger.error(f"Error rendering pose {idx}: {e}")
                continue

        if added_count == 0:
            logger.warning(f"No poses successfully added to {label} cluster")
        else:
            logger.info(f"Added {added_count}/{len(cluster_indices)} poses to visualization")
            self.view.zoomTo()

        return self.view

    def add_individual_poses(
        self,
        poses: List[Chem.Mol],
        engine_type: EngineType,
        label: Optional[str] = None,
    ) -> py3Dmol.view:
        """Add poses from a single engine with engine-specific coloring.
        
        Parameters
        ----------
        poses : List[Chem.Mol]
            RDKit molecules to render.
        engine_type : EngineType
            Engine that produced these poses (for color coding).
        label : str, optional
            Descriptive label.
            
        Returns
        -------
        py3Dmol.view
            Updated view.
        """
        if self.view is None:
            self.create_interactive_view()

        color = _DEFAULT_COLORS.get(engine_type, "#cccccc")
        label = label or engine_type.name
        
        logger.info(f"Adding {len(poses)} {label} poses (color: {color})")
        
        added = 0
        for idx, mol in enumerate(poses):
            if mol is None or mol.GetNumAtoms() == 0:
                logger.debug(f"Skipping invalid pose {idx} from {label}")
                continue
            
            try:
                block = Chem.MolToMolBlock(mol)
                if block:
                    self.view.addModel(block, "sdf")
                    self.view.setStyle(
                        {"model": -1},
                        {"stick": {"color": color, "radius": 0.15}},
                    )
                    added += 1
            except Exception as e:
                logger.debug(f"Error adding {label} pose {idx}: {e}")
        
        logger.info(f"Added {added}/{len(poses)} {label} poses")
        if added > 0:
            self.view.zoomTo()
        
        return self.view

    # ────────────────────────────────────────────────────────────
    # Display & Export
    # ────────────────────────────────────────────────────────────
    def show(self) -> None:
        """Render the WebGL visualization in the notebook.
        
        Raises
        ------
        RuntimeError
            If view is None (not initialized).
        """
        if self.view is None:
            raise RuntimeError(
                "View not initialized. Call create_interactive_view() first."
            )
        
        try:
            self.view.show()
        except Exception as e:
            logger.error(f"Failed to render view: {e}", exc_info=True)
            raise

    def get_html(self) -> str:
        """Export the current view as standalone HTML.
        
        Returns
        -------
        str
            HTML content that can be saved or embedded.
            
        Raises
        ------
        RuntimeError
            If view is None.
        """
        if self.view is None:
            raise RuntimeError("View not initialized")
        
        try:
            # Note: py3Dmol's HTML export is limited; this is a basic wrapper
            logger.info("Exporting view as HTML")
            return str(self.view)
        except Exception as e:
            logger.error(f"Failed to export HTML: {e}")
            raise

    # ────────────────────────────────────────────────────────────
    # Styling & Configuration
    # ────────────────────────────────────────────────────────────
    def set_protein_style(self, style: Dict) -> None:
        """Update protein backbone styling.
        
        Parameters
        ----------
        style : Dict
            py3Dmol style dictionary (e.g., {"cartoon": {"color": "blue"}}).
        """
        if self.view is None:
            logger.warning("View not initialized; style not applied")
            return
        
        try:
            self.view.setStyle(style)
            logger.info(f"Applied protein style: {style}")
        except Exception as e:
            logger.error(f"Failed to apply style: {e}")

    def reset_view(self) -> None:
        """Clear all poses and reset to receptor only."""
        if self.view is None:
            return
        
        try:
            self.view.clear()
            logger.info("View cleared")
            # Reload receptor
            with open(self.receptor_path, "r") as fh:
                receptor_content = fh.read()
            self.view.addModel(receptor_content, self.rec_format)
            self.view.setStyle(VisualizerConfig.PROTEIN_STYLE)
            self.view.zoomTo()
        except Exception as e:
            logger.error(f"Failed to reset view: {e}")
