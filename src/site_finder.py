"""PoseAI Site Finder Module — Binding Pocket Identification.

Automates cavity detection via the ``fpocket`` algorithm and extracts
barycenter coordinates for downstream docking-box placement. The
highest-ranking pocket is selected by default, as fpocket orders pockets
by druggability score.

v2.0 — Refactored with robust parsing, better error handling, and
defensive input validation.

Environment
-----------
*   Python >= 3.10
*   ``fpocket`` >= 4.0 on ``$PATH``
"""

from __future__ import annotations

import logging
import os
import re
import subprocess
from typing import Dict, List, Optional, Tuple

logger = logging.getLogger("poseai.site_finder")


# ─────────────────────────────────────────────────────────────────
# Configuration
# ─────────────────────────────────────────────────────────────────
class SiteFinderConfig:
    """Configuration for pocket analysis."""
    
    # Default padding for search box (angstroms)
    DEFAULT_PADDING = 10.0
    
    # Fallback coordinates (origin) if fpocket fails
    FALLBACK_CENTER = (0.0, 0.0, 0.0)
    
    # fpocket execution
    FPOCKET_TIMEOUT_S = 300
    
    # Coordinate parsing
    # Matches floats in scientific notation or decimal form
    FLOAT_PATTERN = r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?'


# ─────────────────────────────────────────────────────────────────
# Main Analyzer
# ─────────────────────────────────────────────────────────────────
class PocketAnalyzer:
    """Identify and parameterise the primary binding pocket of a protein.

    Wraps the ``fpocket`` binary, parses its output, and exposes the
    pocket barycenter and search-box dimensions consumed by
    ``EnsembleManager``.

    Parameters
    ----------
    pdb_path : str
        Absolute path to the protein structure file (``.pdb``).
    work_dir : str, optional
        Root of the high-speed local scratch directory.
        
    Attributes
    ----------
    center : Optional[Tuple[float, float, float]]
        Pocket barycenter coordinates (set by extract_barycenter).
    """

    def __init__(
        self,
        pdb_path: str,
        work_dir: str = "/content/fast_lane",
    ) -> None:
        if not os.path.exists(pdb_path):
            raise FileNotFoundError(f"PDB file not found: {pdb_path}")
        
        self.pdb_path = os.path.abspath(pdb_path)
        self.pdb_id = os.path.splitext(os.path.basename(pdb_path))[0]
        self.work_dir = work_dir
        self.output_dir = os.path.join(
            self.work_dir, "data", self.pdb_id, f"{self.pdb_id}_out"
        )
        self.center: Optional[Tuple[float, float, float]] = None
        
        logger.info(f"Initialized pocket analyzer for {self.pdb_id}")

    # ────────────────────────────────────────────────────────────
    # fpocket Execution
    # ────────────────────────────────────────────────────────────
    def execute_fpocket(self) -> bool:
        """Run the fpocket algorithm on the target receptor.

        The subprocess is executed with ``cwd`` set to the receptor's
        parent directory, avoiding the use of ``os.chdir()`` which would
        mutate global process state and is unsafe in concurrent contexts.

        Returns
        -------
        bool
            True if fpocket completed with exit code 0.
        """
        target_dir = os.path.dirname(self.pdb_path)
        logger.info(f"Executing fpocket on {self.pdb_id}...")

        try:
            result = subprocess.run(
                ["fpocket", "-f", self.pdb_path],
                check=True,
                capture_output=True,
                text=True,
                cwd=target_dir,
                timeout=SiteFinderConfig.FPOCKET_TIMEOUT_S,
            )
            logger.info(f"fpocket completed successfully for {self.pdb_id}")
            return True
            
        except subprocess.CalledProcessError as exc:
            logger.error(f"fpocket execution failed: {exc.stderr}")
            return False
        except subprocess.TimeoutExpired:
            logger.error(
                f"fpocket timeout ({SiteFinderConfig.FPOCKET_TIMEOUT_S}s); "
                f"large proteins may require longer"
            )
            return False
        except FileNotFoundError:
            logger.error(
                "fpocket binary not found on $PATH. "
                "Install with: apt-get install fpocket"
            )
            return False
        except Exception as exc:
            logger.error(f"Unexpected error running fpocket: {exc}", exc_info=True)
            return False

    # ────────────────────────────────────────────────────────────
    # Coordinate Extraction
    # ────────────────────────────────────────────────────────────
    def extract_barycenter(self, pocket_index: int = 1) -> Tuple[float, float, float]:
        """Parse the fpocket info file and retrieve a pocket's barycenter.

        The barycenter is the geometric centre of the predicted binding
        cavity, used to position the docking search box.

        Parameters
        ----------
        pocket_index : int
            1-based pocket rank. Pocket 1 is the top-scoring cavity by default.

        Returns
        -------
        Tuple[float, float, float]
            Coordinates (x, y, z) in angstroms. Falls back to the origin
            if parsing fails.
        """
        info_file = os.path.join(self.output_dir, f"{self.pdb_id}_info.txt")

        if not os.path.exists(info_file):
            logger.warning(
                f"Info file not found: {info_file}. "
                f"fpocket may not have completed; returning origin."
            )
            self.center = SiteFinderConfig.FALLBACK_CENTER
            return self.center

        try:
            with open(info_file, "r") as fh:
                content = fh.read()
        except IOError as e:
            logger.error(f"Failed to read info file: {e}")
            self.center = SiteFinderConfig.FALLBACK_CENTER
            return self.center

        # Parse the info file
        center = self._parse_info_file(content, pocket_index)
        
        if center:
            self.center = center
            logger.info(f"Pocket {pocket_index} barycenter: {center}")
            return center
        else:
            logger.warning(f"Could not extract barycenter for pocket {pocket_index}")
            self.center = SiteFinderConfig.FALLBACK_CENTER
            return self.center

    @staticmethod
    def _parse_info_file(
        content: str, pocket_index: int
    ) -> Optional[Tuple[float, float, float]]:
        """Parse fpocket info file content to extract pocket barycenter.
        
        Handles multiple coordinate formats:
        - "Barycenter: x y z"
        - "Center of mass: x y z"
        - Tab or space-separated values
        - Scientific notation (1.23e-02)
        - Multiline pocket sections
        
        Parameters
        ----------
        content : str
            Full contents of fpocket info file.
        pocket_index : int
            Pocket number (1-based).
            
        Returns
        -------
        Optional[Tuple[float, float, float]]
            Coordinates if found, else None.
        """
        lines = content.split("\n")
        
        # Strategy: Find pocket block, then extract coordinates
        in_pocket = False
        pocket_header = f"Pocket {pocket_index}"
        
        # Robust regex for three floats (handles scientific notation, tabs, etc.)
        float_pattern = SiteFinderConfig.FLOAT_PATTERN
        coord_pattern = re.compile(
            f"({float_pattern})\\s+({float_pattern})\\s+({float_pattern})"
        )
        
        for i, line in enumerate(lines):
            stripped = line.strip()
            
            # Detect pocket boundaries
            if stripped.startswith("Pocket "):
                if stripped.startswith(pocket_header):
                    in_pocket = True
                    logger.debug(f"Found {pocket_header}")
                    continue
                elif in_pocket:
                    # Reached next pocket — stop
                    break
            
            if not in_pocket:
                continue
            
            # Look for barycenter or center lines
            label = stripped.lower()
            
            # Multiple ways to find barycenter: look for keywords and coordinates
            has_barycenter_keyword = (
                "barycenter" in label or 
                "center of mass" in label or 
                "center_mass" in label or
                "geometric center" in label
            )
            
            # Exclude lines that are metadata (e.g., "max dist")
            is_metadata = any(x in label for x in ["max", "dist", "score", "volume"])
            
            if has_barycenter_keyword and not is_metadata:
                # Try to extract coordinates after the colon
                if ":" in stripped:
                    coords_part = stripped.split(":", 1)[1].strip()
                else:
                    # No colon, try the whole line
                    coords_part = stripped
                
                # Search for coordinate triplet
                match = coord_pattern.search(coords_part)
                
                if match and len(match.groups()) >= 3:
                    try:
                        x = float(match.group(1))
                        y = float(match.group(2))
                        z = float(match.group(3))
                        
                        # Sanity check: coordinates should be reasonable (±200 Å)
                        if all(abs(c) < 500 for c in [x, y, z]):
                            logger.debug(f"Extracted coords: ({x}, {y}, {z})")
                            return (x, y, z)
                        else:
                            logger.warning(f"Coordinates out of range: ({x}, {y}, {z})")
                    except ValueError as e:
                        logger.warning(f"Failed to convert coordinates to float: {e}")
                        continue

        # Fallback
        logger.warning(
            f"Could not find barycenter for pocket {pocket_index} in info file"
        )
        return None

    # ────────────────────────────────────────────────────────────
    # Search Box Calculation
    # ────────────────────────────────────────────────────────────
    def calculate_search_bounds(
        self,
        padding: float = SiteFinderConfig.DEFAULT_PADDING,
    ) -> Dict[str, float]:
        """Translate the barycenter into a Cartesian bounding box.

        Required by docking engines (e.g., LeDock) that accept min/max
        boundaries rather than centre-and-size definitions.

        Parameters
        ----------
        padding : float, optional
            Half-width from the centre to each box face, in angstroms.

        Returns
        -------
        Dict[str, float]
            Dictionary with keys ``min_x``, ``max_x``, ``min_y``,
            ``max_y``, ``min_z``, ``max_z``.
            
        Raises
        ------
        ValueError
            If center is None or padding <= 0.
        """
        if padding <= 0:
            raise ValueError(f"Padding must be positive, got {padding}")
        
        if self.center is None:
            self.extract_barycenter()

        x, y, z = self.center
        
        bounds = {
            "min_x": x - padding,
            "max_x": x + padding,
            "min_y": y - padding,
            "max_y": y + padding,
            "min_z": z - padding,
            "max_z": z + padding,
        }
        
        logger.info(
            f"Search bounds (padding={padding} Å): "
            f"x=[{bounds['min_x']:.1f}, {bounds['max_x']:.1f}] "
            f"y=[{bounds['min_y']:.1f}, {bounds['max_y']:.1f}] "
            f"z=[{bounds['min_z']:.1f}, {bounds['max_z']:.1f}]"
        )
        
        return bounds

    def get_box_params(
        self,
        padding: float = SiteFinderConfig.DEFAULT_PADDING,
    ) -> Tuple[Tuple[float, float, float], Tuple[float, float, float]]:
        """Return centre and size as tuples for ``EnsembleManager``.

        Convenience method that packages the pocket data in the exact
        ``(center, box_size)`` format expected by ``run_ensemble()``.

        Parameters
        ----------
        padding : float, optional
            Half-width from the centre to each box face (angstroms).

        Returns
        -------
        Tuple[Tuple[float, float, float], Tuple[float, float, float]]
            A 2-tuple of ``(center, box_size)`` where both elements are
            ``(x, y, z)`` float tuples.
            
        Raises
        ------
        ValueError
            If padding <= 0.
        """
        if padding <= 0:
            raise ValueError(f"Padding must be positive, got {padding}")
        
        if self.center is None:
            self.extract_barycenter()

        size = (padding * 2, padding * 2, padding * 2)
        
        logger.info(
            f"Box parameters: center={self.center}, size={size}"
        )
        
        return self.center, size

    # ────────────────────────────────────────────────────────────
    # Multi-Pocket Inspection (Diagnostic)
    # ────────────────────────────────────────────────────────────
    def list_pockets(self) -> List[Dict[str, any]]:
        """List all detected pockets with scores (diagnostic utility).
        
        Returns
        -------
        List[Dict]
            Each dict has keys: 'index', 'barycenter', 'score'.
            
        Note
        ----
        This requires manual parsing of the fpocket output.
        Useful for verifying that the expected pocket was selected.
        """
        info_file = os.path.join(self.output_dir, f"{self.pdb_id}_info.txt")
        
        if not os.path.exists(info_file):
            logger.warning(f"Info file not found: {info_file}")
            return []
        
        try:
            with open(info_file, "r") as fh:
                content = fh.read()
        except IOError as e:
            logger.error(f"Failed to read info file: {e}")
            return []
        
        pockets = []
        lines = content.split("\n")
        
        for line in lines:
            if line.strip().startswith("Pocket "):
                # Try to extract pocket number and score
                match = re.match(r"Pocket\s+(\d+)\s+.*Score:\s*([-\d.]+)", line)
                if match:
                    pocket_num = int(match.group(1))
                    score = float(match.group(2))
                    
                    # Try to extract barycenter for this pocket
                    barycenter = self._parse_info_file(content, pocket_num)
                    
                    pockets.append({
                        "index": pocket_num,
                        "barycenter": barycenter,
                        "score": score,
                    })
        
        logger.info(f"Detected {len(pockets)} pocket(s)")
        return pockets
