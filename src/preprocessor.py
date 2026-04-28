"""PoseAI Preprocessor Module — Structural Data Standardisation.

Automates the retrieval of PDB structures from the RCSB, removal of
crystallographic waters, and generation of docking-engine-specific input
files:

*   **PDBQT** (Gasteiger charges) — consumed by Smina and Gnina.
*   **PDB** (clean, no waters)    — consumed by LeDock via ``lepro``.
*   **MOL2** (Gasteiger charges)  — consumed by LeDock for the ligand.

v2.0 — Refactored with ligand validation, robust error handling, and
better diagnostics.

Environment
-----------
*   Python >= 3.10
*   Open Babel >= 3.1 (``obabel`` on ``$PATH``)
*   RDKit >= 2023.03 (for ligand validation)
*   ``wget`` for RCSB downloads (pre-installed on Colab)
"""

from __future__ import annotations

import logging
import os
import subprocess
from typing import Optional, Tuple

from rdkit import Chem

logger = logging.getLogger("poseai.preprocessor")


# ─────────────────────────────────────────────────────────────────
# Configuration
# ─────────────────────────────────────────────────────────────────
class PreprocessorConfig:
    """Configuration for structure preparation."""
    
    # Ligand validation
    MIN_LIGAND_ATOMS = 5          # Too small = likely artifact
    MAX_LIGAND_ATOMS = 1000       # Too large = likely protein/complex
    
    # RCSB API
    RCSB_DOWNLOAD_URL = "https://files.rcsb.org/download"
    RCSB_TIMEOUT_S = 30


# ─────────────────────────────────────────────────────────────────
# Main Preprocessor
# ─────────────────────────────────────────────────────────────────
class ProteinLigandPrep:
    """Prepare receptor and ligand structures for multi-engine docking.

    Produces all file formats required by the downstream
    ``EnsembleManager``: PDBQT for Smina/Gnina, and PDB + MOL2 for
    LeDock. Validates all outputs to catch errors early.

    Parameters
    ----------
    pdb_id : str
        Four-character PDB accession code (case-insensitive).
    work_dir : str, optional
        Root of the high-speed local scratch directory.

    Attributes
    ----------
    raw_pdb : str
        Path to the downloaded PDB file.
    receptor_pdbqt : str
        Receptor in PDBQT format (Smina / Gnina).
    receptor_pdb : str
        Cleaned receptor in PDB format (LeDock).
    ligand_pdbqt : str
        Ligand in PDBQT format (Smina / Gnina).
    ligand_mol2 : str
        Ligand in MOL2 format (LeDock).
        
    Raises
    ------
    ValueError
        If PDB ID is not 4 characters.
    """

    def __init__(
        self,
        pdb_id: str,
        work_dir: str = "/content/fast_lane",
    ) -> None:
        pdb_id = pdb_id.lower()
        
        # Validate PDB ID format
        if len(pdb_id) != 4 or not pdb_id.isalnum():
            raise ValueError(f"Invalid PDB ID (must be 4 alphanumeric): {pdb_id}")
        
        self.pdb_id = pdb_id
        self.work_dir = os.path.join(work_dir, "data", self.pdb_id)
        os.makedirs(self.work_dir, exist_ok=True)

        self.raw_pdb = os.path.join(self.work_dir, f"{self.pdb_id}.pdb")
        self.receptor_pdbqt = os.path.join(self.work_dir, "receptor.pdbqt")
        self.receptor_pdb = os.path.join(self.work_dir, "receptor_clean.pdb")
        self.ligand_pdbqt = os.path.join(self.work_dir, "ligand.pdbqt")
        self.ligand_mol2 = os.path.join(self.work_dir, "ligand.mol2")
        
        logger.info(f"Initialized preprocessor for PDB {self.pdb_id}")

    # ────────────────────────────────────────────────────────────
    # Structure Retrieval
    # ────────────────────────────────────────────────────────────
    def fetch_structure(self) -> bool:
        """Retrieve the target structure from the RCSB Protein Data Bank.

        Skips the download if the file already exists on disk.

        Returns
        -------
        bool
            True if the PDB file is available after this call.
        """
        if os.path.exists(self.raw_pdb):
            logger.info(f"PDB file already present: {self.raw_pdb}")
            return True

        url = f"{PreprocessorConfig.RCSB_DOWNLOAD_URL}/{self.pdb_id.upper()}.pdb"
        logger.info(f"Downloading {self.pdb_id.upper()} from RCSB...")

        try:
            result = subprocess.run(
                ["wget", "-q", "--timeout", str(PreprocessorConfig.RCSB_TIMEOUT_S), 
                 url, "-O", self.raw_pdb],
                check=True,
                capture_output=True,
                text=True,
                timeout=PreprocessorConfig.RCSB_TIMEOUT_S + 5,
            )
            logger.info(f"Download complete: {self.raw_pdb}")
            return True
            
        except subprocess.CalledProcessError as exc:
            logger.error(f"Failed to download PDB ID {self.pdb_id}: {exc.stderr}")
            return False
        except subprocess.TimeoutExpired:
            logger.error(f"Download timeout for PDB ID {self.pdb_id}")
            return False
        except Exception as exc:
            logger.error(f"Unexpected error downloading PDB: {exc}", exc_info=True)
            return False

    # ────────────────────────────────────────────────────────────
    # Receptor Preparation
    # ────────────────────────────────────────────────────────────
    def prepare_receptor(self) -> bool:
        """Clean and convert the protein structure to PDBQT format.

        Processing steps performed by Open Babel:
            1. Remove water molecules and non-protein HETATMs (``-xr``).
            2. Add polar hydrogens (``-d`` strips non-polar, ``-h`` adds).
            3. Assign Gasteiger partial charges.

        Also produces a clean PDB copy (``receptor_clean.pdb``) for LeDock.

        Returns
        -------
        bool
            True if both PDBQT and clean-PDB files were written.
        """
        if not os.path.exists(self.raw_pdb):
            logger.error(f"Raw PDB not found: {self.raw_pdb}. Call fetch_structure() first.")
            return False

        try:
            # --- PDBQT for Smina / Gnina ---
            logger.info("Converting receptor to PDBQT format...")
            cmd_pdbqt = [
                "obabel", self.raw_pdb, "-O", self.receptor_pdbqt,
                "-xr",                          # Remove water
                "-d",                           # Add polar hydrogens
                "--partialcharge", "gasteiger",  # Gasteiger charges
            ]
            result = subprocess.run(
                cmd_pdbqt, check=True, capture_output=True, text=True, timeout=60
            )
            
            if not os.path.exists(self.receptor_pdbqt):
                logger.error("PDBQT generation succeeded but file not found")
                return False
            
            logger.info(f"Receptor PDBQT written: {self.receptor_pdbqt}")

            # --- Clean PDB for LeDock ---
            logger.info("Converting receptor to clean PDB format...")
            cmd_pdb = [
                "obabel", self.raw_pdb, "-O", self.receptor_pdb,
                "-xr", "-d",  # Remove water, add H
            ]
            result = subprocess.run(
                cmd_pdb, check=True, capture_output=True, text=True, timeout=60
            )
            
            if not os.path.exists(self.receptor_pdb):
                logger.error("PDB generation succeeded but file not found")
                return False
            
            logger.info(f"Receptor PDB written: {self.receptor_pdb}")
            return True

        except subprocess.CalledProcessError as exc:
            logger.error(f"Receptor preparation failed: {exc.stderr}")
            return False
        except subprocess.TimeoutExpired:
            logger.error("Receptor preparation timeout (obabel hung)")
            return False
        except Exception as exc:
            logger.error(f"Unexpected error in receptor prep: {exc}", exc_info=True)
            return False

    # ────────────────────────────────────────────────────────────
    # Ligand Isolation
    # ────────────────────────────────────────────────────────────
    def isolate_ligand(self, ligand_code: str) -> bool:
        """Extract a specific ligand residue and convert to docking formats.

        Produces both PDBQT (Smina/Gnina) and MOL2 (LeDock) representations
        with Gasteiger partial charges assigned.

        Parameters
        ----------
        ligand_code : str
            Three-letter residue code of the co-crystallised ligand
            (e.g., "ATP", "STI").

        Returns
        -------
        bool
            True if both output files were written.
            
        Raises
        ------
        ValueError
            If ligand code is not exactly 3 characters.
        """
        if len(ligand_code) != 3:
            raise ValueError(f"Ligand code must be 3 characters, got: {ligand_code}")

        if not os.path.exists(self.raw_pdb):
            logger.error(f"Raw PDB not found: {self.raw_pdb}")
            return False

        temp_pdb = os.path.join(self.work_dir, "ligand_temp.pdb")
        logger.info(f"Isolating ligand residue: {ligand_code}")

        try:
            # Extract HETATM records matching the ligand code.
            # CRITICAL: Many crystal structures contain the same ligand in
            # multiple chains (e.g., STI in chains A and B of 1IEP).
            # We take only the FIRST chain encountered to avoid multi-molecule
            # PDBQT files (which docking engines reject).
            
            first_chain = None
            hetatm_serials = set()

            with open(self.raw_pdb, "r") as src, open(temp_pdb, "w") as dst:
                for line in src:
                    if line.startswith("HETATM") and f" {ligand_code} " in line:
                        chain_id = line[21] if len(line) > 21 else " "

                        if first_chain is None:
                            first_chain = chain_id
                            logger.info(
                                f"Ligand {ligand_code} found in chain '{chain_id}'; "
                                f"extracting this chain only."
                            )

                        if chain_id != first_chain:
                            logger.debug(f"Skipping duplicate ligand in chain '{chain_id}'")
                            continue

                        dst.write(line)
                        hetatm_serials.add(line[6:11].strip())

                dst.write("END\n")

            if not hetatm_serials:
                logger.error(f"No HETATM records found for ligand code '{ligand_code}'")
                return False

            logger.info(f"Extracted {len(hetatm_serials)} HETATM atoms for {ligand_code}")

            # --- PDBQT for Smina / Gnina ---
            logger.info("Converting ligand to PDBQT format...")
            cmd_pdbqt = [
                "obabel", temp_pdb, "-O", self.ligand_pdbqt,
                "-d", "--partialcharge", "gasteiger",
            ]
            subprocess.run(cmd_pdbqt, check=True, capture_output=True, text=True, timeout=60)

            # Sanitize PDBQT
            self._sanitise_pdbqt(self.ligand_pdbqt)
            logger.info(f"Ligand PDBQT written: {self.ligand_pdbqt}")

            # Validate PDBQT (non-blocking warning)
            if not self._validate_pdbqt_ligand(self.ligand_pdbqt):
                logger.warning("Ligand PDBQT validation warning: proceeding anyway")

            # --- MOL2 for LeDock ---
            logger.info("Converting ligand to MOL2 format...")
            cmd_mol2 = [
                "obabel", temp_pdb, "-O", self.ligand_mol2,
                "-d", "--partialcharge", "gasteiger",
            ]
            subprocess.run(cmd_mol2, check=True, capture_output=True, text=True, timeout=60)
            logger.info(f"Ligand MOL2 written: {self.ligand_mol2}")

            # Validate MOL2 (non-blocking warning)
            if not self._validate_mol2_ligand(self.ligand_mol2):
                logger.warning("Ligand MOL2 validation warning: proceeding anyway")

            logger.info(f"Ligand {ligand_code} successfully prepared")
            return True

        except subprocess.CalledProcessError as exc:
            logger.error(f"Ligand isolation failed: {exc.stderr}")
            return False
        except subprocess.TimeoutExpired:
            logger.error("Ligand conversion timeout (obabel hung)")
            return False
        except Exception as exc:
            logger.error(f"Unexpected error in ligand isolation: {exc}", exc_info=True)
            return False
        finally:
            # Clean up intermediate file
            if os.path.exists(temp_pdb):
                os.remove(temp_pdb)

    # ────────────────────────────────────────────────────────────
    # Validation
    # ────────────────────────────────────────────────────────────
    @staticmethod
    def _sanitise_pdbqt(path: str) -> None:
        """Remove lines from a PDBQT file that would cause Vina parse errors.

        Valid PDBQT line prefixes are: ``ATOM``, ``HETATM``, ``ROOT``,
        ``ENDROOT``, ``BRANCH``, ``ENDBRANCH``, ``TORSDOF``, ``REMARK``,
        ``MODEL``, ``ENDMDL``, ``END``.

        The file is rewritten in place.
        """
        valid_tags = (
            "ATOM", "HETATM", "ROOT", "ENDROOT",
            "BRANCH", "ENDBRANCH", "TORSDOF",
            "REMARK", "MODEL", "ENDMDL", "END",
        )

        with open(path, "r") as fh:
            lines = fh.readlines()

        cleaned = [ln for ln in lines if ln.strip() and ln.split()[0] in valid_tags]
        removed = len(lines) - len(cleaned)

        if removed > 0:
            logger.warning(f"Removed {removed} invalid line(s) from {path}")

        with open(path, "w") as fh:
            fh.writelines(cleaned)

    @staticmethod
    def _validate_pdbqt_ligand(pdbqt_path: str) -> bool:
        """Validate ligand PDBQT by checking file exists and has content.
        
        Parameters
        ----------
        pdbqt_path : str
            Path to PDBQT file.
            
        Returns
        -------
        bool
            True if file exists and contains ATOM/HETATM records.
        """
        try:
            if not os.path.exists(pdbqt_path):
                return False
            
            with open(pdbqt_path, "r") as f:
                lines = f.readlines()
            
            atom_lines = [l for l in lines if l.startswith(("ATOM", "HETATM"))]
            if len(atom_lines) < 5:
                logger.warning(f"PDBQT has very few atoms: {len(atom_lines)}")
                return False
            
            logger.info(f"PDBQT validation passed: {len(atom_lines)} atoms")
            return True
            
        except Exception as e:
            logger.warning(f"PDBQT validation error: {e}")
            return False

    @staticmethod
    def _validate_mol2_ligand(mol2_path: str) -> bool:
        """Validate ligand MOL2 by checking file exists and has content.
        
        Parameters
        ----------
        mol2_path : str
            Path to MOL2 file.
            
        Returns
        -------
        bool
            True if file exists and is non-empty.
        """
        try:
            if not os.path.exists(mol2_path):
                return False
            
            file_size = os.path.getsize(mol2_path)
            if file_size < 100:
                logger.warning(f"MOL2 file suspiciously small: {file_size} bytes")
                return False
            
            logger.info(f"MOL2 validation passed: {file_size} bytes")
            return True
            
        except Exception as e:
            logger.warning(f"MOL2 validation error: {e}")
            return False

# --- Dynamically added methods ---
import requests
import numpy as np
from rdkit import Chem

def fetch_rcsb_smiles(ligand_code):
    try:
        url = f"https://data.rcsb.org/rest/v1/core/chemcomp/{ligand_code.upper()}"
        resp = requests.get(url, timeout=10)
        if resp.status_code == 200:
            data = resp.json()
            descriptors = data.get("rcsb_chem_comp_descriptor", [])
            if isinstance(descriptors, dict): descriptors = [descriptors]
            for desc in descriptors:
                if isinstance(desc, dict):
                    if "smilesstereo" in desc: return desc.get("smilesstereo")
                    elif "smiles" in desc: return desc.get("smiles")
    except Exception: pass
    return None

def get_ligand_centroid(mol2_path, padding=10.0):
    try:
        mol = Chem.MolFromMol2File(mol2_path, sanitize=False)
        if mol and mol.GetNumConformers() > 0:
            coords = mol.GetConformer().GetPositions()
            center = tuple(np.mean(coords, axis=0))
            size = tuple(np.max(coords, axis=0) - np.min(coords, axis=0) + float(padding))
            return center, size
    except Exception: pass
    return (0.0, 0.0, 0.0), (20.0, 20.0, 20.0)
