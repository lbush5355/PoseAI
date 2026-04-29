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
import re
import subprocess
from typing import Optional, Tuple

import numpy as np
import requests
from rdkit import Chem

from config import get_config

logger = logging.getLogger("poseai.preprocessor")


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

        cfg = get_config().preprocessor
        url = f"{cfg.rcsb_download_url}/{self.pdb_id.upper()}.pdb"
        logger.info(f"Downloading {self.pdb_id.upper()} from RCSB...")

        try:
            result = subprocess.run(
                ["wget", "-q", "--timeout", str(int(cfg.rcsb_timeout_s)),
                 url, "-O", self.raw_pdb],
                check=True,
                capture_output=True,
                text=True,
                timeout=cfg.rcsb_timeout_s + 5,
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
            obabel_timeout = get_config().preprocessor.obabel_timeout_s
            cmd_pdbqt = [
                "obabel", self.raw_pdb, "-O", self.receptor_pdbqt,
                "-xr",                          # Remove water
                "-d",                           # Add polar hydrogens
                "--partialcharge", "gasteiger",  # Gasteiger charges
            ]
            result = subprocess.run(
                cmd_pdbqt, check=True, capture_output=True, text=True,
                timeout=obabel_timeout,
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
                cmd_pdb, check=True, capture_output=True, text=True,
                timeout=obabel_timeout,
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
            RCSB residue code of the co-crystallised ligand. Per the
            PDB chemical component dictionary, codes are 1-3
            alphanumeric characters (e.g., "STI", "N3", "ZN").

        Returns
        -------
        bool
            True if both output files were written.

        Raises
        ------
        ValueError
            If ligand code is not 1-3 alphanumeric characters.
        """
        if not (1 <= len(ligand_code) <= 3) or not ligand_code.isalnum():
            raise ValueError(
                f"Ligand code must be 1-3 alphanumeric characters, got: {ligand_code!r}"
            )

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
            obabel_timeout = get_config().preprocessor.obabel_timeout_s
            logger.info("Converting ligand to PDBQT format...")
            cmd_pdbqt = [
                "obabel", temp_pdb, "-O", self.ligand_pdbqt,
                "-d", "--partialcharge", "gasteiger",
            ]
            subprocess.run(cmd_pdbqt, check=True, capture_output=True, text=True,
                           timeout=obabel_timeout)

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
            subprocess.run(cmd_mol2, check=True, capture_output=True, text=True,
                           timeout=obabel_timeout)
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
            
            min_atoms = get_config().preprocessor.min_ligand_atoms
            atom_lines = [l for l in lines if l.startswith(("ATOM", "HETATM"))]
            if len(atom_lines) < min_atoms:
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

def fetch_rcsb_smiles(ligand_code: str) -> Optional[str]:
    """Fetch canonical SMILES for a ligand from the RCSB chemcomp API.

    Returns the stereochemistry-aware SMILES when available, otherwise
    the plain SMILES. Returns None if the ligand is not found, the
    response is malformed, or no SMILES descriptor is present. Network
    and parse errors are logged at WARNING and converted to None so
    callers can fall back to a user-supplied SMILES.
    """
    url = f"https://data.rcsb.org/rest/v1/core/chemcomp/{ligand_code.upper()}"

    try:
        resp = requests.get(url, timeout=10)
        resp.raise_for_status()
        data = resp.json()
    except requests.RequestException as e:
        logger.warning(f"RCSB fetch failed for ligand {ligand_code}: {e}")
        return None
    except ValueError as e:
        logger.warning(f"Invalid JSON from RCSB for ligand {ligand_code}: {e}")
        return None

    descriptor = data.get("rcsb_chem_comp_descriptor")
    if not isinstance(descriptor, dict):
        logger.warning(
            f"Unexpected RCSB descriptor shape for ligand {ligand_code}: "
            f"{type(descriptor).__name__}"
        )
        return None

    # Prefer stereo-aware SMILES; fall back to plain SMILES
    smiles = descriptor.get("SMILES_stereo") or descriptor.get("SMILES")
    if smiles:
        return smiles

    logger.warning(f"No SMILES descriptor in RCSB response for ligand {ligand_code}")
    return None


# Common buffers, ions, crystallization additives, and standard amino acids
# — filtered out when picking a drug-like ligand from an RCSB entry's
# non-polymer entities. Amino acids appear as non-polymer entries in entries
# whose actual inhibitor is a polymer chain (e.g. 1A30); without filtering
# they get picked up as bogus 10-atom templates against the real ligand.
_NON_DRUG_LIKE_CODES = {
    "HOH", "DOD", "NA", "K", "MG", "CA", "ZN", "FE", "MN", "NI", "CU", "CO",
    "CD", "HG", "PB", "CL", "BR", "F", "I", "SO4", "PO4", "NO3", "CO3",
    "ACT", "EDO", "GOL", "PEG", "PG4", "PG6", "MES", "TRS", "BME", "DMS",
    "DMF", "MOH", "EOH", "IMD", "FMT", "ACE", "CIT", "MPD", "EPE", "BCT",
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
}


def fetch_rcsb_entry_ligand_codes(pdb_id: str) -> List[str]:
    """Look up non-polymer ligand chemical component IDs for a PDB entry.

    Uses RCSB's GraphQL API to get all non-polymer entities and their
    chemical component codes in a single request. Common buffers, ions,
    and crystallization additives are filtered out so the drug-like
    ligand candidates rise to the top.

    Returns a list of candidate codes (often just one for a typical
    drug-bound structure). Returns an empty list if the entry has no
    non-polymer entities (e.g., peptidomimetic inhibitors that are
    deposited as polymer chains rather than HETATMs — 1A30 is one),
    or if the API call fails.
    """
    query = (
        "query($id: String!) { "
        "entry(entry_id: $id) { "
        "nonpolymer_entities { pdbx_entity_nonpoly { comp_id } } "
        "} }"
    )

    try:
        resp = requests.post(
            "https://data.rcsb.org/graphql",
            json={"query": query, "variables": {"id": pdb_id.upper()}},
            timeout=10,
        )
        resp.raise_for_status()
        data = resp.json()
    except requests.RequestException as e:
        logger.warning(f"RCSB entry lookup failed for {pdb_id}: {e}")
        return []
    except ValueError as e:
        logger.warning(f"RCSB entry lookup returned invalid JSON for {pdb_id}: {e}")
        return []

    entry = (data.get("data") or {}).get("entry") or {}
    entities = entry.get("nonpolymer_entities") or []

    codes: List[str] = []
    for entity in entities:
        comp_id = ((entity or {}).get("pdbx_entity_nonpoly") or {}).get("comp_id")
        if not comp_id:
            continue
        upper = comp_id.upper()
        if upper in _NON_DRUG_LIKE_CODES:
            continue
        if upper not in codes:
            codes.append(upper)

    return codes


def fetch_rcsb_ideal_sdf(
    ligand_code: str, save_dir: str
) -> Optional[str]:
    """Download the RCSB chemical component dictionary's ideal SDF for a ligand.

    The ideal SDF provides a canonical, RDKit-friendly 3D structure with
    bond orders curated by RCSB. Used as a topology template when a
    deposited mol2's bond perception conflicts with engine outputs
    (common for charged peptidomimetic ligands).

    The file is cached on disk after first download. Returns None on
    network failure, missing-ligand 404, or filesystem error so the
    caller can fall back to other topology sources or other candidate
    codes.
    """
    code = ligand_code.upper()
    save_path = os.path.join(save_dir, f"{code}_ideal.sdf")
    if os.path.exists(save_path):
        return save_path

    url = f"https://files.rcsb.org/ligands/download/{code}_ideal.sdf"
    try:
        resp = requests.get(url, timeout=10)
        resp.raise_for_status()
    except requests.RequestException as e:
        logger.warning(f"RCSB ideal SDF fetch failed for {code}: {e}")
        return None

    if not resp.text.strip():
        logger.warning(f"RCSB returned empty ideal SDF for {code}")
        return None

    try:
        os.makedirs(save_dir, exist_ok=True)
        with open(save_path, "w") as f:
            f.write(resp.text)
    except OSError as e:
        logger.error(f"Could not write ideal SDF to {save_path}: {e}")
        return None

    logger.warning(f"Cached ideal SDF for {code} → {save_path}")
    return save_path


def extract_ligand_code_candidates_from_mol2(mol2_path: str) -> List[str]:
    """Generate candidate RCSB ligand codes from a MOL2 file's first ATOM record.

    Mol2 substructure names typically concatenate the residue code with
    a residue number (e.g., "MK11", "STI301", "Z34"), and the boundary
    isn't unambiguous without external context — codes themselves can
    end in digits. We progressively strip trailing digits and yield each
    1-3 character candidate; the caller validates by attempting an
    actual fetch (first 200 OK wins).

    Returns an empty list if the mol2 is missing or has no recognizable
    substructure name.
    """
    if not os.path.exists(mol2_path):
        return []

    subst_name: Optional[str] = None
    try:
        with open(mol2_path) as f:
            in_atom_section = False
            for line in f:
                stripped = line.strip()
                if stripped.startswith("@<TRIPOS>ATOM"):
                    in_atom_section = True
                    continue
                if in_atom_section:
                    if stripped.startswith("@<TRIPOS>"):
                        break
                    fields = stripped.split()
                    if len(fields) >= 8 and fields[7]:
                        subst_name = fields[7]
                        break
    except OSError as e:
        logger.warning(f"Could not read mol2 for residue code extraction: {e}")
        return []

    if not subst_name:
        return []

    candidates: List[str] = []
    base = subst_name.upper()

    if 1 <= len(base) <= 3 and base.isalnum():
        candidates.append(base)

    while base and base[-1].isdigit():
        base = base[:-1]
        if 1 <= len(base) <= 3 and base.isalnum() and base not in candidates:
            candidates.append(base)

    return candidates


def get_ligand_centroid(
    mol2_path: str, padding: float = 10.0
) -> Tuple[Tuple[float, float, float], Tuple[float, float, float]]:
    """Compute docking box center and size from a reference ligand MOL2.

    The center is the centroid of all atom coordinates; the size is the
    axis-aligned bounding box of those coordinates plus ``padding`` on
    each axis.

    Raises
    ------
    ValueError
        If the file is missing, RDKit cannot parse it, or it contains
        no 3D conformer. Raising rather than returning a fallback box
        prevents silently docking against a wrong region of space.
    """
    if not os.path.exists(mol2_path):
        raise ValueError(f"Reference ligand MOL2 not found: {mol2_path}")

    mol = Chem.MolFromMol2File(mol2_path, sanitize=False)
    if mol is None:
        raise ValueError(f"RDKit could not parse MOL2 file: {mol2_path}")
    if mol.GetNumConformers() == 0:
        raise ValueError(f"MOL2 file has no 3D conformer: {mol2_path}")

    coords = mol.GetConformer().GetPositions()
    center = tuple(float(c) for c in np.mean(coords, axis=0))
    size = tuple(
        float(s) for s in np.max(coords, axis=0) - np.min(coords, axis=0) + float(padding)
    )
    return center, size
