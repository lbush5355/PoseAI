"""PoseAI Consensus Module — Standardized Population Clustering.

Uses a Master Topology template (SMILES) to ensure all engine outputs
share an identical chemical graph, enabling symmetry-corrected RMSD 
and HDBSCAN clustering.

v2.0 — Refactored for robustness, flexibility, and error handling.
"""

from __future__ import annotations

import contextlib
import logging
import multiprocessing
import os
import re
import subprocess
import io
from typing import Dict, List, Optional, Tuple

import pandas as pd
import numpy as np
import hdbscan
from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolAlign, AllChem

from config import get_config
from docking import DockingResult, EngineType

logger = logging.getLogger("poseai.consensus")


@contextlib.contextmanager
def _quiet_rdkit():
    """Suppress RDKit's C++ warning/error stream within the block.

    RDKit emits valence and substructure-match warnings directly to
    stderr from C++, bypassing Python's logging module. During pose
    standardization and RMSD computation, these can flood the console
    with hundreds of duplicate messages per target. This context
    manager silences them at entry and re-enables at exit.
    """
    RDLogger.DisableLog("rdApp.error")
    RDLogger.DisableLog("rdApp.warning")
    try:
        yield
    finally:
        RDLogger.EnableLog("rdApp.error")
        RDLogger.EnableLog("rdApp.warning")


# ─────────────────────────────────────────────────────────────────
# Module-level worker for parallel RMSD matrix computation.
# Must be top-level (not a staticmethod) so multiprocessing can pickle it.
# ─────────────────────────────────────────────────────────────────
def _pair_rmsd(args: Tuple[Optional[Chem.Mol], Optional[Chem.Mol]]) -> float:
    """Compute symmetry-aware heavy-atom RMSD for one pose pair.

    Returns 100.0 as a sentinel for any failure mode (None mol, RDKit
    rejection, mismatched atom counts) so the caller can blanket-fill
    the matrix without conditional branches.
    """
    mol_i, mol_j = args
    if mol_i is None or mol_j is None:
        return 100.0
    try:
        return float(rdMolAlign.GetBestRMS(mol_i, mol_j))
    except (ValueError, RuntimeError):
        if mol_i.GetNumAtoms() != mol_j.GetNumAtoms():
            return 100.0
        try:
            coords_i = mol_i.GetConformer().GetPositions()
            coords_j = mol_j.GetConformer().GetPositions()
            diff = coords_i - coords_j
            return float(np.sqrt(np.mean(np.sum(diff ** 2, axis=1))))
        except (ValueError, RuntimeError):
            return 100.0


# ─────────────────────────────────────────────────────────────────
# Main Analyzer
# ─────────────────────────────────────────────────────────────────
class ConsensusAnalyzer:
    """Consensus clustering of multi-engine docking poses.
    
    Automatically detects ligand topology for any molecule by inferring SMILES
    from the extracted ligand structure. This enables generalized clustering
    for any protein-ligand complex without hardcoded SMILES.
    
    Parameters
    ----------
    work_dir : str
        Root working directory for outputs.
    rmsd_threshold : float
        RMSD clustering threshold in angstroms.
    ligand_smiles : str, optional
        SMILES string for the ligand. Used as a fallback master template
        source when no reference structure is available.
    reference_ligand_path : str, optional
        Path to a crystal-structure ligand file (mol2/sdf/pdb). Used as
        a fallback topology source when no ideal SDF is provided. When
        used as topology source, its bond orders come from the deposited
        experimental structure (which can disagree with engine outputs
        for complex ligands — see topology_template_path).
    topology_template_path : str, optional
        Path to an RCSB ideal-SDF file (or any SDF/mol2 with canonical
        RDKit-friendly bond orders). When provided, takes priority as
        the master topology template — its chemistry-only canonical
        graph eliminates the depositor-perception-vs-engine-perception
        mismatches that block AssignBondOrdersFromTemplate on
        peptidomimetic ligands.
    """

    def __init__(
        self,
        work_dir: str = "/content/fast_lane",
        rmsd_threshold: Optional[float] = None,
        ligand_smiles: Optional[str] = None,
        reference_ligand_path: Optional[str] = None,
        topology_template_path: Optional[str] = None,
    ) -> None:
        self.work_dir = work_dir
        self.rmsd_threshold = (
            rmsd_threshold
            if rmsd_threshold is not None
            else get_config().consensus.rmsd_threshold
        )

        self.ligand_smiles = ligand_smiles
        self.reference_ligand_path = reference_ligand_path
        self.topology_template_path = topology_template_path
        self.master_ref: Optional[Chem.Mol] = None

        self.all_poses: List[Chem.Mol] = []
        self.cluster_labels: np.ndarray = np.array([])
        self._metadata: List[Dict] = []

    # ─────────────────────────────────────────────────────────────
    # Input Validation
    # ─────────────────────────────────────────────────────────────
    @staticmethod
    def _validate_smiles(smiles: str) -> None:
        """Validate SMILES string at initialization time.
        
        Raises
        ------
        ValueError
            If SMILES is invalid or produces no atoms.
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES string: {smiles}")
        if mol.GetNumAtoms() == 0:
            raise ValueError(f"SMILES produced zero atoms: {smiles}")
        logger.info(f"Validated ligand SMILES: {mol.GetNumAtoms()} atoms")

    def _infer_ligand_smiles(self, successful_results: List[DockingResult]) -> Optional[str]:
        """Infer SMILES from extracted ligand structure for generalized clustering.
        
        Loads the first successfully docked molecule and generates SMILES from
        its 3D structure. This enables the system to work with ANY ligand without
        requiring hardcoded or user-provided SMILES.
        
        Parameters
        ----------
        successful_results : List[DockingResult]
            Successfully completed docking runs.
            
        Returns
        -------
        Optional[str]
            Inferred SMILES string, or None if inference failed.
        """
        logger.info("Attempting to infer SMILES from first successful docking output...")
        
        # Try to load molecules from each engine until one succeeds
        for result in successful_results:
            try:
                mols = self._load_molecule(result.output_path, result.engine)
                
                if mols:
                    # Use first successfully loaded molecule
                    mol = mols[0]
                    
                    # Add hydrogens for more accurate SMILES
                    mol_with_h = Chem.AddHs(mol)
                    
                    # Generate SMILES from 3D structure
                    smiles = Chem.MolToSmiles(mol_with_h)
                    
                    logger.info(f"✓ Inferred SMILES from {result.engine.name}:")
                    logger.info(f"  {smiles}")
                    logger.info(f"  Atoms: {mol.GetNumAtoms()}")
                    
                    return smiles
                    
            except Exception as e:
                logger.debug(
                    f"Could not infer SMILES from {result.engine.name}: {e}"
                )
                continue
        
        logger.error("Could not infer SMILES from any docking output")
        return None

    # ─────────────────────────────────────────────────────────────
    # Molecule Loading
    # ─────────────────────────────────────────────────────────────
    def _load_molecule(self, path: str, engine: EngineType) -> List[Chem.Mol]:
        """Load molecule(s) from file, handling engine-specific formats.
        
        Parameters
        ----------
        path : str
            Path to ligand file (SDF, PDBQT, or PDB).
        engine : EngineType
            Docking engine that produced the file.
            
        Returns
        -------
        List[Chem.Mol]
            Loaded and sanitized RDKit molecules.
        """
        if not os.path.exists(path):
            logger.warning(f"Molecule file not found: {path}")
            return []
        
        try:
            if engine in (EngineType.GNINA, EngineType.SMINA):
                return self._load_sdf(path)
            elif engine == EngineType.LEDOCK:
                return self._load_dok(path)
            else:
                return self._load_via_obabel(path, engine)
        except Exception as e:
            logger.error(
                f"Failed to load molecule from {engine.name}: {e}",
                exc_info=True
            )
            return []

    @staticmethod
    def _load_sdf(path: str) -> List[Chem.Mol]:
        """Load molecules from SDF file.
        
        Returns
        -------
        List[Chem.Mol]
            Molecules with hydrogens removed, validation errors logged.
        """
        try:
            supplier = Chem.SDMolSupplier(path, removeHs=True, sanitize=False)
            mols = [mol for mol in supplier if mol is not None]
            logger.debug(f"Loaded {len(mols)} molecules from SDF: {path}")
            return mols
        except Exception as e:
            logger.error(f"SDF parsing failed for {path}: {e}")
            return []

    @staticmethod
    def _load_via_obabel(path: str, engine: EngineType) -> List[Chem.Mol]:
        """Load PDBQT via Open Babel conversion to SDF.

        Currently unused in the dispatcher (GNINA/SMINA route through
        _load_sdf, LEDOCK routes through _load_dok). Retained for
        future format additions or experimental engines.
        """
        try:
            input_fmt = "-ipdbqt" if engine == EngineType.SMINA else "-ipdb"

            result = subprocess.run(
                ["obabel", input_fmt, path, "-osdf"],
                capture_output=True,
                text=True,
                check=True,
            )

            supplier = Chem.ForwardSDMolSupplier(
                io.BytesIO(result.stdout.encode()),
                removeHs=True,
                sanitize=False,
            )
            mols = [mol for mol in supplier if mol is not None]
            logger.debug(f"Converted {len(mols)} molecules via obabel from {engine.name}")
            return mols

        except subprocess.CalledProcessError as e:
            logger.error(f"obabel conversion failed for {engine.name}: {e.stderr}")
            return []
        except Exception as e:
            logger.error(f"Unexpected error loading {engine.name} file: {e}")
            return []

    @staticmethod
    def _load_dok(path: str) -> List[Chem.Mol]:
        """Load all poses from a LeDock .dok file.

        LeDock's .dok format concatenates poses as PDB-like blocks
        separated by 'REMARK Cluster N' headers. obabel's -ipdb reader
        only catches the first model and stops, which silently dropped
        all but one pose. This splits the file at REMARK Cluster
        boundaries and converts each pose block independently.

        Returns
        -------
        List[Chem.Mol]
            All poses extracted from the file. Empty list if the file
            is missing, unparseable, or contains no recognizable poses.
        """
        try:
            with open(path, "r") as f:
                content = f.read()
        except OSError as e:
            logger.error(f"Could not read .dok file {path}: {e}")
            return []

        # Split at "REMARK Cluster" markers. The first chunk is anything
        # before the first marker (file header / preamble) and is discarded.
        chunks = re.split(r"^REMARK\s+Cluster\b", content, flags=re.MULTILINE)
        if len(chunks) < 2:
            logger.warning(
                f".dok file has no 'REMARK Cluster' markers; "
                f"falling back to single-block parse: {path}"
            )
            pose_blocks = [content]
        else:
            pose_blocks = ["REMARK Cluster" + chunk for chunk in chunks[1:]]

        mols: List[Chem.Mol] = []
        for idx, block in enumerate(pose_blocks):
            try:
                result = subprocess.run(
                    ["obabel", "-ipdb", "-osdf"],
                    input=block,
                    capture_output=True,
                    text=True,
                    check=True,
                )
            except subprocess.CalledProcessError as e:
                logger.debug(f"obabel failed on .dok pose {idx}: {e.stderr}")
                continue

            supplier = Chem.ForwardSDMolSupplier(
                io.BytesIO(result.stdout.encode()),
                removeHs=True,
                sanitize=False,
            )
            for mol in supplier:
                if mol is not None:
                    mols.append(mol)

        logger.debug(
            f"Extracted {len(mols)} poses from .dok across "
            f"{len(pose_blocks)} cluster block(s): {path}"
        )
        return mols

    # ─────────────────────────────────────────────────────────────
    # Master Template Resolution
    # ─────────────────────────────────────────────────────────────
    def _resolve_master_template(
        self, successful_results: List[DockingResult]
    ) -> Chem.Mol:
        """Build the master topology template, trying sources in priority order.

        Priority:
          1. Ideal SDF at topology_template_path. Canonical RCSB-curated
             chemistry; preferred because its bond perception is
             RDKit-friendly and matches engine outputs more reliably than
             depositor mol2 perception (especially for peptidomimetic
             ligands with quaternary nitrogens / fused aromatics).
          2. Crystal-structure ligand at reference_ligand_path. Falls back
             when no ideal SDF is provided or the fetch failed (custom
             ligand not in RCSB CCD).
          3. User-provided or RCSB-fetched SMILES via self.ligand_smiles.
          4. Stereo-stripped retry of (3).
          5. SMILES inferred from a successful docking output (last resort).

        Raises
        ------
        ValueError
            If all sources fail to produce a parseable, non-empty Mol.
        """
        # 1. Ideal SDF (preferred: canonical chemistry, RDKit-friendly)
        if self.topology_template_path:
            mol = self._try_load_reference(self.topology_template_path)
            if mol is not None:
                logger.info(
                    f"Master template from ideal SDF "
                    f"({mol.GetNumAtoms()} heavy atoms): {self.topology_template_path}"
                )
                return mol
            logger.warning(
                f"Ideal SDF unparseable: {self.topology_template_path}; "
                "falling back to crystal reference"
            )

        # 2. Crystal-structure reference
        if self.reference_ligand_path:
            mol = self._try_load_reference(self.reference_ligand_path)
            if mol is not None:
                logger.info(
                    f"Master template from crystal reference "
                    f"({mol.GetNumAtoms()} heavy atoms): {self.reference_ligand_path}"
                )
                return mol
            logger.warning(
                f"Crystal reference unparseable: {self.reference_ligand_path}; "
                "falling back to SMILES path"
            )

        # 3-4. User/RCSB SMILES with stereo-stripped retry
        if self.ligand_smiles:
            mol = self._try_parse_smiles_with_retry(self.ligand_smiles)
            if mol is not None:
                logger.info(
                    f"Master template from SMILES ({mol.GetNumAtoms()} heavy atoms)"
                )
                return mol
            logger.warning(
                "Provided SMILES unparseable; falling back to 3D inference"
            )

        # 5. Last resort: infer from a successful docking output
        inferred = self._infer_ligand_smiles(successful_results)
        if inferred:
            mol = self._try_parse_smiles_with_retry(inferred)
            if mol is not None:
                self.ligand_smiles = inferred
                logger.info(
                    f"Master template inferred from docking output "
                    f"({mol.GetNumAtoms()} heavy atoms)"
                )
                return mol

        raise ValueError(
            "Could not build master template: ideal SDF, crystal reference, "
            "SMILES, and 3D inference all failed. Provide either a "
            "topology_template_path, reference_ligand_path, or valid ligand_smiles."
        )

    @staticmethod
    def _try_load_reference(path: str) -> Optional[Chem.Mol]:
        """Load mol2/sdf/pdb reference; return heavy-atom Mol or None on failure.

        Uses sanitize=False to match the project's existing convention for
        parsing engine/depositor outputs (see _load_sdf), then attempts a
        non-strict SanitizeMol that tolerates depositor-specific perception
        quirks.
        """
        if not os.path.exists(path):
            return None

        suffix = os.path.splitext(path)[1].lower()
        try:
            if suffix == ".mol2":
                mol = Chem.MolFromMol2File(path, removeHs=True, sanitize=False)
            elif suffix == ".sdf":
                supplier = Chem.SDMolSupplier(path, removeHs=True, sanitize=False)
                mol = next((m for m in supplier if m is not None), None)
            elif suffix in (".pdb", ".pdbqt"):
                mol = Chem.MolFromPDBFile(path, removeHs=True, sanitize=False)
            else:
                logger.debug(f"Unsupported reference suffix {suffix} for {path}")
                return None
        except (ValueError, RuntimeError, OSError) as e:
            logger.debug(f"Reference load failed for {path}: {e}")
            return None

        if mol is None or mol.GetNumAtoms() == 0:
            return None

        try:
            Chem.SanitizeMol(mol)
        except (ValueError, RuntimeError) as e:
            logger.debug(
                f"Reference {path} failed strict sanitize ({e}); using as-loaded"
            )

        return mol

    @staticmethod
    def _try_parse_smiles_with_retry(smiles: str) -> Optional[Chem.Mol]:
        """Parse SMILES; retry with stereo descriptors stripped on failure.

        Stereo descriptors (@, /, \\) are common sources of SMILES parse
        failures when they don't match the actual 3D structure. Stripping
        them produces a stereo-blind topology suitable for heavy-atom
        clustering and RMSD.
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None and mol.GetNumAtoms() > 0:
            return Chem.RemoveHs(mol)

        stripped = smiles.translate(str.maketrans("", "", "@/\\"))
        if stripped != smiles:
            mol = Chem.MolFromSmiles(stripped)
            if mol is not None and mol.GetNumAtoms() > 0:
                logger.warning(
                    "SMILES parsed only after stripping stereo descriptors"
                )
                return Chem.RemoveHs(mol)

        return None

    # ─────────────────────────────────────────────────────────────
    # Consensus Analysis Pipeline
    # ─────────────────────────────────────────────────────────────
    def analyze_ensemble(self, docking_results: List[DockingResult]) -> pd.DataFrame:
        """Cluster docking poses across multiple engines.
        
        Auto-detects ligand topology by inferring SMILES from the extracted
        ligand structure, enabling generalized clustering for any molecule.
        
        Parameters
        ----------
        docking_results : List[DockingResult]
            Results from EnsembleManager.run_ensemble().
            
        Returns
        -------
        pd.DataFrame
            Cluster summary with columns:
            - Cluster: cluster ID (int)
            - Size: number of poses in cluster (int)
            - Engines: list of engines contributing poses
            - Is_Consensus: True if >=2 engines represented (bool)
            - RMSD: mean intra-cluster RMSD (float)
        """
        logger.info("Initiating multi-pose spatial consensus analysis.")
        
        # Validate inputs
        if not docking_results:
            logger.warning("No docking results provided")
            return pd.DataFrame()
        
        successful = [r for r in docking_results if r.success]
        if not successful:
            logger.error("No successful docking runs; cannot cluster")
            return pd.DataFrame()
        
        logger.info(
            f"Processing {len(successful)} successful / {len(docking_results)} total runs"
        )

        # Build master template from the most rigorous available source
        # (crystal mol2 → SMILES → stereo-stripped → 3D inference)
        self.master_ref = self._resolve_master_template(successful)
        
        # Load and standardize topologies
        processed_poses, metadata = self._load_and_standardize(successful)
        
        n = len(processed_poses)
        if n < 2:
            logger.warning(f"Fewer than 2 poses ({n}); skipping clustering")
            return pd.DataFrame()
        
        # Compute RMSD matrix
        logger.info(f"Computing RMSD matrix for {n} poses...")
        matrix = self._compute_rmsd_matrix(processed_poses)
        
        # Cluster via HDBSCAN
        logger.info("Running HDBSCAN clustering...")
        self.all_poses = processed_poses
        self._metadata = metadata
        self.cluster_labels = self._cluster_hdbscan(matrix.astype('float64'))
        
        # Generate summary
        summary_df = self._summarize_clusters(matrix, metadata)
        logger.info(f"Clustering complete: {len(summary_df)} consensus cluster(s)")
        
        return summary_df

    def _load_and_standardize(
        self, successful_results: List[DockingResult]
    ) -> Tuple[List[Chem.Mol], List[Dict]]:
        """Load and standardize topologies across engines.

        Per-pose bond-order failures are aggregated into a single summary
        log line at the end rather than one warning per pose. RDKit's
        C++ stderr stream is also silenced during this phase since its
        valence/substructure messages are duplicates of what we already
        track via fail_counts.

        Returns
        -------
        Tuple[List[Chem.Mol], List[Dict]]
            Processed molecules and metadata.
        """
        processed_poses: List[Chem.Mol] = []
        metadata: List[Dict] = []
        fail_counts: Dict[str, int] = {}
        total_counts: Dict[str, int] = {}

        with _quiet_rdkit():
            for dr in successful_results:
                mols = self._load_molecule(dr.output_path, dr.engine)
                engine = dr.engine.name

                for idx, mol in enumerate(mols):
                    total_counts[engine] = total_counts.get(engine, 0) + 1
                    try:
                        std_mol, ok = self._assign_bond_orders(mol, engine, idx)
                    except Exception as e:
                        logger.error(
                            f"Unexpected error standardizing {engine} pose {idx}: {e}",
                            exc_info=True,
                        )
                        continue

                    if not ok:
                        fail_counts[engine] = fail_counts.get(engine, 0) + 1

                    processed_poses.append(std_mol)
                    metadata.append({
                        "engine": engine,
                        "pose_index": idx,
                        "output_file": dr.output_path,
                    })

        if fail_counts:
            breakdown = ", ".join(
                f"{eng} {fail_counts[eng]}/{total_counts[eng]}"
                for eng in sorted(fail_counts)
            )
            logger.warning(
                f"Bond order standardization fell back to input topology for "
                f"{sum(fail_counts.values())}/{sum(total_counts.values())} "
                f"poses ({breakdown})"
            )

        logger.info(
            f"Standardized {len(processed_poses)} poses from "
            f"{len(successful_results)} engines"
        )
        return processed_poses, metadata

    def _assign_bond_orders(
        self, mol: Chem.Mol, engine_name: str, pose_idx: int
    ) -> Tuple[Chem.Mol, bool]:
        """Assign bond orders from master template, with heavy-atom fallback.

        Stage 1 attempts direct template assignment. Stage 2 strips hydrogens
        from both the mol and template before retrying — handles the common
        case where the engine output has implicit Hs that confuse the bond
        perception in AssignBondOrdersFromTemplate.

        Returns
        -------
        Tuple[Chem.Mol, bool]
            The standardized mol (or the original mol on total failure)
            and a flag indicating whether either stage succeeded.
            Per-pose diagnostic detail is at DEBUG level; the caller
            aggregates failures into a single summary warning.
        """
        try:
            return AllChem.AssignBondOrdersFromTemplate(self.master_ref, mol), True
        except (ValueError, RuntimeError) as e:
            logger.debug(
                f"Direct bond order assignment failed for {engine_name} pose {pose_idx}: {e}"
            )

        try:
            mol_heavy = Chem.RemoveHs(mol, sanitize=False)
            ref_heavy = Chem.RemoveHs(self.master_ref)
            return AllChem.AssignBondOrdersFromTemplate(ref_heavy, mol_heavy), True
        except (ValueError, RuntimeError) as e:
            logger.debug(
                f"Heavy-atom bond order assignment failed for {engine_name} "
                f"pose {pose_idx}: {e}; using input topology as-is"
            )

        return mol, False

    @staticmethod
    def _compute_rmsd_matrix(mols: List[Chem.Mol]) -> np.ndarray:
        """Compute pairwise Heavy-Atom RMSD matrix.

        Parallelizes across CPU cores when pose count exceeds
        cfg.rmsd_parallel_threshold. Below the threshold, the
        process-pool startup cost outweighs the savings, so we go
        sequential. RDKit C++ warnings are silenced for the duration —
        on fork-based platforms (Linux/Colab) the silencing carries
        into Pool workers via inherited state.
        """
        n = len(mols)
        matrix = np.zeros((n, n), dtype=np.float32)

        with _quiet_rdkit():
            # Pre-strip hydrogens for all molecules to ensure strict Heavy-Atom RMSD
            heavy_mols: List[Optional[Chem.Mol]] = []
            for m in mols:
                try:
                    heavy_mols.append(Chem.RemoveHs(m))
                except (ValueError, RuntimeError) as e:
                    logger.debug(f"Could not strip hydrogens for RMSD matrix entry: {e}")
                    heavy_mols.append(None)

            pairs: List[Tuple[int, int]] = [
                (i, j) for i in range(n) for j in range(i + 1, n)
            ]
            if not pairs:
                return matrix

            cfg = get_config().consensus
            threshold = cfg.rmsd_parallel_threshold
            n_workers = cfg.rmsd_n_workers
            if n_workers <= 0:
                n_workers = max(1, (os.cpu_count() or 4) // 2)

            if n < threshold or n_workers <= 1:
                rmsds = [_pair_rmsd((heavy_mols[i], heavy_mols[j])) for i, j in pairs]
            else:
                logger.debug(
                    f"Computing {len(pairs)} pairwise RMSDs across {n_workers} workers"
                )
                with multiprocessing.Pool(processes=n_workers) as pool:
                    rmsds = pool.map(
                        _pair_rmsd,
                        ((heavy_mols[i], heavy_mols[j]) for i, j in pairs),
                    )

            for (i, j), rmsd in zip(pairs, rmsds):
                matrix[i, j] = matrix[j, i] = rmsd

        return matrix

    def _cluster_hdbscan(self, matrix: np.ndarray) -> np.ndarray:
        """Cluster RMSD matrix via HDBSCAN.
        
        Parameters
        ----------
        matrix : np.ndarray
            Precomputed RMSD distance matrix.
            
        Returns
        -------
        np.ndarray
            Cluster labels (-1 for noise).
        """
        cfg = get_config().consensus
        clusterer = hdbscan.HDBSCAN(
            min_cluster_size=cfg.min_cluster_size,
            metric=cfg.hdbscan_metric,
            cluster_selection_epsilon=self.rmsd_threshold,
            allow_single_cluster=False,
        )
        labels = clusterer.fit_predict(matrix)
        
        n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
        n_noise = sum(1 for l in labels if l == -1)
        logger.info(f"HDBSCAN: {n_clusters} clusters, {n_noise} noise points")
        
        return labels

    def _summarize_clusters(
        self, matrix: np.ndarray, metadata: List[Dict]
    ) -> pd.DataFrame:
        """Generate cluster summary DataFrame.
        
        Parameters
        ----------
        matrix : np.ndarray
            RMSD matrix.
        metadata : List[Dict]
            Pose metadata.
            
        Returns
        -------
        pd.DataFrame
            Cluster summary.
        """
        summary = []
        
        for cluster_id in set(self.cluster_labels):
            if cluster_id == -1:  # Skip noise
                continue
            
            # Get indices of poses in this cluster
            idx = np.where(self.cluster_labels == cluster_id)[0]
            
            # Extract engines represented
            cluster_engines = set(metadata[i]["engine"] for i in idx)
            
            # Compute mean intra-cluster RMSD
            cluster_rmsd = np.mean(matrix[np.ix_(idx, idx)])
            
            summary.append({
                "Cluster": int(cluster_id),
                "Size": len(idx),
                "Engines": sorted(list(cluster_engines)),
                "Num_Engines": len(cluster_engines),
                "Is_Consensus": len(cluster_engines) >= 2,
                "RMSD": round(float(cluster_rmsd), 3),
            })
        
        df_sum = pd.DataFrame(summary)
        if df_sum.empty: return df_sum
        return df_sum.sort_values("Size", ascending=False).reset_index(drop=True)

    # ─────────────────────────────────────────────────────────────
    # Confidence Scoring
    # ─────────────────────────────────────────────────────────────
    def get_confidence_score(self, cluster_df: pd.DataFrame) -> float:
        """Compute confidence in consensus (0.0–1.0).
        
        Combines multi-engine agreement and cluster size. Does NOT
        assume a fixed number of engines.
        
        Parameters
        ----------
        cluster_df : pd.DataFrame
            Output from analyze_ensemble().
            
        Returns
        -------
        float
            Confidence score [0.0, 1.0].
        """
        if cluster_df.empty:
            logger.warning("No clusters to score")
            return 0.0
        
        consensus_clusters = cluster_df[cluster_df["Is_Consensus"]]
        if consensus_clusters.empty:
            logger.warning("No multi-engine consensus clusters")
            return 0.0
        
        # Select the largest consensus cluster
        best = consensus_clusters.loc[consensus_clusters["Size"].idxmax()]
        num_engines_in_cluster = best["Num_Engines"]
        cluster_size = best["Size"]
        
        cfg = get_config().consensus

        # Consensus factor: capped at 3 engines (diminishing returns beyond)
        consensus_factor = min(1.0, num_engines_in_cluster / 3.0)

        # Size factor: normalized to ideal cluster size (diminishing returns)
        size_factor = min(1.0, cluster_size / cfg.ideal_cluster_size)

        # Weighted combination
        score = (
            consensus_factor * cfg.consensus_weight +
            size_factor * cfg.cluster_size_weight
        )
        
        return round(float(score), 2)

    # ─────────────────────────────────────────────────────────────
    # Cleanup
    # ─────────────────────────────────────────────────────────────
    def __del__(self) -> None:
        """Clear pose references to release RDKit C++ objects."""
        if hasattr(self, "all_poses"):
            self.all_poses.clear()


def calculate_native_rmsd(pose_mol, reference_mol2_path, reference_smiles=None):
    """
    Calculates the strict in-place heavy-atom RMSD between a predicted pose and a native reference.
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdMolAlign
    import numpy as np

    # Load reference ligand
    raw_ref = Chem.MolFromMol2File(reference_mol2_path, sanitize=False)
    if raw_ref is None:
        raise ValueError(f"Could not load reference ligand: {reference_mol2_path}")

    # Assign bond orders
    if reference_smiles:
        template = Chem.MolFromSmiles(reference_smiles)
        if template:
            ref_mol = AllChem.AssignBondOrdersFromTemplate(template, raw_ref)
        else:
            ref_mol = raw_ref
    else:
        ref_mol = AllChem.AssignBondOrdersFromTemplate(raw_ref, raw_ref)

    # Strip hydrogens for heavy-atom RMSD
    ref_mol = Chem.RemoveHs(ref_mol)
    pose_mol = Chem.RemoveHs(pose_mol)

    # Calculate In-Place RMSD
    rmsd = rdMolAlign.CalcRMS(pose_mol, ref_mol)
    return rmsd
