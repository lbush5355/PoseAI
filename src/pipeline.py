"""PoseAI Pipeline — Single-target end-to-end orchestration.

Cells 3, 4, and 5 of PoseAI.ipynb are thin drivers over this module.
The previous setup inlined the same preprocessing → pocket → docking →
consensus → RMSD flow in two places (Cell 3 single-target, Cell 5 batch
loop). Bug fixes had to be applied twice and could be missed in one of
them — the H-strip and amino-acid-filter regressions both surfaced as
"fixed in one place, broken in the other."

run_from_rcsb() is the single-target entry point: it fetches structure
from RCSB and prepares receptor/ligand via ProteinLigandPrep, then
delegates to _run_pipeline(). run_from_local() is the batch entry point:
it accepts a PDBbind-layout dataset folder and generates PDBQT files via
obabel before delegating. Both produce a TargetResult.
"""

from __future__ import annotations

import logging
import os
from dataclasses import dataclass, field
from typing import List, Optional

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign

from consensus import ConsensusAnalyzer
from docking import DockingResult, EngineType, EnsembleManager
from preprocessor import (
    ProteinLigandPrep,
    extract_ligand_code_candidates_from_mol2,
    fetch_rcsb_entry_ligand_codes,
    fetch_rcsb_ideal_sdf,
    fetch_rcsb_smiles,
    get_ligand_centroid,
)

logger = logging.getLogger("poseai.pipeline")


@dataclass
class PipelineParams:
    """Tunable knobs shared by single-target and batch runs."""

    exhaustiveness: int
    poses_per_engine: int
    box_padding: float
    rmsd_threshold: float
    timeout: int
    engines: List[EngineType] = field(default_factory=list)


@dataclass
class TargetResult:
    """Outcome of running one target end-to-end.

    cluster_df is empty when no consensus clusters were found; native_rmsd
    is None when no clusters formed or the RMSD calculation itself failed.
    Cells consuming this object should branch on `status` rather than
    re-deriving it.

    error_message is set when status="Error" — it carries the exception
    type and message so the failure is visible on the result row (and
    therefore in run_history.csv) rather than only in stderr. Empty
    string for non-Error statuses keeps CSV columns regular.
    """

    target_id: str
    status: str  # Success | Acceptable | Poor | "No Clusters" | Error
    receptor_pdb: str
    receptor_pdbqt: str
    ligand_mol2: str
    cluster_df: pd.DataFrame
    analyzer: Optional[ConsensusAnalyzer] = None
    docking_results: List[DockingResult] = field(default_factory=list)
    top_pose_mol: Optional[Chem.Mol] = None
    best_cluster_id: Optional[int] = None
    best_pose_indices: Optional[np.ndarray] = None
    native_rmsd: Optional[float] = None
    confidence: float = 0.0
    error_message: str = ""


# ──────────────────────────────────────────────────────────────────────
# Public entry points
# ──────────────────────────────────────────────────────────────────────
def run_from_rcsb(
    target_pdb: str,
    *,
    ligand_code: Optional[str] = None,
    ligand_smiles: Optional[str] = None,
    work_dir: str,
    results_dir: str,
    params: PipelineParams,
) -> TargetResult:
    """Single-target run: fetch from RCSB, prepare, dock, score.

    When ``ligand_code`` is None (or empty), the pipeline auto-resolves
    it via fetch_rcsb_entry_ligand_codes (drug-like non-polymer entities,
    amino acids excluded). When ``ligand_smiles`` is None, it is fetched
    from RCSB chemcomp once the code is known.
    """
    target_pdb = target_pdb.lower()

    if not ligand_code:
        ligand_code = _resolve_ligand_code_from_entry(target_pdb)

    if not ligand_smiles:
        ligand_smiles = fetch_rcsb_smiles(ligand_code)
        if ligand_smiles:
            logger.info(f"Fetched SMILES for {ligand_code}: {ligand_smiles}")

    prep = ProteinLigandPrep(target_pdb)
    if not prep.fetch_structure():
        raise RuntimeError(f"Failed to fetch structure for {target_pdb.upper()}")
    if not prep.prepare_receptor():
        raise RuntimeError(f"Failed to prepare receptor for {target_pdb.upper()}")
    if not prep.isolate_ligand(ligand_code):
        raise RuntimeError(
            f"Failed to isolate ligand {ligand_code} from {target_pdb.upper()}"
        )

    ideal_sdf_path = fetch_rcsb_ideal_sdf(ligand_code, prep.work_dir)

    return _run_pipeline(
        target_id=target_pdb.upper(),
        receptor_pdbqt=prep.receptor_pdbqt,
        receptor_pdb=prep.receptor_pdb,
        ligand_pdbqt=prep.ligand_pdbqt,
        ligand_mol2=prep.ligand_mol2,
        ligand_smiles=ligand_smiles,
        ideal_sdf_path=ideal_sdf_path,
        work_dir=work_dir,
        results_dir=results_dir,
        params=params,
    )


def run_from_local(
    target_id: str,
    target_path: str,
    *,
    work_dir: str,
    results_dir: str,
    params: PipelineParams,
) -> TargetResult:
    """Batch run: PDBbind-layout local dataset folder.

    Expects ``<target_path>/<target_id>_protein.pdb`` and
    ``<target_path>/<target_id>_ligand.mol2``. PDBQT files are generated
    via obabel on demand and cached alongside.
    """
    target_id = target_id.lower()
    receptor_file = os.path.join(target_path, f"{target_id}_protein.pdb")
    ligand_file = os.path.join(target_path, f"{target_id}_ligand.mol2")

    if not os.path.exists(receptor_file) or not os.path.exists(ligand_file):
        raise RuntimeError(
            f"Missing _protein.pdb or _ligand.mol2 in {target_path}"
        )

    receptor_pdbqt = os.path.join(target_path, f"{target_id}_protein.pdbqt")
    ligand_pdbqt = os.path.join(target_path, f"{target_id}_ligand.pdbqt")
    if not os.path.exists(receptor_pdbqt):
        os.system(f"obabel {receptor_file} -O {receptor_pdbqt} -xr > /dev/null 2>&1")
    if not os.path.exists(ligand_pdbqt):
        os.system(f"obabel {ligand_file} -O {ligand_pdbqt} -p 7.4 > /dev/null 2>&1")

    # Fall back to a SMILES inferred from the local mol2 — the consensus
    # SMILES path expects something parseable when no ideal SDF is found.
    native_mol = Chem.MolFromMol2File(ligand_file, sanitize=False)
    ligand_smiles = Chem.MolToSmiles(native_mol) if native_mol else None

    ideal_sdf_path = _resolve_ideal_sdf_from_local(target_id, target_path, ligand_file)

    return _run_pipeline(
        target_id=target_id.upper(),
        receptor_pdbqt=receptor_pdbqt,
        receptor_pdb=receptor_file,
        ligand_pdbqt=ligand_pdbqt,
        ligand_mol2=ligand_file,
        ligand_smiles=ligand_smiles,
        ideal_sdf_path=ideal_sdf_path,
        work_dir=work_dir,
        results_dir=results_dir,
        params=params,
    )


# ──────────────────────────────────────────────────────────────────────
# Shared inner orchestration
# ──────────────────────────────────────────────────────────────────────
def _run_pipeline(
    *,
    target_id: str,
    receptor_pdbqt: str,
    receptor_pdb: str,
    ligand_pdbqt: str,
    ligand_mol2: str,
    ligand_smiles: Optional[str],
    ideal_sdf_path: Optional[str],
    work_dir: str,
    results_dir: str,
    params: PipelineParams,
) -> TargetResult:
    """Dock → consensus → RMSD. Inputs already prepared by caller."""
    os.makedirs(results_dir, exist_ok=True)

    center, size = get_ligand_centroid(ligand_mol2, padding=params.box_padding)

    mgr = EnsembleManager()
    docking_results = mgr.run_ensemble(
        receptor_pdbqt, ligand_pdbqt, center, size, results_dir,
        exhaustiveness=params.exhaustiveness,
        num_modes=params.poses_per_engine,
        receptor_pdb=receptor_pdb,
        ligand_mol2=ligand_mol2,
        n_ledock_poses=params.poses_per_engine,
        timeout=params.timeout,
        engines=params.engines,
    )

    # Crystal mol2 is the primary master template; ideal SDF is the automatic
    # fallback if the crystal mol2 causes >50% bond-order assignment failure
    # (see ConsensusAnalyzer.analyze_ensemble and cfg.bond_order_fallback_threshold).
    if ideal_sdf_path:
        logger.info(f"Ideal SDF available as fallback template: {ideal_sdf_path}")
    analyzer = ConsensusAnalyzer(
        work_dir=work_dir,
        rmsd_threshold=params.rmsd_threshold,
        ligand_smiles=ligand_smiles,
        reference_ligand_path=ligand_mol2,
        fallback_topology_path=ideal_sdf_path,
    )
    cluster_df = analyzer.analyze_ensemble(docking_results)

    base = dict(
        target_id=target_id,
        receptor_pdb=receptor_pdb,
        receptor_pdbqt=receptor_pdbqt,
        ligand_mol2=ligand_mol2,
        cluster_df=cluster_df if cluster_df is not None else pd.DataFrame(),
        analyzer=analyzer,
        docking_results=docking_results,
    )

    if cluster_df is None or cluster_df.empty:
        return TargetResult(status="No Clusters", **base)

    # Post-clustering extraction can fail in ways that are pose-data-specific
    # (e.g., RDKit's valence checker rejects a pose with an over-valent N from
    # an engine output, raising ValueError out of Chem.RemoveHs). Without this
    # wrap such an exception escapes run_from_local() and the caller loses the
    # cluster_df, the analyzer, and the docking_results — all useful for
    # post-mortem inspection. Catching here is NOT a silent suppression: we
    # log type(e).__name__ and the message, attach the same string to
    # TargetResult.error_message so it appears in run_history.csv, and set
    # status="Error" so _grade_rmsd's contract is preserved.
    try:
        best_cluster_id = int(cluster_df.iloc[0]["Cluster"])
        best_indices = np.where(analyzer.cluster_labels == best_cluster_id)[0]
        top_pose_mol = Chem.RemoveHs(analyzer.all_poses[best_indices[0]])

        native_rmsd = _compute_native_rmsd(analyzer, ligand_mol2, top_pose_mol)
        confidence = analyzer.get_confidence_score(cluster_df)
        status = _grade_rmsd(native_rmsd)
    except Exception as e:
        err = f"{type(e).__name__}: {e}"
        logger.error(
            f"Post-clustering pose extraction failed for {target_id} ({err}); "
            f"reporting Status=Error with cluster_df preserved for inspection",
            exc_info=True,
        )
        return TargetResult(
            status="Error",
            error_message=err,
            **base,
        )

    return TargetResult(
        status=status,
        top_pose_mol=top_pose_mol,
        best_cluster_id=best_cluster_id,
        best_pose_indices=best_indices,
        native_rmsd=native_rmsd,
        confidence=confidence,
        **base,
    )


# ──────────────────────────────────────────────────────────────────────
# Helpers
# ──────────────────────────────────────────────────────────────────────
def _resolve_ligand_code_from_entry(target_pdb: str) -> str:
    """Pick a ligand code from the PDB entry's drug-like non-polymer entities."""
    candidates = fetch_rcsb_entry_ligand_codes(target_pdb)
    if not candidates:
        raise RuntimeError(
            f"No drug-like ligand candidates found for {target_pdb.upper()}. "
            f"Provide LIGAND_CODE explicitly."
        )
    code = candidates[0]
    logger.info(f"Auto-resolved ligand code for {target_pdb.upper()}: {code}")
    return code


def _resolve_ideal_sdf_from_local(
    target_id: str, target_path: str, ligand_file: str
) -> Optional[str]:
    """Find an ideal SDF for a batch target.

    Tries RCSB GraphQL first (authoritative non-polymer codes), then
    falls back to candidates extracted from the mol2 substructure name
    when GraphQL returns nothing (multi-residue ligand entries).
    Returns the first ideal SDF that successfully downloads, or None.
    """
    candidates = fetch_rcsb_entry_ligand_codes(target_id)
    if not candidates:
        candidates = extract_ligand_code_candidates_from_mol2(ligand_file)
    for code in candidates:
        path = fetch_rcsb_ideal_sdf(code, target_path)
        if path:
            return path
    return None


def _compute_native_rmsd(
    analyzer: ConsensusAnalyzer,
    ligand_mol2: str,
    top_pose_mol: Chem.Mol,
) -> Optional[float]:
    """Heavy-atom in-place RMSD vs the native crystal structure.

    Same H-strip gotcha as ConsensusAnalyzer._try_load_reference:
    Chem.MolFromMol2File with sanitize=False does not reliably remove
    hydrogens, and an H-laden native_mol breaks the
    AssignBondOrdersFromTemplate substructure match against the heavy-
    atom-only master_ref. Force RemoveHs explicitly.

    Exception handling: AssignBondOrdersFromTemplate and CalcRMS surface
    the same underlying topology mismatch ("No sub-structure match
    found between the reference and probe mol") through different
    exception classes across RDKit versions — we have observed it as
    ValueError, RuntimeError, and Boost.Python-wrapped variants. We
    catch broadly here AND log type(e).__name__ explicitly so failures
    remain auditable: a reviewer can grep the run log to see whether a
    particular target hit AssignBondOrdersFromTemplate or CalcRMS, what
    exception type was raised, and which fallback path executed. This
    is not a silent suppression — the function still returns None on
    CalcRMS failure, which maps to Status="Error" via _grade_rmsd().
    """
    native_mol = Chem.MolFromMol2File(ligand_mol2, sanitize=False)
    if native_mol is None:
        logger.warning(f"Could not load native mol2 for RMSD: {ligand_mol2}")
        return None
    try:
        native_mol = Chem.RemoveHs(native_mol)
    except (ValueError, RuntimeError) as e:
        logger.debug(f"RemoveHs failed on native_mol ({e}); using as-loaded")

    try:
        ref_mol = Chem.RemoveHs(
            AllChem.AssignBondOrdersFromTemplate(analyzer.master_ref, native_mol)
        )
    except Exception as e:
        logger.warning(
            f"AssignBondOrdersFromTemplate failed "
            f"({type(e).__name__}: {e}); falling back to master_ref "
            f"directly for RMSD reference (RMSD will compare predicted "
            f"pose against master_ref coordinates instead of bond-order-"
            f"reconciled native_mol)"
        )
        ref_mol = analyzer.master_ref

    try:
        return float(rdMolAlign.CalcRMS(top_pose_mol, ref_mol))
    except Exception as e:
        logger.error(
            f"RMSD calculation failed ({type(e).__name__}: {e}); "
            f"returning None — target will be reported as Status=Error"
        )
        return None


def _grade_rmsd(rmsd: Optional[float]) -> str:
    """Map a heavy-atom RMSD to the project's grading scale."""
    if rmsd is None:
        return "Error"
    if rmsd <= 2.0:
        return "Success"
    if rmsd <= 3.0:
        return "Acceptable"
    return "Poor"
