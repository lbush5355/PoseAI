"""PoseAI Run History — Persistent per-target performance tracking.

Appends one row per pipeline run to a CSV on Drive so that statistics
accumulate across sessions. Use append_run() after each TargetResult
in Cell 5 and per_target_stats() to produce the presentation table.

Schema (run_history.csv):
  run_id               -- caller-supplied tag (e.g., "batch_20260429_01")
  timestamp            -- ISO-8601 UTC
  target_id            -- uppercase PDB code (e.g., "1IEP")
  status               -- Success | Acceptable | Poor | No Clusters | Error
  native_rmsd          -- float Å, empty on Error/No Clusters
  confidence           -- float 0-1, empty on Error/No Clusters
  cluster_size         -- int, empty on Error/No Clusters
  num_engines          -- int engines in best cluster, empty on Error/No Clusters
  error_message        -- populated only on Status=Error
  exhaustiveness       -- docking exhaustiveness knob at time of run
  poses_per_engine     -- poses requested per engine at time of run
"""

from __future__ import annotations

import csv
import logging
import os
from datetime import datetime, timezone
from typing import TYPE_CHECKING, Optional

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    from pipeline import TargetResult

logger = logging.getLogger("poseai.run_history")

_FIELDS = [
    "run_id",
    "timestamp",
    "target_id",
    "status",
    "native_rmsd",
    "confidence",
    "cluster_size",
    "num_engines",
    "error_message",
    "exhaustiveness",
    "poses_per_engine",
]


def append_run(
    result: "TargetResult",
    csv_path: str,
    run_id: str,
    exhaustiveness: Optional[int] = None,
    poses_per_engine: Optional[int] = None,
) -> None:
    """Append one result row to ``csv_path``, creating it with a header if new.

    Safe to call even when ``result.native_rmsd`` is None or
    ``result.cluster_df`` is empty — those fields are stored as empty strings
    so every row has the same column count.
    """
    cluster_size: Optional[int] = None
    num_engines: Optional[int] = None

    if result.best_pose_indices is not None and len(result.best_pose_indices) > 0:
        cluster_size = int(len(result.best_pose_indices))

    if (
        result.cluster_df is not None
        and not result.cluster_df.empty
        and result.best_cluster_id is not None
    ):
        row_match = result.cluster_df[
            result.cluster_df["Cluster"] == result.best_cluster_id
        ]
        if not row_match.empty and "Num_Engines" in row_match.columns:
            num_engines = int(row_match.iloc[0]["Num_Engines"])

    row = {
        "run_id": run_id,
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "target_id": result.target_id,
        "status": result.status,
        "native_rmsd": "" if result.native_rmsd is None else f"{result.native_rmsd:.4f}",
        "confidence": "" if result.confidence == 0.0 and result.status in ("Error", "No Clusters")
                      else f"{result.confidence:.4f}",
        "cluster_size": "" if cluster_size is None else str(cluster_size),
        "num_engines": "" if num_engines is None else str(num_engines),
        "error_message": result.error_message,
        "exhaustiveness": "" if exhaustiveness is None else str(exhaustiveness),
        "poses_per_engine": "" if poses_per_engine is None else str(poses_per_engine),
    }

    write_header = not os.path.exists(csv_path)
    try:
        with open(csv_path, "a", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=_FIELDS)
            if write_header:
                writer.writeheader()
            writer.writerow(row)
        logger.info(f"Run history: appended {result.target_id} ({result.status}) → {csv_path}")
    except OSError as e:
        logger.error(f"Failed to append run history for {result.target_id}: {e}")


def per_target_stats(csv_path: str) -> pd.DataFrame:
    """Aggregate per-target success statistics from ``csv_path``.

    Returns a DataFrame with one row per target, columns:
      target_id, Runs, Success, Acceptable, Poor, No_Clusters, Error,
      Success_Rate, Pass_Rate, Mean_RMSD_Success, Best_RMSD

    Pass_Rate = (Success + Acceptable) / Runs.
    Mean_RMSD_Success and Best_RMSD are computed over Success-status rows only;
    NaN when the target has never produced a Success result.

    Returns an empty DataFrame when ``csv_path`` does not exist.
    """
    if not os.path.exists(csv_path):
        logger.warning(f"Run history CSV not found: {csv_path}")
        return pd.DataFrame(
            columns=[
                "target_id", "Runs", "Success", "Acceptable", "Poor",
                "No_Clusters", "Error", "Success_Rate", "Pass_Rate",
                "Mean_RMSD_Success", "Best_RMSD",
            ]
        )

    try:
        df = pd.read_csv(csv_path)
    except Exception as e:
        logger.error(f"Failed to read run history from {csv_path}: {e}")
        return pd.DataFrame()

    if df.empty:
        return pd.DataFrame()

    df["native_rmsd_num"] = pd.to_numeric(df["native_rmsd"], errors="coerce")

    records = []
    for target_id, grp in df.groupby("target_id"):
        status_counts = grp["status"].value_counts()
        runs = len(grp)
        success = int(status_counts.get("Success", 0))
        acceptable = int(status_counts.get("Acceptable", 0))
        poor = int(status_counts.get("Poor", 0))
        no_clusters = int(status_counts.get("No Clusters", 0))
        error = int(status_counts.get("Error", 0))

        success_rmsds = grp.loc[grp["status"] == "Success", "native_rmsd_num"].dropna()
        mean_rmsd = float(success_rmsds.mean()) if not success_rmsds.empty else float("nan")
        best_rmsd = float(success_rmsds.min()) if not success_rmsds.empty else float("nan")

        records.append(
            {
                "target_id": target_id,
                "Runs": runs,
                "Success": success,
                "Acceptable": acceptable,
                "Poor": poor,
                "No_Clusters": no_clusters,
                "Error": error,
                "Success_Rate": round(success / runs, 3),
                "Pass_Rate": round((success + acceptable) / runs, 3),
                "Mean_RMSD_Success": round(mean_rmsd, 3) if not np.isnan(mean_rmsd) else float("nan"),
                "Best_RMSD": round(best_rmsd, 3) if not np.isnan(best_rmsd) else float("nan"),
            }
        )

    return pd.DataFrame(records).sort_values("target_id").reset_index(drop=True)
