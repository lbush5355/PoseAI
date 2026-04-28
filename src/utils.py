"""PoseAI Utilities Module — Environment Management and Data Integrity.

Centralises filesystem operations, logging configuration, and dependency
verification for the PoseAI molecular-docking pipeline.  All other modules
obtain their loggers from the ``"poseai"`` namespace hierarchy established
here.

Typical usage from a Colab notebook::

    from utils import PoseAIUtils

    env = PoseAIUtils()
    env.sync_to_fast_lane("/content/drive/MyDrive/PoseAI", "/content/fast_lane")
    env.check_environment(["smina", "gnina", "ledock", "lepro", "fpocket", "obabel"])
"""

from __future__ import annotations

import logging
import os
import shutil
from datetime import datetime
from typing import List


class PoseAIUtils:
    """Environment manager and utility provider for the PoseAI pipeline.

    Responsibilities include configuring the unified logging system,
    synchronising source code and binaries to the high-speed local scratch
    disk, persisting results back to Google Drive, and verifying that all
    required computational binaries are available.

    Args:
        log_dir: Directory for rotating log files.
        log_level: Minimum severity level for both file and console handlers.
    """

    def __init__(
        self,
        log_dir: str = "/content/fast_lane/logs",
        log_level: int = logging.INFO,
    ) -> None:
        self.log_dir = log_dir
        os.makedirs(self.log_dir, exist_ok=True)
        self._setup_logging(log_level)

    # ------------------------------------------------------------------
    # Logging configuration
    # ------------------------------------------------------------------
    def _setup_logging(self, level: int) -> None:
        """Configure the ``poseai`` root logger with file and console handlers.

        Every module in the pipeline (``docking``, ``consensus``,
        ``site_finder``, ``preprocessor``, ``visualizer``) creates its
        logger via ``logging.getLogger("poseai.<module>")``.  Configuring
        the ``"poseai"`` parent here ensures that a single call controls
        formatting, level, and output destinations for the entire pipeline.
        """
        log_file = os.path.join(
            self.log_dir,
            f"poseai_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log",
        )

        formatter = logging.Formatter(
            "%(asctime)s | %(name)-24s | %(levelname)-8s | %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )

        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(formatter)

        console_handler = logging.StreamHandler()
        console_handler.setFormatter(formatter)

        # Attach handlers to the "poseai" namespace root — all child
        # loggers (poseai.docking, poseai.consensus, etc.) inherit these.
        root_logger = logging.getLogger("poseai")
        root_logger.setLevel(level)

        # Prevent duplicate handlers on repeated instantiation.
        if not root_logger.handlers:
            root_logger.addHandler(file_handler)
            root_logger.addHandler(console_handler)

        # Stop propagation to the root logger, which may carry a default
        # StreamHandler from prior basicConfig() calls or Colab's own
        # logging setup.  Without this, every message prints twice.
        root_logger.propagate = False

        self.logger = root_logger
        self.logger.info("Logging initialised.  Log file: %s", log_file)

    # ------------------------------------------------------------------
    # Filesystem synchronisation
    # ------------------------------------------------------------------
    def sync_to_fast_lane(self, drive_path: str, fast_lane_path: str) -> bool:
        """Synchronise source modules and binaries to the high-speed disk.

        Copying from a virtualised Google Drive mount to the local scratch
        disk eliminates I/O bottlenecks during execution.

        Args:
            drive_path: Source path on the mounted Google Drive.
            fast_lane_path: Destination path on the local scratch disk.

        Returns:
            ``True`` if synchronisation completed without error.
        """
        if not os.path.exists(drive_path):
            self.logger.error("Sync failed: source path does not exist: %s", drive_path)
            return False

        try:
            if os.path.isdir(drive_path):
                shutil.copytree(drive_path, fast_lane_path, dirs_exist_ok=True)
            else:
                os.makedirs(os.path.dirname(fast_lane_path), exist_ok=True)
                shutil.copy2(drive_path, fast_lane_path)
            self.logger.info("Sync complete: %s -> %s", drive_path, fast_lane_path)
            return True
        except Exception as exc:
            self.logger.error("Sync error: %s", exc)
            return False

    def persist_results(self, fast_lane_results: str, drive_results: str) -> bool:
        """Archive computational results to persistent Google Drive storage.

        Args:
            fast_lane_results: Local directory containing simulation outputs.
            drive_results: Destination directory on Google Drive.

        Returns:
            ``True`` if archival completed without error.
        """
        os.makedirs(drive_results, exist_ok=True)
        self.logger.info("Archiving results to persistent storage...")

        try:
            shutil.copytree(fast_lane_results, drive_results, dirs_exist_ok=True)
            self.logger.info("Results secured at %s", drive_results)
            return True
        except Exception as exc:
            self.logger.error("Persistence failure: %s", exc)
            return False

    # ------------------------------------------------------------------
    # Dependency verification
    # ------------------------------------------------------------------
    def check_environment(self, required_bins: List[str]) -> bool:
        """Verify that all required computational binaries are on ``$PATH``.

        Args:
            required_bins: Binary names to check
                (e.g., ``["smina", "gnina", "fpocket"]``).

        Returns:
            ``True`` if every binary is found and executable.
        """
        missing = [b for b in required_bins if shutil.which(b) is None]

        if missing:
            self.logger.warning("Missing dependencies: %s", ", ".join(missing))
            return False

        self.logger.info(
            "Environment check passed: all %d dependencies located.",
            len(required_bins),
        )
        return True

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------
    @staticmethod
    def clean_filename(name: str) -> str:
        """Sanitise a string for safe filesystem usage.

        Retains only alphanumeric characters, spaces, dots, underscores,
        and hyphens.  Spaces are replaced with underscores.
        """
        safe = "".join(
            c for c in name if c.isalnum() or c in (" ", ".", "_", "-")
        )
        return safe.strip().replace(" ", "_")
