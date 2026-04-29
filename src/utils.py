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
import shlex
import shutil
from datetime import datetime
from typing import List, Tuple

logger = logging.getLogger("poseai.utils")

# x86-64 ELF e_machine value (EM_X86_64). Used to reject binaries built for
# the wrong architecture (a real risk if Colab assigns a non-x86_64 host).
_EM_X86_64 = 0x3E


def _read_elf_header(path: str) -> Tuple[bytes, int, int]:
    """Return (magic, ei_class, e_machine) from a candidate ELF file.

    Raises OSError if the file is shorter than a full ELF header (52 bytes
    minimum for ELF32, 64 for ELF64; 20 is enough to extract the fields we
    need: magic, EI_CLASS, e_machine).
    """
    with open(path, "rb") as f:
        header = f.read(20)
    if len(header) < 20:
        raise OSError(f"file too short to be an ELF: {len(header)} bytes")
    return header[:4], header[4], int.from_bytes(header[0x12:0x14], "little")


def _looks_like_x86_64_elf(path: str) -> bool:
    """True iff path is an ELF with e_machine == EM_X86_64."""
    try:
        magic, _ei_class, e_machine = _read_elf_header(path)
    except OSError:
        return False
    return magic == b"\x7fELF" and e_machine == _EM_X86_64


def _remove_if_exists(path: str) -> None:
    try:
        os.remove(path)
    except FileNotFoundError:
        pass


def download_and_verify_binary(
    url: str, dest_path: str, name: str, *, skip_if_valid: bool = True
) -> None:
    """Download a binary and verify it is an x86-64 ELF before keeping it.

    Cell 1's previous wget calls would silently write whatever the server
    returned — HTML error pages, rate-limit responses, partial content.
    chmod +x then ran on garbage and the existence check happily skipped
    re-download on every subsequent session, leaving the runtime in a
    broken state until manual cleanup. This helper closes that loophole:
    on any failure mode (bad rc, short file, wrong magic, wrong arch),
    it deletes the file and raises a clear error.

    Parameters
    ----------
    url : str
        Source URL. Must support ``-L`` redirects (sourceforge, github).
    dest_path : str
        Filesystem destination.
    name : str
        Human-readable name for log lines and error messages.
    skip_if_valid : bool
        When True (default) and ``dest_path`` already holds a valid x86-64
        ELF, no download is attempted. A present-but-corrupt file is always
        re-downloaded regardless of this flag.
    """
    if skip_if_valid and os.path.exists(dest_path) and _looks_like_x86_64_elf(dest_path):
        logger.debug(f"{name} already present and valid: {dest_path}")
        return

    if os.path.exists(dest_path):
        logger.warning(f"{name} present but invalid; re-downloading from {url}")
        _remove_if_exists(dest_path)

    logger.info(f"Downloading {name} from {url}")
    rc = os.system(
        f"wget -q -L {shlex.quote(url)} -O {shlex.quote(dest_path)}"
    )
    if rc != 0:
        _remove_if_exists(dest_path)
        raise RuntimeError(f"{name} download failed (wget rc={rc}): {url}")

    try:
        magic, _ei_class, e_machine = _read_elf_header(dest_path)
    except OSError as e:
        _remove_if_exists(dest_path)
        raise RuntimeError(
            f"{name} download produced no usable file ({e}): {url}"
        )

    if magic != b"\x7fELF":
        _remove_if_exists(dest_path)
        raise RuntimeError(
            f"{name} is not an ELF binary (magic={magic!r}). Source URL "
            f"likely returned an HTML error or rate-limit page: {url}"
        )

    if e_machine != _EM_X86_64:
        _remove_if_exists(dest_path)
        raise RuntimeError(
            f"{name} is wrong architecture (e_machine=0x{e_machine:x}); "
            f"expected x86-64 (0x{_EM_X86_64:x}): {url}"
        )

    os.chmod(dest_path, 0o755)
    logger.info(f"{name} verified: x86-64 ELF at {dest_path}")


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
