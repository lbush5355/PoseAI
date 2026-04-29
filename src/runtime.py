"""PoseAI Runtime Bootstrap — Colab session setup.

Cell 1 of PoseAI.ipynb is a thin orchestration layer over this module.
The split is deliberate: anything that benefits from versioned, reviewed,
testable code lives here; the cell itself contains only the few steps
that must run before src/ is importable (Drive mount + git clone).
Once the clone completes and src/ is on sys.path, control hands off to
setup_environment().
"""

from __future__ import annotations

import logging
import os
from dataclasses import dataclass

from utils import download_and_verify_binary

logger = logging.getLogger("poseai.runtime")


@dataclass
class RuntimeContext:
    """Paths and version metadata for the active Colab runtime."""

    fast_lane: str
    src_dir: str
    bin_dir: str
    git_ref: str
    head_sha: str


# Engine binary download URLs. Each URL must serve a valid x86-64 ELF;
# download_and_verify_binary rejects HTML error pages and wrong-arch
# binaries, deletes the bad file, and raises.
_ENGINE_BINARIES = (
    ("smina",  "https://sourceforge.net/projects/smina/files/smina.static/download"),
    ("gnina",  "https://github.com/gnina/gnina/releases/download/v1.1/gnina"),
    ("ledock", "https://www.lephar.com/download/ledock_linux_x86"),
    ("lepro",  "https://www.lephar.com/download/lepro_linux_x86"),
)

_FPOCKET_INSTALL_PATH = "/usr/local/bin/fpocket"
_FPOCKET_REPO_URL = "https://github.com/Discngine/fpocket.git"


def setup_environment(
    fast_lane: str,
    *,
    git_ref: str,
    head_sha: str,
) -> RuntimeContext:
    """Idempotent post-clone setup.

    Assumes Cell 1 has already mounted Drive, ensured the fast_lane root,
    cloned the repo, copied src/* into <fast_lane>/src, and added
    <fast_lane>/src to sys.path. From there: install pip dependencies,
    install i386 support libs (historically required for LeDock), build
    fpocket if not present, and download/verify all four engine binaries.

    Re-runnable: each step short-circuits when its output is already
    present and valid, so calling this twice in one session is safe.
    """
    src_dir = os.path.join(fast_lane, "src")
    bin_dir = os.path.join(fast_lane, "bin")
    os.makedirs(bin_dir, exist_ok=True)
    os.makedirs(os.path.join(fast_lane, "logs"), exist_ok=True)

    _install_python_deps()
    _install_i386_support_libs()
    _build_fpocket_if_missing()
    _download_and_verify_all_binaries(bin_dir)

    return RuntimeContext(
        fast_lane=fast_lane,
        src_dir=src_dir,
        bin_dir=bin_dir,
        git_ref=git_ref,
        head_sha=head_sha,
    )


def _install_python_deps() -> None:
    logger.info("Installing Python dependencies (rdkit, openbabel-wheel, py3Dmol, hdbscan)...")
    os.system("pip install -q rdkit openbabel-wheel py3Dmol hdbscan -q")


def _install_i386_support_libs() -> None:
    logger.info("Installing i386 support libraries...")
    os.system("sudo dpkg --add-architecture i386")
    os.system("sudo apt-get update -qq")
    os.system(
        "sudo apt-get install -y -qq "
        "libc6:i386 libncurses5:i386 libstdc++6:i386 zlib1g:i386"
    )


def _build_fpocket_if_missing() -> None:
    if os.path.exists(_FPOCKET_INSTALL_PATH):
        logger.debug(f"fpocket already installed at {_FPOCKET_INSTALL_PATH}")
        return
    logger.info("Building fpocket from source...")
    os.system(f"git clone -q {_FPOCKET_REPO_URL} /tmp/fpocket")
    os.system("cd /tmp/fpocket && make -s && sudo make install -s")


def _download_and_verify_all_binaries(bin_dir: str) -> None:
    logger.info("Verifying docking engine binaries...")
    for name, url in _ENGINE_BINARIES:
        download_and_verify_binary(url, os.path.join(bin_dir, name), name)
    logger.info("All engine binaries verified.")
