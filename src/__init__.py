"""
PoseAI: A Multi-Engine Consensus Molecular Docking Pipeline

A modular computational framework for executing and harmonizing ligand-binding
simulations across multiple docking scoring functions.

Modules:
    - preprocessor: Structural retrieval and formatting
    - site_finder: Binding site identification
    - docking: Multi-threaded docking orchestration
    - consensus: Pose analysis and clustering
    - visualizer: 3D visualization generation
    - config: Pipeline configuration management
    - utils: Utility functions

Version: 1.0.0
Author: Liam Bush
Email: liam.bush@ucdenver.edu
"""

__version__ = "1.0.0"
__author__ = "Liam Bush"
__email__ = "liam.bush@ucdenver.edu"

# Import core modules for easy access
try:
    from src import (
        config,
        consensus,
        docking,
        preprocessor,
        site_finder,
        utils,
        visualizer,
    )
    __all__ = [
        "config",
        "consensus",
        "docking",
        "preprocessor",
        "site_finder",
        "utils",
        "visualizer",
    ]
except ImportError:
    __all__ = []
