"""PoseAI Configuration Module — Centralized Settings Management.

Consolidates all hardcoded constants, magic numbers, and configuration
parameters into a single, easily-modifiable module. Replaces scattered
constants throughout the pipeline with a unified configuration system.

v2.0 — Supports environment variable overrides and validation.
"""

from __future__ import annotations

import os
import logging
from dataclasses import dataclass, field
from typing import Optional, Dict, Any

logger = logging.getLogger("poseai.config")


# ─────────────────────────────────────────────────────────────────
# Consensus Configuration
# ─────────────────────────────────────────────────────────────────
@dataclass
class ConsensusConfig:
    """Configuration for consensus clustering analysis."""
    
    # RMSD clustering threshold (angstroms)
    rmsd_threshold: float = 2.0
    
    # Confidence score weighting
    consensus_weight: float = 0.7        # Multi-engine agreement factor
    cluster_size_weight: float = 0.3     # Cluster population factor
    
    # HDBSCAN parameters
    min_cluster_size: int = 2
    hdbscan_metric: str = "precomputed"
    
    # Ideal cluster size (poses) for scoring normalization
    ideal_cluster_size: int = 10
    
    # Default ligand SMILES (Imatinib/STI)
    # Override with ligand_smiles parameter in ConsensusAnalyzer
    default_ligand_smiles: str = (
        "Cc1ccc(cc1Nc2nccc(n2)c3ccc(cc3)CN4CCN(CC4)C)Nc5ccc(cc5)c6cc(cnc6)C"
    )
    
    def validate(self) -> bool:
        """Validate configuration parameters.
        
        Returns
        -------
        bool
            True if all parameters are valid.
            
        Raises
        ------
        ValueError
            If any parameter is invalid.
        """
        if self.rmsd_threshold <= 0:
            raise ValueError(f"rmsd_threshold must be positive, got {self.rmsd_threshold}")
        
        if not (0 < self.consensus_weight < 1):
            raise ValueError(f"consensus_weight must be in (0, 1), got {self.consensus_weight}")
        
        if not (0 < self.cluster_size_weight < 1):
            raise ValueError(f"cluster_size_weight must be in (0, 1), got {self.cluster_size_weight}")
        
        total_weight = self.consensus_weight + self.cluster_size_weight
        if not (0.99 < total_weight < 1.01):  # Allow small floating point error
            raise ValueError(
                f"Weights must sum to 1.0, got {total_weight} "
                f"(consensus={self.consensus_weight}, size={self.cluster_size_weight})"
            )
        
        if self.min_cluster_size < 2:
            raise ValueError(f"min_cluster_size must be >=2, got {self.min_cluster_size}")
        
        if self.ideal_cluster_size <= 0:
            raise ValueError(f"ideal_cluster_size must be positive, got {self.ideal_cluster_size}")
        
        logger.info("✓ Consensus configuration validated")
        return True


# ─────────────────────────────────────────────────────────────────
# Docking Configuration
# ─────────────────────────────────────────────────────────────────
@dataclass
class DockingConfig:
    """Configuration for ensemble docking execution."""
    
    # Hardware allocation (fraction of total CPUs)
    gnina_cpu_fraction: float = 0.17            # GPU does heavy lifting
    gnina_cpu_fraction_no_gpu: float = 0.25     # Fallback if no CUDA
    smina_cpu_fraction: float = 0.42
    reserved_cpus: float = 1.0                  # Min CPUs reserved for system
    
    # Timeout enforcement (seconds)
    default_engine_timeout: float = 3600.0      # 1 hour per engine
    pool_timeout_factor: float = 1.5            # Pool timeout multiplier
    
    # Process management
    stream_buffer_limit: int = 50_000           # Max log lines per engine
    poll_interval: float = 0.5                  # Seconds between future checks
    
    # Process pool configuration
    max_workers: Optional[int] = None           # Auto-detect if None
    
    def validate(self) -> bool:
        """Validate configuration parameters.
        
        Returns
        -------
        bool
            True if valid.
            
        Raises
        ------
        ValueError
            If parameters are invalid.
        """
        if not (0 < self.gnina_cpu_fraction < 1):
            raise ValueError(f"gnina_cpu_fraction must be in (0, 1), got {self.gnina_cpu_fraction}")
        
        if not (0 < self.smina_cpu_fraction < 1):
            raise ValueError(f"smina_cpu_fraction must be in (0, 1), got {self.smina_cpu_fraction}")
        
        if self.reserved_cpus < 0:
            raise ValueError(f"reserved_cpus must be non-negative, got {self.reserved_cpus}")
        
        if self.default_engine_timeout <= 0:
            raise ValueError(f"timeout must be positive, got {self.default_engine_timeout}")
        
        if self.pool_timeout_factor < 1.0:
            raise ValueError(f"pool_timeout_factor must be >=1.0, got {self.pool_timeout_factor}")
        
        if self.stream_buffer_limit <= 0:
            raise ValueError(f"stream_buffer_limit must be positive, got {self.stream_buffer_limit}")
        
        logger.info("✓ Docking configuration validated")
        return True


# ─────────────────────────────────────────────────────────────────
# Preprocessor Configuration
# ─────────────────────────────────────────────────────────────────
@dataclass
class PreprocessorConfig:
    """Configuration for structure preparation."""
    
    # Ligand validation (atoms)
    min_ligand_atoms: int = 5
    max_ligand_atoms: int = 1000
    
    # RCSB API
    rcsb_download_url: str = "https://files.rcsb.org/download"
    rcsb_timeout_s: float = 30.0
    
    # Open Babel execution
    obabel_timeout_s: float = 60.0
    
    def validate(self) -> bool:
        """Validate configuration.
        
        Returns
        -------
        bool
            True if valid.
        """
        if self.min_ligand_atoms <= 0:
            raise ValueError(f"min_ligand_atoms must be positive, got {self.min_ligand_atoms}")
        
        if self.max_ligand_atoms <= self.min_ligand_atoms:
            raise ValueError(
                f"max_ligand_atoms must be > min_ligand_atoms "
                f"({self.max_ligand_atoms} vs {self.min_ligand_atoms})"
            )
        
        if self.rcsb_timeout_s <= 0:
            raise ValueError(f"rcsb_timeout_s must be positive, got {self.rcsb_timeout_s}")
        
        if self.obabel_timeout_s <= 0:
            raise ValueError(f"obabel_timeout_s must be positive, got {self.obabel_timeout_s}")
        
        logger.info("✓ Preprocessor configuration validated")
        return True


# ─────────────────────────────────────────────────────────────────
# Site Finder Configuration
# ─────────────────────────────────────────────────────────────────
@dataclass
class SiteFinderConfig:
    """Configuration for pocket analysis."""
    
    # Search box
    default_padding: float = 10.0               # Half-width (angstroms)
    fallback_center: tuple = field(default_factory=lambda: (0.0, 0.0, 0.0))
    
    # fpocket execution
    fpocket_timeout_s: float = 300.0
    
    # Coordinate parsing
    float_pattern: str = r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?'
    
    def validate(self) -> bool:
        """Validate configuration.
        
        Returns
        -------
        bool
            True if valid.
        """
        if self.default_padding <= 0:
            raise ValueError(f"default_padding must be positive, got {self.default_padding}")
        
        if self.fpocket_timeout_s <= 0:
            raise ValueError(f"fpocket_timeout_s must be positive, got {self.fpocket_timeout_s}")
        
        if len(self.fallback_center) != 3:
            raise ValueError(f"fallback_center must have 3 coordinates, got {self.fallback_center}")
        
        logger.info("✓ Site finder configuration validated")
        return True


# ─────────────────────────────────────────────────────────────────
# Visualizer Configuration
# ─────────────────────────────────────────────────────────────────
@dataclass
class VisualizerConfig:
    """Configuration for 3D visualization."""
    
    # View dimensions (pixels)
    default_width: int = 800
    default_height: int = 600
    
    # Atom limits
    max_atoms_for_rendering: int = 50_000
    min_atoms_for_ligand: int = 5
    
    # Colors (hex)
    consensus_color: str = "#f1c40f"    # Gold
    smina_color: str = "#2ca02c"        # Green
    gnina_color: str = "#1f77b4"        # Blue
    ledock_color: str = "#d62728"       # Red
    
    def validate(self) -> bool:
        """Validate configuration.
        
        Returns
        -------
        bool
            True if valid.
        """
        if self.default_width <= 0 or self.default_height <= 0:
            raise ValueError(f"View dimensions must be positive")
        
        if self.max_atoms_for_rendering <= 0:
            raise ValueError(f"max_atoms_for_rendering must be positive")
        
        if self.min_atoms_for_ligand <= 0:
            raise ValueError(f"min_atoms_for_ligand must be positive")
        
        logger.info("✓ Visualizer configuration validated")
        return True


# ─────────────────────────────────────────────────────────────────
# Master Configuration
# ─────────────────────────────────────────────────────────────────
@dataclass
class PoseAIConfig:
    """Master configuration for the entire PoseAI pipeline.
    
    Can be overridden via environment variables (prefix: POSEAI_).
    Example: POSEAI_CONSENSUS__RMSD_THRESHOLD=2.5
    """
    
    consensus: ConsensusConfig = field(default_factory=ConsensusConfig)
    docking: DockingConfig = field(default_factory=DockingConfig)
    preprocessor: PreprocessorConfig = field(default_factory=PreprocessorConfig)
    site_finder: SiteFinderConfig = field(default_factory=SiteFinderConfig)
    visualizer: VisualizerConfig = field(default_factory=VisualizerConfig)
    
    # Logging
    log_level: int = logging.INFO
    log_dir: str = "/content/fast_lane/logs"

    def validate_all(self) -> bool:
        """Validate all sub-configurations.
        
        Returns
        -------
        bool
            True if all configurations are valid.
        """
        self.consensus.validate()
        self.docking.validate()
        self.preprocessor.validate()
        self.site_finder.validate()
        self.visualizer.validate()
        logger.info("✓ All configurations validated successfully")
        return True

    @classmethod
    def from_environment(cls) -> PoseAIConfig:
        """Load configuration from environment variables.
        
        Supports overrides via POSEAI_{MODULE}__{PARAM} format.
        Example: POSEAI_CONSENSUS__RMSD_THRESHOLD=3.0
        
        Returns
        -------
        PoseAIConfig
            Configuration with environment overrides applied.
        """
        config = cls()
        
        # Parse environment for overrides
        for key, value in os.environ.items():
            if not key.startswith("POSEAI_"):
                continue
            
            parts = key[7:].lower().split("__")  # Remove prefix and split
            if len(parts) != 2:
                continue
            
            module_name, param_name = parts
            
            # Map module names to config objects
            module_map = {
                "consensus": config.consensus,
                "docking": config.docking,
                "preprocessor": config.preprocessor,
                "site_finder": config.site_finder,
                "visualizer": config.visualizer,
            }
            
            module = module_map.get(module_name)
            if module and hasattr(module, param_name):
                # Type conversion
                field_type = type(getattr(module, param_name))
                try:
                    setattr(module, param_name, field_type(value))
                    logger.info(f"Override from env: {module_name}.{param_name} = {value}")
                except ValueError as e:
                    logger.warning(f"Failed to parse {key}={value}: {e}")
        
        return config

    def to_dict(self) -> Dict[str, Any]:
        """Export configuration as dictionary.
        
        Returns
        -------
        Dict[str, Any]
            Flattened configuration dictionary.
        """
        return {
            "consensus": vars(self.consensus),
            "docking": vars(self.docking),
            "preprocessor": vars(self.preprocessor),
            "site_finder": vars(self.site_finder),
            "visualizer": vars(self.visualizer),
            "log_level": self.log_level,
            "log_dir": self.log_dir,
        }

    def __str__(self) -> str:
        """Pretty-print configuration."""
        lines = ["PoseAI Configuration:"]
        for module_name in ["consensus", "docking", "preprocessor", "site_finder", "visualizer"]:
            module = getattr(self, module_name)
            lines.append(f"  {module_name.upper()}:")
            for key, value in vars(module).items():
                lines.append(f"    {key}: {value}")
        return "\n".join(lines)


# ─────────────────────────────────────────────────────────────────
# Global Singleton
# ─────────────────────────────────────────────────────────────────
_global_config: Optional[PoseAIConfig] = None


def get_config() -> PoseAIConfig:
    """Get or create the global PoseAI configuration.
    
    Returns
    -------
    PoseAIConfig
        The global configuration object (lazy-initialized).
    """
    global _global_config
    if _global_config is None:
        _global_config = PoseAIConfig.from_environment()
        _global_config.validate_all()
    return _global_config


def reset_config() -> None:
    """Reset the global configuration (for testing)."""
    global _global_config
    _global_config = None


# ─────────────────────────────────────────────────────────────────
# Usage Examples
# ─────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    # Initialize configuration
    config = get_config()
    
    # Access parameters
    print(f"RMSD threshold: {config.consensus.rmsd_threshold} Å")
    print(f"Default padding: {config.site_finder.default_padding} Å")
    print(f"Docking timeout: {config.docking.default_engine_timeout} s")
    
    # Override programmatically
    config.consensus.rmsd_threshold = 2.5
    config.validate_all()
    
    # Export
    print("\n" + str(config))
