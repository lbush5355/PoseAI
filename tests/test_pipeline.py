
"""
Comprehensive test suite for PoseAI pipeline.

Tests cover module imports, package integration, and data structure validation.
"""

import sys
import os
from pathlib import Path

import pytest

# Add src to path for imports
src_path = Path(__file__).parent.parent / "src"
if str(src_path) not in sys.path:
    sys.path.insert(0, str(src_path))


class TestImports:
    """Test that all modules can be imported successfully."""

    def test_import_config(self):
        """Test config module import."""
        try:
            import config
            assert config is not None
        except ImportError as e:
            pytest.fail(f"Failed to import config: {e}")

    def test_import_preprocessor(self):
        """Test preprocessor module import."""
        try:
            import preprocessor
            assert preprocessor is not None
        except ImportError as e:
            pytest.fail(f"Failed to import preprocessor: {e}")

    def test_import_site_finder(self):
        """Test site_finder module import."""
        try:
            import site_finder
            assert site_finder is not None
        except ImportError as e:
            pytest.fail(f"Failed to import site_finder: {e}")

    def test_import_docking(self):
        """Test docking module import."""
        try:
            import docking
            assert docking is not None
        except ImportError as e:
            pytest.fail(f"Failed to import docking: {e}")

    def test_import_consensus(self):
        """Test consensus module import."""
        try:
            import consensus
            assert consensus is not None
        except ImportError as e:
            pytest.fail(f"Failed to import consensus: {e}")

    def test_import_visualizer(self):
        """Test visualizer module import."""
        try:
            import visualizer
            assert visualizer is not None
        except ImportError as e:
            pytest.fail(f"Failed to import visualizer: {e}")

    def test_import_utils(self):
        """Test utils module import."""
        try:
            import utils
            assert utils is not None
        except ImportError as e:
            pytest.fail(f"Failed to import utils: {e}")


class TestPackageIntegration:
    """Test package-level integration and metadata."""

    def test_package_version(self):
        """Test that package has version metadata."""
        try:
            import src
            assert hasattr(src, '__version__')
            assert isinstance(src.__version__, str)
            assert len(src.__version__) > 0
        except (ImportError, AttributeError) as e:
            pytest.fail(f"Package metadata missing: {e}")

    def test_package_author(self):
        """Test that package has author metadata."""
        try:
            import src
            assert hasattr(src, '__author__')
            assert isinstance(src.__author__, str)
        except (ImportError, AttributeError) as e:
            pytest.fail(f"Author metadata missing: {e}")

    def test_package_all_exports(self):
        """Test that __all__ is properly defined."""
        try:
            import src
            assert hasattr(src, '__all__')
            assert isinstance(src.__all__, list)
            assert len(src.__all__) > 0
        except (ImportError, AttributeError) as e:
            pytest.fail(f"__all__ exports missing: {e}")


class TestPipelineIntegration:
    """Test overall pipeline integration."""

    def test_pipeline_structure(self):
        """Test that core pipeline structure exists."""
        expected_dirs = ['src', 'data', 'tests']
        for dir_name in expected_dirs:
            dir_path = Path(__file__).parent.parent / dir_name
            assert dir_path.exists(), f"Expected directory '{dir_name}' not found"
            assert dir_path.is_dir(), f"'{dir_name}' is not a directory"

    def test_src_modules_exist(self):
        """Test that all expected source modules exist."""
        expected_modules = [
            'config.py',
            'preprocessor.py',
            'site_finder.py',
            'docking.py',
            'consensus.py',
            'visualizer.py',
            'utils.py',
        ]
        src_dir = Path(__file__).parent.parent / 'src'
        for module in expected_modules:
            module_path = src_dir / module
            assert module_path.exists(), f"Module '{module}' not found in src/"
            assert module_path.is_file(), f"'{module}' is not a file"


class TestDataStructures:
    """Test data directory structure and files."""

    def test_data_directory_exists(self):
        """Test that data directory exists."""
        data_dir = Path(__file__).parent.parent / 'data'
        assert data_dir.exists(), "Data directory not found"
        assert data_dir.is_dir(), "Data is not a directory"

    def test_requirements_file_exists(self):
        """Test that requirements.txt exists."""
        req_file = Path(__file__).parent.parent / 'requirements.txt'
        assert req_file.exists(), "requirements.txt not found"
        assert req_file.is_file(), "requirements.txt is not a file"

    def test_environment_yml_exists(self):
        """Test that environment.yml exists."""
        env_file = Path(__file__).parent.parent / 'environment.yml'
        assert env_file.exists(), "environment.yml not found"
        assert env_file.is_file(), "environment.yml is not a file"

    def test_readme_exists(self):
        """Test that README.md exists."""
        readme = Path(__file__).parent.parent / 'README.md'
        assert readme.exists(), "README.md not found"
        assert readme.is_file(), "README.md is not a file"

    def test_notebook_exists(self):
        """Test that Jupyter notebook exists."""
        notebook = Path(__file__).parent.parent / 'PoseAI.ipynb'
        assert notebook.exists(), "PoseAI.ipynb not found"
        assert notebook.is_file(), "PoseAI.ipynb is not a file"


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
