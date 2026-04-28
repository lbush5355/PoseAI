# Contributing to PoseAI

Thank you for your interest in contributing to PoseAI! This document provides guidelines and instructions for contributing to the project.

## Getting Started

### Prerequisites
- Python 3.8 or higher
- Git
- Conda (optional, but recommended)

### Local Development Setup

1. **Fork and Clone the Repository**
   ```bash
   git clone https://github.com/YOUR_USERNAME/PoseAI.git
   cd PoseAI
   ```

2. **Create a Development Environment**
   
   Using Conda (recommended):
   ```bash
   conda env create -f environment.yml
   conda activate poseai
   ```
   
   Or using pip:
   ```bash
   pip install -r requirements.txt
   pip install pytest pytest-cov flake8 black isort
   ```

3. **Verify Installation**
   ```bash
   pytest tests/ -v
   ```

## Development Workflow

### 1. Create a Feature Branch
```bash
git checkout -b feature/your-feature-name
```

### 2. Make Your Changes
- Write clean, well-documented code
- Follow PEP 8 style guidelines
- Add docstrings to all functions and classes
- Include type hints where appropriate

### 3. Run Tests Locally
```bash
# Run all tests
pytest tests/ -v

# Run with coverage report
pytest tests/ -v --cov=src --cov-report=html

# Run specific test file
pytest tests/test_pipeline.py -v
```

### 4. Check Code Quality
```bash
# Format with black
black src tests

# Sort imports with isort
isort src tests

# Lint with flake8
flake8 src tests
```

### 5. Commit Your Changes
```bash
git add .
git commit -m "feat: add new feature description"
```

Use conventional commits:
- `feat:` for new features
- `fix:` for bug fixes
- `docs:` for documentation
- `test:` for tests
- `refactor:` for code refactoring

### 6. Push and Create a Pull Request
```bash
git push origin feature/your-feature-name
```

Then create a PR on GitHub with:
- Clear title summarizing the change
- Description of what was changed and why
- Reference to any related issues

## Code Standards

### Docstring Format
Use Google-style docstrings:
```python
def function_name(param1: str, param2: int) -> bool:
    """
    Brief description of the function.
    
    Longer description if needed, explaining the purpose and any 
    important details about the function's behavior.
    
    Args:
        param1: Description of param1
        param2: Description of param2
        
    Returns:
        Description of return value
        
    Raises:
        ValueError: When value is invalid
    """
    pass
```

### Type Hints
Include type hints for better code clarity:
```python
from typing import List, Dict, Optional

def process_data(
    data: List[Dict[str, float]],
    threshold: float = 0.5
) -> Optional[List[str]]:
    """Process data and return filtered results."""
    pass
```

## Testing Guidelines

- Write unit tests for new functions
- Maintain or improve code coverage
- Test edge cases and error conditions
- Use descriptive test names: `test_<function>_<condition>`

Example:
```python
def test_preprocessor_strips_water_molecules(self):
    """Test that water molecules are properly stripped."""
    # test implementation
    pass
```

## Module-Specific Guidelines

### preprocessor.py
- Validate input file formats
- Handle missing or corrupted files gracefully
- Log all preprocessing steps

### site_finder.py
- Document coordinate systems (PDB vs internal)
- Include validation of binding site centroids
- Add tests for edge cases (e.g., multiple ligands)

### docking.py
- Ensure subprocess calls are properly managed
- Handle engine-specific parameter variations
- Log docking progress and errors

### consensus.py
- Document clustering algorithms and parameters
- Include RMSD calculation validation
- Test symmetry-aware scoring

### visualizer.py
- Test HTML output generation
- Validate 3D coordinate transformations
- Document visualization parameters

## Reporting Issues

When reporting bugs, please include:
1. Python version and OS
2. Steps to reproduce
3. Expected vs actual behavior
4. Error messages or stack traces
5. Relevant code snippets

## Questions?

- Check existing issues and discussions
- Open a new GitHub Discussion for questions
- Contact the maintainers via email

## License

By contributing to PoseAI, you agree that your contributions will be licensed under the same license as the project.

Thank you for contributing! 🎉
