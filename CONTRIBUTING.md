# Contributing to PoseAI

Thank you for your interest in contributing to PoseAI. This document covers the development workflow, coding standards, and testing requirements.

## Important: Runtime vs. Development Environment

The PoseAI pipeline runs **exclusively on Google Colab** (Linux x86-64). The docking engine binaries (Smina, Gnina, LeDock) are Linux ELF files and cannot be invoked on macOS or Windows. Local setup below is for editing and linting source code only — all functional testing must be run in Colab.

## Getting Started

### Prerequisites
- Python 3.10 or higher
- Git
- Conda (optional, but recommended)

### Local Development Setup

1. **Fork and Clone the Repository**
   ```bash
   git clone https://github.com/lbush5355/PoseAI.git
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
   pip install pytest pytest-cov flake8 black isort mypy
   ```

3. **Verify Installation**
   ```bash
   pytest tests/ -v
   ```
   Note: the test suite currently covers imports and structure only. Functional validation requires a Colab session.

## Development Workflow

### 1. Create a Feature Branch
```bash
git checkout -b feature/your-feature-name
```

### 2. Make Your Changes

Follow the coding standards below. All changes to `src/` are subject to the non-negotiable rules in the next section.

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

# Type-check with mypy
mypy src
```

### 5. Commit Your Changes
```bash
git add <specific files>
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

These rules apply to all modules in `src/` and are enforced on every change.

### Non-Negotiable Rules

- **No bare except clauses.** Always catch specific exceptions (`ValueError`, `RuntimeError`, `IOError`, etc.) with a logged diagnostic. `except Exception as e: logger.error(...)` is acceptable. `except:` and `except Exception: pass` are never acceptable.
- **No magic numbers.** All constants go in `config.py` as named dataclass fields.
- **No `os.chdir()`.** Use the `cwd=` parameter in subprocess calls instead.
- **No hardcoded paths to `/content/`.** All paths come from constructor arguments or config.
- **No duplicate imports.** Each dependency imported once at the top of the file only.
- **No silent exception swallowing.** Every caught exception must be logged at minimum.

### Type Hints

Type hints are **required** on all function signatures and return types — not optional.

```python
def process_poses(
    mols: List[Chem.Mol],
    threshold: float,
) -> Optional[pd.DataFrame]:
```

### Comments and Docstrings

Default to writing no comments. Only add a comment when the **why** is non-obvious: a hidden constraint, a subtle invariant, a workaround for a specific bug. One short line maximum — no multi-paragraph docstrings or multi-line comment blocks. Do not describe what the code does; well-named identifiers already do that.

### String Formatting

Use f-strings exclusively. No `%` formatting or `.format()`.

### Logging

Use the `poseai.*` namespace hierarchy:
- `logging.getLogger("poseai.consensus")`
- `logging.getLogger("poseai.docking")`
- etc.

Standard levels: DEBUG for internal state, INFO for stage milestones, WARNING for recoverable failures or fallbacks, ERROR for unrecoverable failures.

## Testing Guidelines

- All tests in `tests/` using pytest
- Mock all external binaries (fpocket, smina, gnina, ledock, obabel)
- Mock all RCSB network calls
- Test fixtures in `tests/fixtures/`
- Test files named `test_<module>.py`
- Every public method needs at least one test
- Every error path needs at least one test
- Use descriptive test names: `test_<function>_<condition>`

Example:
```python
def test_load_dok_returns_all_poses_for_multi_cluster_file():
    # test implementation
    pass
```

## Module-Specific Guidelines

### preprocessor.py
- Validate input file formats at system boundaries (user-supplied files, RCSB responses)
- Log all preprocessing steps at INFO level
- `get_ligand_centroid()` is the active pocket detection method — raise `ValueError` on parse failure; do not fall back silently

### site_finder.py
- `PocketAnalyzer` (fpocket-based) is implemented but not currently the active code path; `get_ligand_centroid()` in `preprocessor.py` is active
- Changes here should not affect the active pipeline unless explicitly switching the detection strategy via config

### docking.py
- Ensure all subprocess calls use `cwd=` not `os.chdir()`
- Handle engine-specific parameter variations through `PipelineParams` / `DockingConfig`
- Log docking progress and errors; never silently drop a failed engine run

### consensus.py
- Document clustering parameters and their effect on the grading scale
- The master template fallback priority order is load-bearing — changes require updating CLAUDE.md
- Test symmetry-aware RMSD via `rdMolAlign.GetBestRMS`

### visualizer.py
- Test HTML output generation
- Validate that the receptor file exists before initializing (constructor raises `FileNotFoundError` if absent)

## Reporting Issues

When reporting bugs, please include:
1. Python version and OS
2. Steps to reproduce
3. Expected vs. actual behavior
4. Error messages or stack traces
5. Relevant log output (the `poseai.*` logger output, if available)

## Questions?

- Check existing issues and discussions on GitHub
- Open a new GitHub Discussion for questions

## License

By contributing to PoseAI, you agree that your contributions will be licensed under the MIT License as described in the `LICENSE` file.
