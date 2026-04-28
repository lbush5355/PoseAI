# CLAUDE.md -- PoseAI Project Instructions

## Overview

PoseAI is an automated pipeline for ensemble molecular docking. It
integrates multiple docking engines (Gnina, Smina, LeDock) and employs
a consensus clustering strategy via HDBSCAN to identify the most probable
binding poses. Confidence is quantified via an Ensemble Confidence Score
based on spatial agreement and multi-engine representation.

Production readiness (v1.x) is defined as: the pipeline successfully
processes multiple gold standard PDB/ligand pairs and identifies the
correct consensus binding mode for each, validated against the crystal
structure via heavy-atom RMSD.

---

## Directory Structure

### Repository Layout (Claude Code's working context)

```
PoseAI/
  src/
    __init__.py        -- Package metadata and module re-exports
    config.py          -- Centralized dataclass configuration (PoseAIConfig,
                          get_config()); env-var overrides via POSEAI_*
    consensus.py       -- ConsensusAnalyzer: topology standardization,
                          pairwise RMSD matrix, HDBSCAN clustering,
                          confidence scoring. Master template resolved
                          via _resolve_master_template() in priority
                          order: topology_template_path (ideal SDF) >
                          reference_ligand_path (crystal mol2/sdf/pdb) >
                          SMILES > stereo-stripped retry > 3D-inferred.
                          calculate_native_rmsd(): standalone function,
                          strict in-place heavy-atom RMSD between predicted
                          pose and crystal structure (no alignment).
                          Returns float in angstroms.
    docking.py         -- EnsembleManager, HardwareProfile, DockingResult,
                          EngineType: parallel Smina/Gnina/LeDock execution
                          via ProcessPoolExecutor with timeout and
                          process-group signal handling
    preprocessor.py    -- ProteinLigandPrep: RCSB fetch, water removal,
                          ligand isolation by residue code (single chain),
                          PDBQT/MOL2 generation via Open Babel.
                          get_ligand_centroid(): module-level function,
                          active pocket detection method.
                          fetch_rcsb_smiles(): module-level function,
                          dynamically fetches ligand SMILES from RCSB.
                          fetch_rcsb_ideal_sdf(): downloads the RCSB
                          ideal SDF for a ligand code; cached on disk.
                          fetch_rcsb_entry_ligand_codes(): GraphQL
                          lookup of non-polymer ligand codes for a PDB
                          entry, filtered to drug-like candidates.
                          extract_ligand_code_candidates_from_mol2():
                          generates fallback candidate codes from a
                          mol2 substructure name (residue-number stripping).
    site_finder.py     -- PocketAnalyzer: fpocket wrapper, pocket barycenter
                          parsing, docking box parameters. Currently not the
                          active pocket detection code path.
    utils.py           -- PoseAIUtils: logging setup under poseai.* namespace,
                          Drive/scratch sync, PATH verification
    visualizer.py      -- DockingVisualizer: py3Dmol rendering, engine-color-
                          coded poses, HTML export
  tests/
    test_pipeline.py   -- pytest suite (currently imports and structure only)
  data/
    1iep/              -- Sample PDB target and preprocessed files
                          (tracked in repo as reference structure)
  PoseAI.ipynb         -- Colab notebook entry point
  environment.yml      -- Conda environment spec
  requirements.txt     -- Pip dependencies
  pyproject.toml       -- Package metadata and dev tooling config
  CLAUDE.md            -- This file
  README.md            -- Reproducibility-focused project readme
  CONTRIBUTING.md      -- Contribution guidelines
  LICENSE              -- MIT license
```

### What Is and Is Not Tracked in GitHub

Tracked: all src/ modules, tests/, data/ sample structures,
         notebook, environment specs, documentation
Not tracked: results/, batch_results/, dataset/, fast_lane/, bin/
             (binaries are downloaded fresh each Colab session)

### Runtime Environment (not Claude Code's working context)

The pipeline runs exclusively on Google Colab using high-performance
hardware (NVIDIA A100 GPU, high vCPU count) unavailable on a local
machine. Docking engine binaries are Linux x86_64 ELF files and cannot
run on macOS or Windows.

Development workflow is in transition:
- Previously: active development in Colab with Gemini, pushing stable
  milestones to GitHub
- Going forward: potentially shifting to local development with Claude
  Code, pushing to GitHub, pulling into Colab for execution and testing

All execution and validation must happen in Colab regardless of where
code is edited. Code is not considered stable until it passes batch
validation in the Colab runtime.

Runtime paths (for context only):
- /content/drive/MyDrive/PoseAI/src/ -- canonical source on Drive
- /content/fast_lane/                -- ephemeral high-speed scratch disk
- /content/fast_lane/bin/            -- docking engine binaries

---

## Core Modules

| Module          | Primary Class / Functions          | Responsibility                              |
|-----------------|------------------------------------|---------------------------------------------|
| config.py       | PoseAIConfig, get_config()         | All tunable parameters, env-var overrides   |
| consensus.py    | ConsensusAnalyzer,                 | Pose clustering, confidence scoring,        |
|                 | calculate_native_rmsd()            | crystal structure validation                |
| docking.py      | EnsembleManager, HardwareProfile   | Parallel engine dispatch, timeout safety    |
| preprocessor.py | ProteinLigandPrep,                 | Structure fetch, format conversion,         |
|                 | get_ligand_centroid(),             | pocket center calculation                   |
|                 | fetch_rcsb_smiles()                |                                             |
| site_finder.py  | PocketAnalyzer                     | fpocket wrapper (not active code path)      |
| utils.py        | PoseAIUtils                        | Logging, sync, dependency verification      |
| visualizer.py   | DockingVisualizer                  | 3D rendering, HTML export                   |

---

## Pipeline Execution Flow

1. **Environment Setup** (PoseAI.ipynb)
   - Mount Google Drive, sync src/ modules to Fast Lane scratch disk
   - Install Python dependencies (rdkit, openbabel-wheel, py3Dmol, hdbscan)
   - Build fpocket from source, download Linux x86_64 engine binaries
     (smina, gnina, ledock, lepro) to fast_lane/bin/
   - Verify all binaries are valid ELF files matching host architecture
   - Runs once per Colab session only

2. **Preprocessing** (preprocessor.py)
   - Fetch PDB from RCSB, strip waters, assign Gasteiger charges
   - Isolate ligand by 1-3 alphanumeric residue code (single chain only)
   - Generate receptor PDBQT (Smina/Gnina) and clean PDB (LeDock)
   - Generate ligand PDBQT (Smina/Gnina) and MOL2 (LeDock)
   - Fetch ligand SMILES dynamically via fetch_rcsb_smiles() — used
     as a fallback template source by ConsensusAnalyzer when no
     crystal mol2 is available

3. **Pocket Detection** (preprocessor.py)
   - Active method: get_ligand_centroid() in preprocessor.py
   - Calculates search box center from reference ligand coordinates
   - PocketAnalyzer (fpocket) exists in site_finder.py but is not
     the current active code path

4. **Ensemble Docking** (docking.py)
   - Dispatch Smina, Gnina, and LeDock in parallel via
     ProcessPoolExecutor
   - HardwareProfile allocates threads adaptively based on CPU count
     and CUDA availability
   - Timeout enforced at both engine and pool level via
     as_completed(timeout=pool_timeout)
   - Returns List[DockingResult] with success/failure per engine

5. **Consensus Clustering** (consensus.py)
   - Load poses from all successful engines
   - Resolve master template via _resolve_master_template():
     1. Ideal SDF at topology_template_path (canonical RCSB chemistry)
     2. Crystal mol2/sdf/pdb at reference_ligand_path
     3. User/RCSB SMILES, with stereo-stripped retry on failure
     4. SMILES inferred from a successful docking output (last resort)
   - Standardize topology against the master template
   - Compute pairwise heavy-atom RMSD matrix (O(n^2)); parallelized
     above cfg.rmsd_parallel_threshold
   - Cluster via HDBSCAN, identify multi-engine consensus clusters
   - Score confidence based on engine agreement and cluster size

6. **Validation** (consensus.py)
   - Single-target run (Cell 4 of PoseAI.ipynb): call
     calculate_native_rmsd() with explicit reference_smiles to compute
     strict in-place heavy-atom RMSD against the crystal mol2.
   - Batch validation (Cell 5 of PoseAI.ipynb): build ref_mol via
     AssignBondOrdersFromTemplate(analyzer.master_ref, native_mol).
     This combines crystal positions (from native_mol) with the
     master template's canonical bond perception (from ideal SDF when
     available, else crystal mol2). Both ref_mol and the predicted
     poses end up with isomorphic graphs, ensuring rdMolAlign.CalcRMS
     succeeds without topology-mismatch artifacts.
   - No alignment applied in either path; comparison is in-place.
   - Grading scale:
     - < 2.0 angstroms: Success (near-native pose)
     - 2.0 - 3.0 angstroms: Acceptable
     - > 3.0 angstroms: Poor (pipeline failed on this target)

7. **Visualization** (visualizer.py)
   - Render receptor and consensus cluster in py3Dmol
   - Engine-specific color coding for pose origin
   - HTML export for sharing results outside Colab

8. **Batch Validation** (PoseAI.ipynb)
   - Environment setup (Stage 1) runs once per session only
   - Stages 2-7 repeat for each PDB/ligand pair in the batch
   - Aggregates RMSD grades (Success/Acceptable/Poor) across targets
   - Results persisted to Drive before session ends
   - Pipeline is milestone-ready when majority of gold standard
     targets score Success or Acceptable with zero NaN results

---

## Known Issues and Priorities

### Critical (fix before v1.x)

1. **No functional unit tests**
   - tests/test_pipeline.py currently covers imports and structure only
   - No functional test coverage on any module

### Medium

2. **Hardcoded /3.0 in get_confidence_score() (consensus.py)**
   - consensus_factor = min(1.0, num_engines_in_cluster / 3.0)
   - Should route through config as cfg.max_engines_for_consensus
   - Inconsistent with rest of scoring which already uses get_config()

3. **Ligand validation in isolate_ligand() is shallow**
   - _validate_pdbqt_ligand() and _validate_mol2_ligand() only check
     file existence and atom line count
   - Neither calls Chem.MolFromPDBQTFile / Chem.MolFromMol2File to
     confirm RDKit can actually parse the output
   - A corrupted but non-empty file will pass validation silently

4. **calculate_native_rmsd() bond order fallback is unreliable**
   - When reference_smiles is None:
     AllChem.AssignBondOrdersFromTemplate(raw_ref, raw_ref) uses the
     molecule as its own template
   - May not correctly assign bond orders for all ligand types
   - Should log a warning or require explicit SMILES for reliable
     validation results

### Low

5. **RMSD matrix is O(n^2) with no parallelization**
   - Acceptable for current batch sizes but will bottleneck at scale

6. **Active pocket detection method is not the documented one**
   - get_ligand_centroid() in preprocessor.py is active
   - PocketAnalyzer (fpocket) in site_finder.py is documented but
     not currently called
   - Preferred method should be made explicit in config

7. **No checkpointing**
   - Pipeline crash at any stage loses all work for that target
   - Intermediate results should be saved after each stage

### Resolved

- Hardcoded STI SMILES: SMILES now fetched dynamically via
  fetch_rcsb_smiles() based on LIGAND_CODE
- Pool-level timeout: enforced via as_completed(timeout=pool_timeout)
  in docking.py
- Receptor existence check: DockingVisualizer raises FileNotFoundError
  on init if receptor is missing
- 17 bare except clauses in consensus.py: replaced with specific
  exceptions (ValueError, RuntimeError) and logged diagnostics. The
  5-level fallback pyramid in _load_and_standardize was flattened to
  a 2-stage helper (_assign_bond_orders) since the deeper levels
  were no-ops -- Chem.RemoveHs is idempotent. The dead "del mol"
  loop in __del__ was removed.
- preprocessor.py "Dynamically added" tail block: fetch_rcsb_smiles()
  and get_ligand_centroid() are now proper module-level functions
  with type hints and namespace logging. Imports moved to top of file
  (numpy, requests added; duplicate Chem re-import removed). Both
  silent except-pass clauses replaced with specific exceptions and
  logged diagnostics. get_ligand_centroid() now raises ValueError on
  parse failure instead of silently returning a (0,0,0) box -- the
  prior fallback could mask docking-against-wrong-region failures as
  Poor RMSD results in batch validation.
- Hardcoded len == 3 ligand-code check in isolate_ligand: relaxed to
  accept 1-3 alphanumeric characters per the PDB chemical component
  dictionary spec (unblocks N3, ZN, MG, CA, etc.).
- fetch_rcsb_smiles() field name bug: RCSB chemcomp API returns a
  dict (not list) with uppercase keys SMILES_stereo and SMILES; code
  was looking for lowercase smilesstereo and smiles, so the fetch
  always returned None even for well-known ligands like STI.
- SMILES robustness redesign in ConsensusAnalyzer: added
  reference_ligand_path constructor param. New _resolve_master_template
  tries crystal mol2/sdf/pdb first, then SMILES (with stereo-stripped
  retry), then 3D inference as last resort. Eliminates the
  SMILES-vs-3D-mismatch class of RDKit failures (tautomers, aromaticity
  perception, charge states) for pose-validation runs by using the
  deposited experimental structure directly. Notebook Cell 5 RMSD
  comparison simplified to use analyzer.master_ref instead of the
  AssignBondOrdersFromTemplate(native_mol, native_mol) self-template
  trick, giving deterministic same-source comparison.
- Engine imbalance in ensemble docking: SMINA --energy_range default
  of 3 kcal/mol silently dropped most poses at high exhaustiveness
  (now passes --energy_range 10 explicitly via DockingConfig.SMINA_ENERGY_RANGE);
  LeDock .dok output was parsed by obabel -ipdb which only catches the
  first MODEL (now split at "REMARK Cluster" markers in _load_dok and
  each pose block is converted independently). Both fixes restore
  genuine three-engine consensus in batch validation.
- RMSD matrix parallelization: _compute_rmsd_matrix uses
  multiprocessing.Pool above cfg.rmsd_parallel_threshold. New
  rmsd_parallel_threshold and rmsd_n_workers fields on ConsensusConfig.
- Log noise reduction: aggregated per-pose bond-order failures into a
  single summary WARNING with engine breakdown counts; new
  _quiet_rdkit() context manager silences RDKit's C++ stderr stream
  during _load_and_standardize and _compute_rmsd_matrix.
- Ideal-SDF topology source: new fetch_rcsb_ideal_sdf() and
  fetch_rcsb_entry_ligand_codes() (GraphQL lookup) helpers in
  preprocessor.py; ConsensusAnalyzer accepts topology_template_path
  which takes highest priority in master template resolution.
  Addresses bond-perception mismatch on complex peptidomimetic
  ligands (1HSG/MK1, 1IEP/STI) where obabel's mol2 perception
  conflicts with engine-output perception, blocking
  AssignBondOrdersFromTemplate for 100% of poses. Cell 5 RMSD
  reference rebuilt as AssignBondOrdersFromTemplate(master_ref,
  native_mol) so crystal positions are kept while topology comes
  from the canonical ideal SDF. Cell 5 ligand-code resolution is
  now GraphQL-primary with mol2-candidate fallback (see
  extract_ligand_code_candidates_from_mol2). Master-template-source
  log line is at WARNING level for diagnostic visibility in Colab.
- fpocket coordinate regex: robust scientific notation pattern sourced
  from config

---

## Coding Standards

These apply to all modules in src/. Claude Code must enforce these
on every change it makes and flag violations it encounters in existing
code.

### Non-Negotiable Rules

- **No bare except clauses.** Always catch specific exceptions
  (ValueError, RuntimeError, IOError, etc.) with a logged diagnostic.
  `except Exception as e: logger.error(...)` is acceptable.
  `except:` and `except Exception: pass` are never acceptable.
- **No magic numbers.** All constants go in config.py as named
  dataclass fields or class-level constants.
- **No os.chdir().** Use cwd= parameter in subprocess calls instead.
- **No hardcoded paths to /content/.** All paths come from constructor
  arguments or config.
- **No duplicate imports.** Each dependency imported once at the top
  of the file only.
- **No silent exception swallowing.** Every caught exception must be
  logged at minimum.

### Code Style

- **Type hints** on all function signatures and return types.
- **f-strings** for all string formatting. No % or .format().

### Logging

Use the poseai.* namespace hierarchy throughout:
- logging.getLogger("poseai.consensus")
- logging.getLogger("poseai.docking")
- logging.getLogger("poseai.preprocessor")
- etc.

Standard levels:
- DEBUG: internal state, subprocess output
- INFO: stage milestones, major method entry/exit
- WARNING: recoverable failures, fallbacks activated
- ERROR: unrecoverable failures, caught exceptions

### Testing

- All tests in tests/ using pytest
- Mock all external binaries (fpocket, smina, gnina, ledock, obabel)
- Mock all RCSB network calls
- Test fixtures in tests/fixtures/
- Test files named test_<module>.py
- Every public method needs at least one test
- Every error path needs at least one test

---

## External Dependencies

### Python (see requirements.txt and environment.yml)

Core: rdkit, openbabel-wheel, py3Dmol, hdbscan, numpy, pandas, requests
Dev:  pytest, pytest-cov

### Docking Engine Binaries

These are Linux x86_64 ELF binaries. They cannot run on macOS or
Windows and are never committed to the repository (bin/ is gitignored).
Downloaded fresh each Colab session by PoseAI.ipynb setup cells.

| Binary  | Source                                   | Purpose              |
|---------|------------------------------------------|----------------------|
| fpocket | Built from github.com/Discngine/fpocket  | Pocket detection     |
| smina   | SourceForge (smina.static)               | Force-field docking  |
| gnina   | github.com/gnina/gnina releases v1.1     | CNN-scored docking   |
| ledock  | lephar.com/download/ledock_linux_x86     | Stochastic docking   |
| lepro   | lephar.com/download/lepro_linux_x86      | LeDock preprocessing |
| obabel  | pip install openbabel-wheel              | Format conversion    |

### Binary Verification

After download, each binary is verified as a valid ELF file matching
the host architecture. An HTML error page or ARM binary will be caught
before the pipeline runs. Use utils.py check_environment() to verify
all binaries are on PATH before executing.

---

## Gold Standard Validation Targets

### Dataset

PDBbind demo set. Files sit in dataset/ (gitignored -- not in repo).
Must be present in the Colab environment before batch validation runs.

### File Structure

Each target has its own subdirectory:

```
dataset/
  1ett/
    1ett_ligand.mol2
    1ett_ligand.pdbqt
    1ett_ligand.pdb
    1ett_protein.pdb
    1ett_protein.pdbqt
    1ett_pocket.pdb
```

### Current Batch Targets

1ett, 1owe, 1a30, 1stp, 1hxw, 1iep, 1fjs, 1hsg

### Grading Scale

- < 2.0 angstroms RMSD: Success (near-native pose)
- 2.0 - 3.0 angstroms RMSD: Acceptable
- > 3.0 angstroms RMSD: Poor (pipeline failed on this target)
- NaN: Pipeline error (crash, failed docking, or clustering failure)

### Current Benchmark (last recorded run, 2026-04-28)

- Success:    4/8 (1FJS 1.19 Å, 1STP 0.62 Å, 1OWE 0.23 Å, 1ETT 0.67 Å)
- Poor:       2/8 (1HXW 4.24 Å, 1A30 7.78 Å)
- NaN:        2/8 (1HSG, 1IEP — bond-perception mismatch; ideal-SDF
              fix shipped in commit 9e10d35 but not yet validated)
- Run was at POSES_PER_ENGINE=20 and EXHAUSTIVENESS=8 (faster batch
  defaults). 1ETT recovered from Poor (4.00 → 0.67) after the ideal-SDF
  topology change, suggesting other Poor targets may also benefit.

### Known Limitation: Multi-Residue Ligands

Targets where the ligand is deposited as a polymer chain (peptidic
inhibitor) rather than a non-polymer entity cannot use the ideal-SDF
topology path — RCSB's GraphQL endpoint returns null
nonpolymer_entities. 1A30 is the canonical example in the current
batch. These targets fall back to the crystal mol2 master template
and are subject to the same bond-perception risks the ideal SDF was
designed to mitigate. Supporting them properly requires composing the
inhibitor from its constituent residue components — out of scope for
v1.x.

### v1.x Release Criteria

Pipeline is considered production-ready when:
- Minimum 6/8 targets score Success or Acceptable
- Zero NaN results (no unhandled pipeline errors)
- All Critical known issues resolved
- tests/test_pipeline.py contains functional coverage of all modules
