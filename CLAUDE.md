# CLAUDE.md -- PoseAI Project Instructions

## What This File Is

This is a [Claude Code](https://claude.ai/code) configuration file. Claude Code
is an AI-assisted development tool; this file provides it with persistent project
context — architecture, module responsibilities, coding standards, known issues,
and benchmark state — so that AI-assisted development remains coherent across
sessions. It is committed to the repository as the single source of truth for
that context. It is not user documentation; see README.md for that.

This project was developed for an AI in Chemistry & Biochemistry course at the
University of Colorado Denver, where AI-assisted development was an explicit
component of the coursework.

---

## Overview

PoseAI is an automated pipeline for ensemble molecular docking. It
integrates multiple docking engines (Gnina, Smina, LeDock) and employs
a consensus clustering strategy via HDBSCAN to identify the most probable
binding poses. Confidence is quantified via an Ensemble Confidence Score
based on spatial agreement and multi-engine representation.

Current version: v0.9 pre-release. Production readiness (v1.x) is
defined as: the pipeline successfully processes multiple gold standard
PDB/ligand pairs and identifies the correct consensus binding mode for
each, validated against the crystal structure via heavy-atom RMSD.

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
                          order: reference_ligand_path (crystal mol2/
                          sdf/pdb) > SMILES > stereo-stripped retry >
                          3D-inferred. If the crystal mol2 causes a
                          bond-order fail rate above
                          cfg.bond_order_fallback_threshold (default
                          0.5), analyze_ensemble() swaps master_ref to
                          fallback_topology_path (ideal SDF) and re-runs
                          standardization in-memory without re-docking.
                          _load_and_standardize() returns a 3-tuple
                          (poses, metadata, fail_rate). Cluster ranking
                          is deterministic: sort by Num_Engines DESC,
                          Size DESC, intra-RMSD ASC.
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
    runtime.py         -- Colab session bootstrap. setup_environment()
                          installs pip deps, i386 libs, builds fpocket if
                          missing, and downloads + ELF-validates engine
                          binaries. Cell 1 of PoseAI.ipynb is a thin
                          orchestration layer over this module.
    site_finder.py     -- PocketAnalyzer: fpocket wrapper, pocket barycenter
                          parsing, docking box parameters. Currently not the
                          active pocket detection code path.
    utils.py           -- PoseAIUtils: logging setup under poseai.* namespace,
                          Drive/scratch sync, PATH verification.
                          download_and_verify_binary(): module-level
                          function used by runtime.py to download an
                          engine binary and reject anything that isn't a
                          valid x86-64 ELF.
    pipeline.py        -- run_from_rcsb(): single-target entry point
                          (RCSB fetch + ProteinLigandPrep). run_from_local():
                          batch entry point (PDBbind-layout folder, obabel
                          PDBQT generation). Both delegate to _run_pipeline()
                          which owns docking → consensus → RMSD. Returns
                          TargetResult (status, native_rmsd, confidence,
                          cluster_df, analyzer, docking_results,
                          error_message).
    run_history.py     -- append_run(): appends one TargetResult row to
                          run_history.csv on Drive (run_id, timestamp,
                          target_id, status, native_rmsd, confidence,
                          cluster_size, num_engines, error_message,
                          exhaustiveness, poses_per_engine). Called by
                          Cell 5 after each target. per_target_stats():
                          aggregates per-target Success_Rate, Pass_Rate,
                          Mean_RMSD_Success, Best_RMSD across all recorded
                          runs; displayed at end of Cell 5 for presentation.
    visualizer.py      -- DockingVisualizer: py3Dmol rendering, engine-color-
                          coded poses, HTML export
  PoseAI.ipynb         -- Colab notebook entry point
  requirements.txt     -- Pip dependencies
  pyproject.toml       -- Package metadata and dev tooling config
  CLAUDE.md            -- This file
  README.md            -- Project documentation
  LICENSE              -- MIT license
```

### What Is and Is Not Tracked in GitHub

Tracked: all src/ modules, notebook, requirements.txt,
         pyproject.toml, documentation
Not tracked: results/, batch_results/, dataset/, fast_lane/, bin/
             (binaries are downloaded fresh each Colab session)

### Runtime Environment (not Claude Code's working context)

The pipeline runs exclusively on Google Colab using high-performance
hardware (NVIDIA A100 GPU, high vCPU count) unavailable on a local
machine. Docking engine binaries are Linux x86_64 ELF files and cannot
run on macOS or Windows.

Development workflow:
- Active development is local with Claude Code, committed to GitHub.
- Stage 1 of PoseAI.ipynb clones lbush5355/PoseAI fresh into
  /tmp/poseai_repo on every Colab session and copies src/*.py into
  /content/fast_lane/src. GitHub is the single source of truth — Drive's
  src/ is no longer used.
- All execution and validation must happen in Colab regardless of where
  code is edited. Code is not considered stable until it passes batch
  validation in the Colab runtime.

Notebook architecture:
- Cells 1, 3, and 4 are thin orchestration drivers. The pipeline logic
  they invoke lives in src/runtime.py (Stage 1 setup) and src/pipeline.py
  (Stage 3-4 single-target run; created in Phase 3 of the streamlining
  refactor, in progress). Cell 5 (batch validation) is intentionally an
  exception — it is a test harness, not a v1.x release artifact, and may
  inline batch-specific logic. To prevent single-target and batch flows
  from drifting apart, both must call the same pipeline helpers under
  the hood.

Runtime paths (for context only):
- /tmp/poseai_repo/                  -- fresh git clone, source of src/ modules
- /content/fast_lane/src/            -- working copy of src/ modules (synced from clone)
- /content/fast_lane/                -- ephemeral high-speed scratch disk
- /content/fast_lane/bin/            -- docking engine binaries
- /content/drive/MyDrive/PoseAI/     -- Drive root, used for dataset/, results/,
                                        batch_results/ persistence (not src/)

---

## Core Modules

| Module          | Primary Class / Functions          | Responsibility                              |
|-----------------|------------------------------------|---------------------------------------------|
| config.py       | PoseAIConfig, get_config()         | All tunable parameters, env-var overrides   |
| consensus.py    | ConsensusAnalyzer,                 | Pose clustering, confidence scoring,        |
|                 | calculate_native_rmsd()            | crystal structure validation                |
| docking.py      | EnsembleManager, HardwareProfile   | Parallel engine dispatch, timeout safety    |
| pipeline.py     | run_from_rcsb(), run_from_local(), | End-to-end orchestration; single-target     |
|                 | PipelineParams, TargetResult       | and batch entry points                      |
| preprocessor.py | ProteinLigandPrep,                 | Structure fetch, format conversion,         |
|                 | get_ligand_centroid(),             | pocket center calculation                   |
|                 | fetch_rcsb_smiles()                |                                             |
| run_history.py  | append_run(), per_target_stats()   | Persistent per-target run statistics CSV    |
| runtime.py      | setup_environment(),               | Cell 1 bootstrap: deps, fpocket, ELF-       |
|                 | RuntimeContext                     | validated engine binaries                   |
| site_finder.py  | PocketAnalyzer                     | fpocket wrapper (not active code path)      |
| utils.py        | PoseAIUtils,                       | Logging, sync, dependency verification,     |
|                 | download_and_verify_binary()       | x86-64 ELF download with validation         |
| visualizer.py   | DockingVisualizer                  | 3D rendering, HTML export                   |

---

## Pipeline Execution Flow

1. **Environment Setup** (PoseAI.ipynb Cell 1 + runtime.py)
   - Cell 1 mounts Drive, clones the repo from GitHub into
     /tmp/poseai_repo, copies src/* into /content/fast_lane/src, and
     adds src/ to sys.path (these steps must stay inline because
     runtime.py is itself part of the cloned src/).
   - Cell 1 then calls runtime.setup_environment(), which installs pip
     dependencies, installs i386 support libs, builds fpocket if
     missing, and downloads each engine binary via
     utils.download_and_verify_binary() — rejecting HTML error pages
     and wrong-architecture binaries before they hit subprocess.Popen.
   - Each step is idempotent; safe to re-run within a session.

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
     1. Crystal mol2/sdf/pdb at reference_ligand_path (primary)
     2. User/RCSB SMILES, with stereo-stripped retry on failure
     3. SMILES inferred from a successful docking output (last resort)
   - Standardize topology against the master template
   - If bond-order fail rate > cfg.bond_order_fallback_threshold (0.5),
     swap master_ref to fallback_topology_path (ideal SDF) and re-run
     standardization in-memory — no re-docking required (~2 s)
   - Compute pairwise heavy-atom RMSD matrix (O(n^2)); parallelized
     above cfg.rmsd_parallel_threshold
   - Cluster via HDBSCAN, identify multi-engine consensus clusters
   - Rank clusters: Num_Engines DESC, Size DESC, intra-RMSD ASC
   - Score confidence based on engine agreement and cluster size

6. **Validation** (pipeline.py + consensus.py)
   - Handled inside pipeline._run_pipeline() via _compute_native_rmsd().
   - Builds ref_mol via AssignBondOrdersFromTemplate(analyzer.master_ref,
     native_mol), combining crystal coordinates with the canonical bond
     perception from whichever master template was selected.
   - Post-clustering block is wrapped in try/except so valence errors
     from engine-output poses (e.g., over-valent N) set Status=Error
     and populate TargetResult.error_message rather than crashing the
     caller; cluster_df and analyzer are preserved for inspection.
   - No alignment applied; comparison is strict in-place heavy-atom RMSD.
   - Grading scale:
     - < 2.0 angstroms: Success (near-native pose)
     - 2.0 - 3.0 angstroms: Acceptable
     - > 3.0 angstroms: Poor (pipeline failed on this target)

7. **Visualization** (visualizer.py)
   - Render receptor and consensus cluster in py3Dmol
   - Engine-specific color coding for pose origin
   - HTML export for sharing results outside Colab

8. **Batch Validation** (PoseAI.ipynb Cell 5 + pipeline.py)
   - Environment setup (Stage 1) runs once per session only
   - Cell 5 calls pipeline.run_from_local() per target and calls
     run_history.append_run() after each result, writing one row to
     run_history.csv on Drive for cross-session tracking
   - Displays per_target_stats() (Success_Rate, Pass_Rate, Best_RMSD)
     after the batch table so cumulative performance is visible without
     leaving the notebook
   - Results and HTML overlays persisted to batch_results/ on Drive
   - Pipeline is milestone-ready when majority of gold standard
     targets score Success or Acceptable with zero NaN results

---

## Known Issues and Priorities

### Critical (fix before v1.x)

1. **No functional unit tests**
   - No test suite exists; all pipeline correctness is validated through
     the Colab batch run against gold-standard targets
   - Post-v1.x goal: build a pytest suite under tests/ with mocked
     external binaries and RCSB calls, covering at minimum consensus.py
     cluster ranking, _load_dok() multi-pose recovery, Class C retry
     trigger, and _grade_rmsd() boundary conditions

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

4. **_compute_native_rmsd() ref_mol construction fails for complex ligands
   (post-presentation project)**
   - When the crystal mol2 topology is incompatible with the ideal SDF
     (the same condition that triggers Class C auto-retry), both
     AssignBondOrdersFromTemplate(master_ref, native_mol) and the CalcRMS
     substructure match fail. The fallback ref_mol = analyzer.master_ref
     is the ideal SDF at RCSB canonical coordinates, not the crystal
     binding site — making the RMSD scientifically meaningless even if
     CalcRMS could run. Affects 1HSG and 1IEP in the current batch.
   - Root cause: the crystal mol2 for complex peptidomimetic/flexible
     ligands has different aromatic/charge perception than the ideal SDF,
     blocking template assignment in both directions.
   - Fix plan (post-presentation): when AssignBondOrdersFromTemplate fails,
     attempt to recover native coordinates by overlaying the crystal mol2
     atom positions onto the ideal SDF graph via atom-by-atom coordinate
     copy after a maximum common substructure match; fall back to
     Status=Error with HTML overlay as the validation path if that also
     fails. Until fixed, 1HSG/1IEP report Status=Error with note that
     clustering succeeded — use the generated HTML overlay for manual
     validation.

### Low

5. **No checkpointing**
   - Pipeline crash at any stage loses all work for that target
   - Intermediate results should be saved after each stage

---

## Post-v1.x Scientific Goals

These represent the primary directions for advancing PoseAI beyond
redocking validation into genuine scientific contribution.

### De Novo Binding Site Discovery (highest priority)
   - Current pipeline requires a co-crystallized ligand to define the
     docking box via get_ligand_centroid(). This restricts it to
     redocking against known structures.
   - Plan: use fpocket (already built and installed in Cell 1) to find
     candidate pockets from receptor surface geometry alone. Match
     ligand pharmacophoric features (HBD/HBA/hydrophobic/aromatic via
     RDKit MolChemicalFeatures) against each pocket's chemical
     descriptors to select the best pocket before docking.
   - PocketAnalyzer in site_finder.py already implements the fpocket
     wrapper and output parser. The pharmacophore-to-pocket matching
     layer is the remaining work.
   - Do NOT remove or refactor site_finder.py — it is the foundation
     for this capability.

### Pharmacophore-Aware Confidence Scoring
   - Current ECS measures only spatial convergence across engines.
     It cannot distinguish a geometrically convergent pose that
     satisfies key interactions from one that does not.
   - Plan: add a pharmacophore satisfaction component to the confidence
     score — verify HBD/HBA complementarity, hydrophobic burial, etc.
     between the consensus pose and the receptor binding site.
   - Requires the de novo mode pocket analysis as a prerequisite
     (need pocket residue features to compare against ligand features).

### Expanded Benchmarking
   - Current 8-target set is sufficient for development but too small
     for scientific claims.
   - Target: CASF-2016 (285 complexes, community standard for scoring
     power / ranking power / docking power evaluation).

### Engine Diversity Upgrade: DiffDock
   - Smina and Gnina are both forks of AutoDock Vina; they share the
     same iterated local search sampling algorithm and differ only at
     the rescoring stage. Two of three current engines therefore sample
     pose space identically, weakening the independence assumption that
     makes cross-engine spatial agreement meaningful.
   - Plan: replace Smina with DiffDock (end-to-end diffusion model,
     no explicit force field, trained on PDBbind). Target three-engine
     set: DiffDock (generative ML) + Gnina (physics sampling + CNN
     rescoring) + LeDock (simulated annealing). Three orthogonal
     search paradigms, genuinely independent failure modes.
   - EnsembleManager and EngineType are designed to accommodate this
     with minimal structural changes.
   - This upgrade directly supports the publication claim that spatial
     agreement across architecturally distinct engines is a stronger
     signal than score aggregation.

### Publication Pathway
   The core publishable claim is: density-based clustering of the full
   cross-engine pose pool with engine-count-primary ranking outperforms
   fixed-threshold clustering (MetaDOCK) and score aggregation (dockECR)
   on CASF-2016 docking power (SR₁ = % targets where top-ranked pose
   has RMSD < 2.0 Å). This is a specific, testable, controlled claim.
   Two supporting experiments make it defensible:
   1. HDBSCAN vs. fixed 2.0 Å cutoff on the same pose pool: same
      engines, same poses, different clustering — shows density-based
      adaptation improves SR₁.
   2. Engine-count ranking vs. best-score-in-cluster on the same
      clusters: same clustering, different selection — shows
      score-agnostic ranking is more robust across protein families.
   The de novo mode (fpocket + pharmacophore matching) is a second
   contribution, likely a separate paper or a major extension section.
   Both require CASF-2016 validation to make performance claims.

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

No test suite currently exists (post-v1.x goal — see Known Issues).
When implemented, tests should:
- Live in tests/ using pytest
- Mock all external binaries (fpocket, smina, gnina, ledock, obabel)
- Mock all RCSB network calls
- Use fixtures in tests/fixtures/
- Name files test_<module>.py
- Cover every public method and every error path

---

## External Dependencies

### Python (see requirements.txt and pyproject.toml)

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

Each binary download in Cell 1 routes through
utils.download_and_verify_binary(), which reads the ELF header after
wget and rejects anything whose magic bytes are not \x7fELF or whose
e_machine field is not EM_X86_64 (0x3E). Failures delete the file and
raise — there is no path by which an HTML error page or wrong-arch
binary can persist on disk past Cell 1. A corrupt file from a previous
session is detected and re-downloaded automatically. Use utils.py
check_environment() to additionally verify all binaries are on PATH.

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

### Current Benchmark (last recorded run, 2026-04-30)

- Success:    3/8 (1OWE 0.23 Å, 1STP 0.65 Å, 1FJS 1.63 Å)
- Poor:       3/8 (1ETT 3.70 Å, 1HXW 5.55 Å, 1A30 7.78 Å)
- Error:      2/8 (1HSG, 1IEP — clustering succeeded; RMSD validation
              fails at CalcRMS due to ref_mol construction failure;
              see Known Issues item 4)
- Run was at POSES_PER_ENGINE=20 and EXHAUSTIVENESS=8. Class C
  auto-retry confirmed: 1HSG/1HXW/1IEP all triggered the ideal-SDF
  fallback (100% crystal mol2 fail rate → 0-7% post-retry). Clustering
  is working for all 8 targets. The 2 remaining Errors are isolated to
  the RMSD validation step, not docking or consensus.

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
- Functional pytest suite in tests/ covers all modules with mocked
  external binaries and RCSB calls
