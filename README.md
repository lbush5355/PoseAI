# PoseAI: A Multi-Engine Consensus Molecular Docking Pipeline

**Status: v0.9 pre-release** — The pipeline is fully functional end-to-end. The v1.x milestone requires a minimum of 6/8 gold-standard targets scoring Success or Acceptable, zero Error results, and functional unit test coverage. Current benchmark: 3/8 Success, 3/8 Poor, 2/8 Error (see [Benchmark Performance](#current-benchmark-performance)).

PoseAI is a modular computational framework designed to execute and harmonize ligand-binding simulations across multiple docking scoring functions. By integrating **Gnina** (Vina-family sampler with CNN rescoring), **LeDock** (simulated annealing / independent architecture), and **Smina** (Vina-family empirical force-field baseline), the pipeline collects the full pose output from all engines into a shared pool and applies HDBSCAN clustering on a pairwise heavy-atom RMSD matrix to identify binding modes where independently operating engines converge in 3D space. Cross-engine agreement is quantified as an **Ensemble Confidence Score** combining multi-engine representation (70%) and cluster population (30%).

The system runs exclusively on **Google Colab** to guarantee reproducibility, abstracting away dependency management and system architecture differences by automatically provisioning a standardized Ubuntu environment with ELF-validated engine binaries.

Developed as a course project for CHEM 4640 at the University of Colorado Denver, Spring 2026.

---

## System Architecture

The pipeline is composed of distinct functional modules in the `src/` directory:

*   **`runtime.py`**: Colab session bootstrap. Installs Python dependencies, adds 32-bit architecture support (required for the LeDock/lepro legacy binaries), builds fpocket from source if absent, and downloads each docking engine binary via `download_and_verify_binary()`, which reads the ELF header and rejects HTML error pages and wrong-architecture files before they reach subprocess.
*   **`preprocessor.py`**: Automates structural retrieval from RCSB, solvent/ion stripping, and engine-specific format generation (PDBQT, mol2) via Open Babel. Provides `get_ligand_centroid()` (active binding site detection via crystal ligand heavy-atom coordinates), `fetch_rcsb_smiles()`, `fetch_rcsb_ideal_sdf()`, and `fetch_rcsb_entry_ligand_codes()` (GraphQL lookup of drug-like non-polymer entities for a given PDB entry).
*   **`docking.py`**: Orchestrates parallel subprocess execution of Gnina, Smina, and LeDock via `ProcessPoolExecutor`, with adaptive CPU allocation based on GPU availability and independent per-engine and pool-level timeout enforcement.
*   **`consensus.py`**: Loads and topology-standardizes poses across all engines using a 5-level master template fallback (crystal mol2/sdf/pdb → RCSB ideal SDF → SMILES → stereo-stripped SMILES retry → 3D inference from docking output). Computes a pairwise symmetry-corrected heavy-atom RMSD matrix (parallelized above a configurable pose-count threshold), clusters via HDBSCAN, and ranks clusters deterministically by engine count, cluster size, and intra-cluster RMSD.
*   **`pipeline.py`**: Unified end-to-end orchestration. `run_from_rcsb()` handles single-target runs from RCSB; `run_from_local()` handles batch validation over PDBbind-layout datasets. Both return a `TargetResult` dataclass carrying status, native RMSD, confidence score, cluster DataFrame, analyzer state, and docking results.
*   **`run_history.py`**: Persistent run tracking. `append_run()` writes one row per target to `run_history.csv` on Google Drive; `per_target_stats()` aggregates cross-session success rate, pass rate, and best RMSD per target for presentation and evaluation.
*   **`visualizer.py`**: Generates interactive 3D native-overlay visualizations using py3Dmol, with engine-color-coded poses (Smina: green, Gnina: blue, LeDock: red, consensus cluster: gold). Exports HTML reports for review outside Colab.
*   **`config.py`**: Centralized typed dataclass configuration with named fields for all tunable parameters. Supports runtime overrides via `POSEAI_<MODULE>__<PARAM>` environment variables without modifying source. See [Configuration](#configuration) for details.

---

## Pipeline Implementation

### Phase 1: Data Integration and Structure Preparation

*   **Multi-format ingestion**: Unified pose extraction from Gnina/Smina SDF outputs and LeDock `.dok` files. LeDock's `.dok` format concatenates poses as PDB-like blocks separated by `REMARK Cluster` markers; the parser splits at those boundaries and converts each block independently via Open Babel, recovering all poses that a naive single-pass `obabel -ipdb` call silently drops.
*   **Structure preparation**: Automated RCSB fetch, water/ion removal, ligand isolation by 1–3 character alphanumeric residue code (per the PDB Chemical Component Dictionary spec), and PDBQT/mol2 generation with Gasteiger charges. When no ligand code is supplied, the code auto-resolves it from RCSB's GraphQL non-polymer entity endpoint, filtered to drug-like candidates.
*   **Binding site detection**: Active method is `get_ligand_centroid()` in `preprocessor.py`, which computes the docking search box center from the co-crystallized ligand's heavy-atom coordinates. Raises `ValueError` on parse failure to prevent silent docking-against-wrong-region errors that would otherwise appear as Poor RMSD without a diagnostic.

### Phase 2: Unsupervised Learning and Consensus Clustering

*   **Topology standardization**: All engine outputs are reconciled to a canonical bond-order graph before RMSD computation to eliminate atom-count and bond-perception mismatches. The master template is resolved by priority: (1) crystal mol2/sdf/pdb at the reference path, (2) RCSB ideal SDF, (3) user/RCSB SMILES, (4) stereo-stripped SMILES retry, (5) SMILES inferred from a successful docking output. If the crystal mol2 causes more than 50% bond-order assignment failures — characteristic of peptidomimetic ligands where depositor and engine bond perception conflict — standardization is automatically retried against the RCSB ideal SDF in memory without re-docking (~2 seconds).
*   **RMSD computation**: Pairwise symmetry-corrected heavy-atom RMSD via `rdMolAlign.GetBestRMS`. Parallelized across CPU cores when pose count exceeds a configurable threshold (default: 30 poses); sequential below it where multiprocessing overhead exceeds the savings.
*   **Clustering**: HDBSCAN with a precomputed distance matrix and `cluster_selection_epsilon` equal to the RMSD threshold (default: 2.0 Å). Clusters are ranked deterministically: engine count descending, cluster size descending, intra-cluster RMSD ascending. The three-key sort prevents non-deterministic best-cluster selection when two clusters tie on size.
*   **Confidence scoring**: Ensemble Confidence Score = 0.70 × (engines in best cluster / 3) + 0.30 × (cluster size / ideal size), bounded to [0, 1]. Weights, engine denominator, and ideal size are all configurable.

### Phase 3: Scientific Validation and Visualization

*   **Redocking benchmarks**: Validated against the PDBbind gold-standard set (see [Benchmark Performance](#current-benchmark-performance)). Native RMSD is strict in-place heavy-atom RMSD between the top consensus pose and the crystal ligand — no alignment applied. The reference molecule is constructed by combining crystal coordinates with the canonical bond graph from whichever master template was selected, so the topology source is always auditable in the run log.
*   **Cross-session tracking**: `run_history.csv` accumulates one row per run across Colab sessions. `per_target_stats()` surfaces per-target success rate, pass rate, and best RMSD at the end of each batch run without requiring log inspection.
*   **Visualization**: HTML reports with engine-color-coded poses and consensus cluster highlight are generated for every target and persisted to `batch_results/` on Drive.

---

## Future Roadmap (Post-v1.x)

The following development directions are aimed at advancing PoseAI from a validated redocking pipeline toward a tool capable of genuine scientific contribution to structure-based drug discovery. Most existing consensus docking pipelines assume a known binding site and evaluate poses purely by geometric convergence; the priorities below address both of those limitations.

### Scientific Capabilities

**De novo binding site discovery.** The current pipeline requires a co-crystallized ligand to define the docking search box via `get_ligand_centroid()`, restricting it to redocking against targets with known structures. The planned extension uses fpocket — already built and installed in Cell 1 — to identify candidate binding pockets from receptor surface geometry alone via Voronoi tessellation, then selects the most compatible pocket by matching the ligand's pharmacophoric features (hydrogen bond donors/acceptors, hydrophobic regions, aromatic rings via RDKit's `MolChemicalFeatures`) against each pocket's chemical descriptors. `PocketAnalyzer` in `site_finder.py` already implements the fpocket wrapper and output parser; the pharmacophore-to-pocket matching layer is the remaining work. This enables docking against novel receptor targets with no crystallographic data — the primary scientific limitation of the current approach.

**Pharmacophore-aware confidence scoring.** The current Ensemble Confidence Score measures only spatial convergence — how tightly poses cluster and how many engines agree. It cannot distinguish a geometrically convergent pose that satisfies key chemical interactions from one that does not, a known failure mode in consensus docking. Adding a pharmacophore satisfaction component — verifying that HBD groups face HBA residues, hydrophobic groups are buried in hydrophobic pockets, etc. — would make the confidence score chemically meaningful rather than purely geometric.

**Expanded benchmarking against CASF-2016.** The current 8-target validation set is appropriate for development but too small to make performance claims against the field. The CASF-2016 benchmark (285 protein-ligand complexes from PDBbind with standardized scoring, ranking, and docking power evaluations) is the community standard for evaluating docking pipelines. Validation against CASF-2016 would establish PoseAI's scientific standing and make its performance directly comparable to published methods.

**Engine diversity upgrade: DiffDock integration.** Smina and Gnina share the same underlying pose sampling algorithm — both are forks of AutoDock Vina using iterated local search — and differ only at the rescoring stage. This means two of three current engines sample pose space the same way, limiting the architectural independence that makes cross-engine spatial agreement meaningful. The planned upgrade replaces Smina with **DiffDock**, an end-to-end diffusion model trained on PDBbind that treats binding pose generation as a generative modeling problem with no explicit force field or hand-engineered scoring function. The resulting three-engine set — DiffDock (generative ML), Gnina (physics sampling + learned rescoring), LeDock (stochastic physics) — represents three genuinely orthogonal search paradigms. The core hypothesis is that engines failing by different mechanisms will fail on different targets, making their spatial agreement a stronger and more discriminating signal than agreement between Vina-family methods. The `EnsembleManager` and `EngineType` enum are designed to accommodate this with minimal structural changes.

**Multi-residue ligand support.** Targets where the inhibitor is deposited as a polymer chain (e.g., 1A30) return no results from RCSB's non-polymer entity endpoint, preventing ideal SDF retrieval. Full support requires composing the ligand topology from its constituent PDB Chemical Component Dictionary residues.

**Receptor flexibility.** All current docking uses a rigid receptor model. For targets with known conformational variability — kinases, GPCRs, allosteric sites — docking against an ensemble of receptor conformations generated from MD snapshots or rotamer sampling would improve accuracy and is a natural extension of the existing ensemble philosophy.

### Technical Fixes

**RMSD recovery for topology-incompatible ligands.** For 1HSG and 1IEP, `AssignBondOrdersFromTemplate` fails in both directions, blocking native RMSD calculation even though clustering succeeds. The planned fix is an MCS-based coordinate-copy overlay: match the crystal mol2 atom positions onto the ideal SDF graph via maximum common substructure, then compute RMSD against the reconciled reference.

### Engineering

**Functional test coverage.** No test suite currently exists. Priority targets: `_summarize_clusters()` deterministic sort, `_load_dok()` multi-pose recovery, `analyze_ensemble()` Class C retry trigger, and `_grade_rmsd()` boundary conditions — all mockable without running docking engines.

**Pipeline checkpointing.** A crash at any stage discards all work for that target. Saving intermediate results after preprocessing, docking, and clustering would allow restarts without full re-runs.

**`max_engines_for_consensus` config parameter.** The confidence score denominator is hardcoded to 3.0 in `get_confidence_score()`. Adding a fourth engine without updating source silently caps its contribution.

### Path to Publication

PoseAI makes two specific methodological claims that are testable against published literature on a community-standard benchmark.

**Core claim.** Score-agnostic density-based consensus — HDBSCAN on a full all-pairs heavy-atom RMSD matrix with engine-count-primary cluster ranking — identifies near-native binding modes more reliably than fixed-threshold clustering with score-based selection. The existing implementation embodies this claim on 8 targets; the question is whether it holds at the scale needed for a peer-reviewed result.

**Experiment 1: clustering method.** Run the same engine pool (Gnina + LeDock + Smina) on CASF-2016 (285 protein-ligand complexes) under two conditions: HDBSCAN with adaptive density detection, and a fixed 2.0 Å RMSD cutoff applied to the same pose pool. Success rate and mean RMSD across all 285 targets isolates the contribution of the clustering algorithm from all other pipeline choices.

**Experiment 2: ranking criterion.** On the same CASF-2016 runs, compare engine-count-primary cluster selection against best-score-in-cluster selection (using Gnina's CNN score as the representative). A statistically significant improvement by engine-count-primary ranking supports the score-agnostic independence argument independently of Experiment 1.

**Comparison baselines.** MetaDOCK (fixed 2.5 Å threshold, best-scored pose in largest cluster) and dockECR (exponential consensus ranking with RMSD component) are re-implemented or rerun on the same CASF-2016 set. Both are open methodologies with published code, making a same-dataset comparison feasible.

**Engine upgrade prerequisite.** The DiffDock engine swap (Smina → DiffDock) is a prerequisite for the strongest version of the core claim, since the current Gnina/Smina pair shares Vina-family sampling. With DiffDock replacing Smina, the three engines represent genuinely orthogonal paradigms, and cross-engine agreement becomes a stronger independence signal. Experiments 1 and 2 can be run on the current three-engine set as a baseline, with the DiffDock results as the headline result.

**De novo mode as a second contribution.** The pharmacophore-to-pocket matching extension — fpocket for pocket detection, RDKit `MolChemicalFeatures` for ligand pharmacophore profiling — is architecturally independent of the consensus clustering work and addresses a different problem (no crystallographic reference available). If the de novo binding site experiments produce compelling results, this warrants a separate follow-up contribution rather than inclusion in the primary consensus validation study.

**What is needed before submission.** Access to the full CASF-2016 dataset (available from PDBbind with academic registration), DiffDock integration and validation on the current benchmark set, calibration of the Ensemble Confidence Score against a held-out subset, and functional unit test coverage to support the claims in the methods section.

---

## Key Engineering Challenges

The following problems required non-trivial root-cause diagnosis and drove significant design decisions in the current implementation.

**LeDock multi-pose loss.** LeDock writes all docked poses into a single `.dok` file as concatenated PDB-like blocks separated by `REMARK Cluster N` headers. A naive `obabel -ipdb` call treats the file as a single-model PDB and silently returns only the first pose, making LeDock appear to produce one result per run regardless of `n_poses`. The fix splits the file at `REMARK Cluster` boundaries and converts each block independently, recovering the full pose population that consensus clustering depends on.

**Topology standardization failure for peptidomimetic ligands (Class C retry).** For targets with flexible or peptidomimetic inhibitors (1HSG, 1IEP, 1HXW), Open Babel's perception of the crystal mol2 bond orders conflicts with how docking engines perceive the same ligand. This caused `AssignBondOrdersFromTemplate` to fail for 100% of poses, leaving the entire ensemble unstandardized. The resolution was a two-stage automatic retry: if the crystal mol2 master template causes a bond-order failure rate above 50%, the pipeline swaps the template to the RCSB ideal SDF and re-runs standardization in memory in roughly 2 seconds without re-docking. Post-retry failure rates dropped from 100% to 0–7% for all three affected targets.

**Non-deterministic cluster selection** (1OWE: 0.23 Å vs. 8.86 Å across runs). When two clusters had equal pose counts, sorting by size alone left the tiebreaker undefined, causing the correct near-native cluster and a decoy cluster to alternate as the top-ranked result between runs. Fixed by a deterministic three-key sort — engine count descending, cluster size descending, intra-cluster RMSD ascending — which is chemically motivated without peeking at the crystal structure.

**Binary validation: HTML error pages as executables.** Downloaded engine binaries were passed directly to `subprocess.Popen` without format verification. When a download URL returned an HTML error page or a partial file, the resulting "binary" crashed at process launch with `[Errno 8] Exec format error` — a symptom that gave no indication of the actual cause. The fix reads the ELF magic bytes (`\x7fELF`) and `e_machine` field (`0x3E` for x86-64) immediately after each download and deletes the file and raises on any mismatch, so a bad download is caught before it can silently poison a session.

**RCSB SMILES fetch returning None for all ligands.** The RCSB Chemical Component Dictionary API returns a JSON dict with uppercase keys (`SMILES_stereo`, `SMILES`). The original fetch code queried lowercase key names, so the lookup always returned `None` even for well-characterized ligands like imatinib (STI), forcing the pipeline to fall through to 3D inference as the master template source. Once identified, the fix was a one-line key correction, but diagnosing it required tracing the template resolution fallback chain to its origin.

---

## Current Benchmark Performance

**Benchmark run: 2026-04-30** | `POSES_PER_ENGINE=20`, `EXHAUSTIVENESS=8`

| Target | Status  | Native RMSD | Notes |
|--------|---------|-------------|-------|
| 1OWE   | Success | 0.23 Å      | |
| 1STP   | Success | 0.65 Å      | |
| 1FJS   | Success | 1.63 Å      | |
| 1ETT   | Poor    | 3.70 Å      | |
| 1HXW   | Poor    | 5.55 Å      | Class C retry triggered; clustering succeeded |
| 1A30   | Poor    | 7.78 Å      | Multi-residue ligand; no ideal SDF via GraphQL |
| 1HSG   | Error   | —           | Clustering succeeded; RMSD blocked by topology mismatch |
| 1IEP   | Error   | —           | Clustering succeeded; RMSD blocked by topology mismatch |

**Summary**: 3/8 Success, 3/8 Poor, 2/8 Error. Class C auto-retry (crystal mol2 → ideal SDF fallback) triggered for 1HSG, 1HXW, and 1IEP; post-retry bond-order failure rate dropped from 100% to 0–7% for all three. The 2 Error results are isolated to the RMSD validation step — docking and consensus clustering succeeded for all 8 targets.

**v1.x release criteria**: Minimum 6/8 targets scoring Success or Acceptable, zero NaN/Error results, all critical known issues resolved, and functional test coverage in `tests/`.

### Grading Scale

| RMSD        | Grade      | Interpretation                      |
|-------------|------------|-------------------------------------|
| < 2.0 Å     | Success    | Near-native pose recovered          |
| 2.0–3.0 Å   | Acceptable | Plausible binding mode              |
| > 3.0 Å     | Poor       | Pipeline failed on this target      |
| — (None)    | Error      | RMSD calculation failed (see logs)  |

---

## Known Limitations

### Multi-Residue Ligands
Targets where the inhibitor is deposited as a polymer chain rather than a non-polymer entity (e.g., 1A30/TIQ) return no entries from RCSB's GraphQL non-polymer entity endpoint. These targets fall back to the crystal mol2 master template and are subject to the same depositor-vs-engine bond-perception conflict the ideal SDF path was designed to mitigate. Out of scope for v1.x.

### RMSD Validation for Complex Peptidomimetic Ligands
For 1HSG (MK1) and 1IEP (STI), the crystal mol2 and RCSB ideal SDF have incompatible bond perception in both the crystal-to-engine and engine-to-crystal assignment directions. `AssignBondOrdersFromTemplate` fails for all poses, and the fallback reference (ideal SDF at RCSB canonical coordinates rather than the crystal binding site) makes the resulting RMSD scientifically meaningless even when `CalcRMS` can run. The pipeline reports `Status=Error` with the cluster DataFrame and analyzer preserved for inspection. The generated HTML overlay is the recommended manual validation path for these targets until the MCS-based coordinate-copy fallback is implemented.

### No Functional Unit Tests
`tests/test_pipeline.py` currently validates imports and module structure only. All pipeline correctness is verified through the Colab batch run, and no automated regression protection exists for individual modules.

### Confidence Score Engine Denominator
`get_confidence_score()` normalizes engine count against a hardcoded value of 3.0. Adding a fourth docking engine to the ensemble without updating `config.py` will silently cap its contribution to the confidence score at the same level as a three-engine result.

---

## Configuration

All tunable parameters are defined as typed dataclass fields in `src/config.py` with validated defaults. The global singleton is lazy-initialized on first `get_config()` call and can be overridden at runtime without modifying source:

```bash
# POSEAI_<MODULE>__<PARAM>=<value>
export POSEAI_CONSENSUS__RMSD_THRESHOLD=2.5
export POSEAI_CONSENSUS__BOND_ORDER_FALLBACK_THRESHOLD=0.3
export POSEAI_DOCKING__DEFAULT_ENGINE_TIMEOUT=7200
```

Key parameters for batch validation:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `consensus.rmsd_threshold` | 2.0 Å | HDBSCAN epsilon and cluster-membership cutoff |
| `consensus.bond_order_fallback_threshold` | 0.5 | Fail-rate threshold triggering Class C ideal-SDF retry |
| `consensus.consensus_weight` | 0.7 | Multi-engine agreement weight in confidence score |
| `consensus.cluster_size_weight` | 0.3 | Cluster population weight in confidence score |
| `consensus.rmsd_parallel_threshold` | 30 | Pose count above which RMSD matrix is parallelized |
| `docking.default_engine_timeout` | 3600 s | Per-engine wall-clock timeout before process-group kill |
| `preprocessor.rcsb_timeout_s` | 30 s | HTTP timeout for RCSB structure and SDF downloads |

`EXHAUSTIVENESS` and `POSES_PER_ENGINE` are set as notebook variables in Cell 5 and passed through `PipelineParams` rather than `PoseAIConfig`.

---

## Reproducibility Instructions (For Evaluation & Grading)

To ensure complete reproducibility of the pipeline and batch validation results, follow these steps in order.

### 1. Runtime Environment

PoseAI is designed exclusively for **Google Colab** with a **GPU runtime**. An NVIDIA A100 is recommended to reproduce published benchmark times; a T4 will work but Gnina's CNN-scoring path will be slower. The docking engine binaries are Linux x86-64 ELF files and cannot run on macOS or Windows. All Python source modules are cloned automatically from GitHub at session start — no manual file uploads are required.

1. Open `PoseAI.ipynb` in Google Colab.
2. Go to **Runtime > Change runtime type** and select a **GPU** hardware accelerator (required for Gnina's CNN-scoring path).
3. Execute cells sequentially from **Step 1 to Step 5**.

### 2. Cloud Storage Setup (dataset only)

Source code is cloned from GitHub automatically. Only the dataset must be placed on Drive:

1. Log into your Google account and open Google Drive.
2. Create a root directory named `PoseAI` at `MyDrive/PoseAI`.
3. Upload the `dataset/` directory to `MyDrive/PoseAI/dataset/`.

> **Note:** Do **not** upload `src/` to Drive. Step 1 clones the latest `src/` modules directly from GitHub (`lbush5355/PoseAI`, branch `main`) into the Colab session on every startup.

### 3. Expected Dataset Structure

For batch validation to function, the `dataset/` directory must follow standard PDBbind naming conventions. Each target requires a subfolder with a protein PDB and a ligand mol2. PDBQT files are generated automatically by the pipeline if absent.

```text
MyDrive/PoseAI/
└── dataset/
    ├── 1ett/
    │   ├── 1ett_protein.pdb      (required)
    │   └── 1ett_ligand.mol2      (required)
    ├── 1owe/
    │   ├── 1owe_protein.pdb
    │   └── 1owe_ligand.mol2
    └── ...
```

The full set of gold-standard targets is: `1ett`, `1owe`, `1a30`, `1stp`, `1hxw`, `1iep`, `1fjs`, `1hsg`.

### 4. Automated Dependency Management

**Step 1** of the notebook automatically:

*   Mounts your Google Drive to access the `PoseAI` directory.
*   Clones the `PoseAI` repository from GitHub and copies `src/*.py` into the Colab session.
*   Installs required Python libraries (`rdkit`, `py3Dmol`, `hdbscan`, etc.).
*   Adds 32-bit (`i386`) architecture support to the Colab Ubuntu instance — required because the LeDock and lepro binaries are 32-bit legacy ELFs.
*   Downloads Smina, Gnina, LeDock, lepro, and fpocket; validates each download as a well-formed x86-64 ELF before allowing it to be invoked.

### 5. Validation Protocol

**Step 5: Local Dataset Batch Validation** benchmarks the pipeline against all targets in `dataset/`. For each target it runs ensemble docking, consensus clustering, and RMSD grading, then outputs:

*   `batch_summary.csv` — per-target status and RMSD for this run
*   `run_history.csv` — cumulative per-target statistics across all sessions (success rate, best RMSD)
*   3D `.html` visualization files in `batch_results/<TARGET>/`

Per-target aggregate statistics (success rate, pass rate, best RMSD) are also displayed inline at the end of Step 5 without leaving the notebook.

---

## Related Work

PoseAI was developed independently. The tools described below were identified post-hoc for methodological comparison to situate the work within the broader consensus docking literature; none of them informed the design or implementation of PoseAI.

Consensus docking approaches fall into two paradigms in the literature.

**Score and rank aggregation** is the dominant paradigm. Tools including dockECR (Gimeno et al., 2021), DockingPie (Paiardi et al., 2022), and DockM8 run each engine independently and combine per-molecule scores or ranks using exponential consensus ranking, Z-score normalization, or majority vote. These methods never compare poses spatially; an engine's contribution is its numerical score, not the geometry of its predicted binding mode.

**Pose-based spatial consensus** methods compare actual 3D pose coordinates across engines. dockECR includes an RMSD-Based Scoring component that computes pairwise RMSD between each engine's single best pose, using spatial agreement as a secondary confidence signal. **MetaDOCK** (Ramírez & Caballero, 2023) is the most direct precedent for PoseAI's approach: it pools the top-5 poses per engine into a shared set (15 total), applies fixed 2.5 Å RMSD-threshold grouping to the joint pool, and selects the best-scored pose from the largest cluster. **VoteDock** (Plewczynski et al., 2011) similarly pooled poses from seven engines and applied hierarchical clustering, predating density-based methods.

PoseAI extends this direction in two specific respects. First, it applies HDBSCAN to the full all-pairs heavy-atom RMSD matrix across the entire cross-engine pose pool — a density-based algorithm that adapts to the natural cluster structure of the pose distribution without requiring a predetermined distance cutoff or cluster count. MDSCAN (Ferruz et al., 2022) is the only prior work applying HDBSCAN to an RMSD distance matrix in structural biology, in the context of molecular dynamics trajectory clustering rather than docking. Second, PoseAI ranks clusters by the **number of contributing engines** as the primary selection criterion rather than by pose score. This makes the selection explicitly score-agnostic: the consensus binding mode is defined by where independently operating engines converge in 3D space, not by what any individual scoring function assigns.

On engine diversity, **ESSENCE-Dock** (Sánchez-Murcia et al., 2024) makes the strongest published argument for combining algorithmically distinct engines, pairing DiffDock (end-to-end diffusion model), Gnina (CNN-augmented Vina sampler), and LeadFinder (genetic algorithm). The current PoseAI engine set — Gnina and Smina (both Vina-family) plus LeDock (simulated annealing) — provides partial architectural diversity. Replacing Smina with DiffDock is a targeted roadmap item to achieve three fully orthogonal search paradigms and strengthen the core independence assumption that underlies the spatial consensus approach.

---

## Development Notes

The pipeline runs exclusively on **Google Colab** (Linux x86-64). Local setup via `pip install -r requirements.txt` supports editing and linting only — functional testing requires a Colab session with GPU runtime.

**Running tests locally:**
```bash
pytest tests/ -v
```
Note: the test suite currently covers imports and structure only. All functional validation runs through Colab Cell 5.

**Commit conventions:** `feat:` / `fix:` / `docs:` / `test:` / `refactor:` prefixes. Stage specific files rather than `git add -A`.

**Non-negotiable coding rules for all `src/` changes:**
- No bare `except` — always catch specific exceptions with a logged diagnostic
- No magic numbers — all constants in `config.py` as named dataclass fields
- No `os.chdir()` — use `cwd=` in subprocess calls
- No hardcoded `/content/` paths — all paths from constructor arguments or config
- No duplicate imports — one import per dependency at the top of the file
- No silent exception swallowing — every caught exception must be logged

Type hints are required on all function signatures. f-strings only — no `%` or `.format()`. See `CLAUDE.md` for the full coding standards reference.

---

## References and Citations

The PoseAI framework integrates several peer-reviewed docking engines and bioinformatics libraries. Please cite the following primary literature when utilizing this pipeline for research or analysis:

**Docking Engines**
*   **Smina:** Koes, D. R., Baumgartner, M. P., & Camacho, C. J. (2013). Lessons learned from optimizing docking scoring functions. *Journal of Chemical Information and Modeling*, 53(8), 1893–1904. https://doi.org/10.1021/ci300604z
*   **LeDock:** Zhao, H., & Caflisch, A. (2013). Molecular docking by simulated annealing and minimization. *European Journal of Medicinal Chemistry*, 61, 155–172. https://doi.org/10.1016/j.ejmech.2013.01.057
*   **Gnina:** McNutt, A., Li, Y., Meli, R., Aggarwal, R., Koes, D. R. (2025). GNINA 1.3: the next increment in molecular docking with deep learning. *Journal of Cheminformatics*. https://pubmed.ncbi.nlm.nih.gov/39837943/

**Software and Libraries**
*   **Biopython:** Cock, P. J., et al. (2009). Biopython: Freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics*, 25(11), 1422–1423. https://doi.org/10.1093/bioinformatics/btp163
*   **Open Babel:** O'Boyle, N. M., et al. (2011). Open Babel: An open chemical toolbox. *Journal of Cheminformatics*, 3(1), 33. https://doi.org/10.1186/1758-2946-3-33
*   **RDKit:** RDKit: Open-source cheminformatics. https://www.rdkit.org
*   **HDBSCAN:** McInnes, L., Healy, J., & Astels, S. (2017). hdbscan: Hierarchical density based clustering. *Journal of Open Source Software*, 2(11), 205. https://doi.org/10.21105/joss.00205
*   **py3Dmol:** Rego, N., & Koes, D. R. (2015). 3Dmol.js: Molecular visualization with WebGL. *Bioinformatics*, 31(8), 1322–1324. https://doi.org/10.1093/bioinformatics/btu829
*   **fpocket:** Le Guilloux, V., Schmidtke, P., & Tuffery, P. (2009). Fpocket: An open source platform for ligand pocket detection. *BMC Bioinformatics*, 10, 168. https://doi.org/10.1186/1471-2105-10-168
*   **NumPy:** Harris, C. R., et al. (2020). Array programming with NumPy. *Nature*, 585, 357–362. https://doi.org/10.1038/s41586-020-2649-2
*   **pandas:** The Pandas Development Team (2020). pandas-dev/pandas: Pandas. Zenodo. https://doi.org/10.5281/zenodo.3509134

**Structural Data Sources**
*   **PDBbind:** Liu, Z., et al. (2017). Forging the Basis for Developing Protein-Ligand Interaction Scoring Functions. *Accounts of Chemical Research*, 50(2): 302-309.
*   **RCSB Protein Data Bank:** Berman, H. M., et al. (2000). The Protein Data Bank. *Nucleic Acids Research*, 28(1), 235–242.
*   **RCSB PDB API:** Bittrich, S., et al. (2023). RCSB Protein Data Bank: Powerful new tools for exploring 3D structures of biological macromolecules for basic and applied research and education in fundamental biology, biomedicine, biotechnology, bioengineering and energy sciences. *Nucleic Acids Research*, 51(D1), D488–D501. https://doi.org/10.1093/nar/gkac1019
*   **CASF-2016:** Su, M., et al. (2019). Comparative Assessment of Scoring Functions: The CASF-2016 Update. *Journal of Chemical Information and Modeling*, 59(2), 895–913. https://doi.org/10.1021/acs.jcim.8b00545

**Related Consensus Docking Methods**
*   **dockECR:** Gimeno, A., et al. (2021). Open consensus docking and ranking protocol for virtual screening of small molecules. *Journal of Molecular Structure*, 1229, 129519. https://doi.org/10.1016/j.molstruc.2020.129519
*   **DockingPie:** Paiardi, G., et al. (2022). DockingPie: a consensus docking plugin for PyMOL. *Bioinformatics*, 38(17), 4233–4234. https://doi.org/10.1093/bioinformatics/btac452
*   **MetaDOCK:** Ramírez, D., & Caballero, J. (2023). MetaDOCK: A Combinatorial Molecular Docking Approach. *ACS Omega*, 8(6), 5718–5731. https://doi.org/10.1021/acsomega.2c07784
*   **ESSENCE-Dock:** Sánchez-Murcia, P. A., et al. (2024). ESSENCE-Dock: A Consensus-Based Approach to Enhance Virtual Screening Enrichment in Drug Discovery. *Journal of Chemical Information and Modeling*, 64(6), 1829–1843. https://doi.org/10.1021/acs.jcim.3c01617
*   **MDSCAN:** Ferruz, N., et al. (2022). MDSCAN: RMSD-based HDBSCAN clustering of long molecular dynamics. *Bioinformatics*, 38(23), 5191–5192. https://doi.org/10.1093/bioinformatics/btac666
*   **DiffDock:** Corso, G., et al. (2023). DiffDock: Diffusion Steps, Twists, and Turns for Molecular Docking. *International Conference on Learning Representations (ICLR)*. arXiv:2210.01776.
