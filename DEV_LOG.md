# PoseAI Development Log

## Date: April 2, 2026
**Current Phase:** Phase 2: Automated Structural Preparation and Pocket Detection
**Status:** Complete

---

### 1. Technical Accomplishments

#### 1.1 Structural Pre-processing Logic (`src/preprocessor.py`)
A standardized pipeline was established to transform raw Protein Data Bank (PDB) files into software-specific inputs required by the consensus docking engines.

* **Remote Data Retrieval:** Integrated the `urllib` library to programmatically fetch structural data from the RCSB PDB servers.
* **Chemical Refinement:** Implemented automated filtration of water molecules (HOH), ions, and non-protein crystallographic artifacts to ensure search space integrity.
* **Software-Specific Formatting:**
    * **AutoDock 4:** Automated PDBQT generation with Gasteiger charge assignment via the `prepare_receptor4.py` utility.
    * **LeDock:** Integrated the `lepro` binary to generate the required `pro.pdb` structural files and `dock.in` configuration templates.

#### 1.2 Automated Coordinate Detection (`src/detector.py`)
Developed a module to identify target binding site coordinates without manual intervention.

* **Ligand Centroid Extraction:** Utilized the Biopython `PDBParser` to identify co-crystallized hetero-atoms (HETATM).
* **Selection Algorithm:** Implemented an exclusion list to bypass common solvents (GOL, SO4, CL, DMS) and focus on bioactive small molecules (e.g., STI in PDB 1IEP).
* **Spatial Calculation:** Applied NumPy-based vector averaging to calculate the 3D centroid of the ligand, establishing the Cartesian coordinates for the search volume.

---

### 2. Technical Issues and Resolutions

#### 2.1 Shell Path Parsing (Exit Status 127)
* **Issue:** Subprocess executions failed when the project directory path contained spaces (e.g., `/Spring 2026/`).
* **Resolution:** Re-implemented all `subprocess.run` calls to include explicit double-quoting of all absolute and relative file paths to ensure correct shell interpretation.

#### 2.2 Binary Integrity (HTML Redirect Error)
* **Issue:** The `lepro` binary was initially identified as an HTML 301 Redirect file rather than a Mach-O executable, causing shell syntax errors.
* **Resolution:** Re-acquired the binary using `curl -L` to follow server-side redirects. Verified the **Mach-O 64-bit x86_64** architecture using the `file` utility to confirm compatibility with the Rosetta 2 translation layer.

#### 2.3 Syntax Normalization
* **Issue:** Encountered `SyntaxError` at f-string declarations caused by hidden non-breaking space characters (ASCII 160) introduced during code migration.
* **Resolution:** Overwrote all source files using terminal-native heredocs (`cat <<EOF`) to ensure standard UTF-8 encoding.

---

### 3. Validation Results: PDB 1IEP
The pipeline successfully processed the Human Abl Kinase (1IEP) benchmark system:

* **Protein Preparation:** Generated `1iep_clean.pdb`.
* **AutoDock Input:** Generated `1iep_clean.pdbqt`.
* **LeDock Input:** Generated `pro.pdb`.
* **Calculated Centroid:** `X: 14.081, Y: 75.123, Z: 37.365`.

---

### 4. Planned Activities: Phase 3
* Implementation of the `DockingManager` for multi-engine execution.
* Development of the configuration writer for LeDock’s boundary-based search box.
* Integration of a unified scoring parser for binding affinity normalization.