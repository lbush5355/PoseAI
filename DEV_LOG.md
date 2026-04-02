# PoseAI Development Log

## Date: April 2, 2026
**Current Phase:** Phase 1: Environment and Software Installation
**Status:** Complete

---

## Environment Configuration
To address architecture incompatibilities between legacy bioinformatics software and Apple Silicon (M2), a dedicated translation environment was established.

* **Environment Name:** `poseai_intel`
* **Architecture:** `osx-64` (Intel Translation via Rosetta 2)
* **Python Version:** 3.10
* **Package Manager:** Mamba/Conda

---

## Software Specifications and Algorithmic Variance
The project utilizes three distinct mathematical approaches to molecular docking to ensure a comprehensive consensus model.

| Software | Command | Algorithm Category | Search Method |
| :--- | :--- | :--- | :--- |
| **Smina** | `smina` | Empirical/Machine Learning | Gradient Descent |
| **AutoDock 4** | `autodock4` | Evolutionary | Lamarckian Genetic Algorithm |
| **LeDock** | `./bin/ledock` | Physics-based | Simulated Annealing |



---

## Technical Issues and Resolutions

### 1. Architecture Compatibility
* **Issue:** Several Bioconda packages, including `autodock`, do not have native ARM64 (M2) builds.
* **Resolution:** Configured the environment to `osx-64` mode using `CONDA_SUBDIR=osx-64`. This forces the installation of Intel-based packages which are then executed via Rosetta 2.

### 2. File Integrity Errors
* **Issue:** `curl` requests resulted in 404 error pages saved as 9-byte text files instead of the intended executable binaries.
* **Resolution:** Performed manual downloads via browser for `ledock` and `gnina` (v1.1) to verify file sizes exceeded 70MB and ensured executable integrity.

### 3. Operating System Security Permissions
* **Issue:** Unsigned scientific binaries were blocked by macOS security protocols for unidentified developers.
* **Resolution:** Applied `chmod +x` permissions and manually authorized binaries within System Settings > Privacy & Security.

---

## Planned Activities: Phase 2 (Data Preparation)
* Standardize Protein (PDB) preparation including solvent removal and addition of polar hydrogens.
* Convert Ligands (SDF) to PDBQT format and calculate Gasteiger charges.
* Define XYZ search coordinates for consistent docking across all three platforms.