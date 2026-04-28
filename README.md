# PoseAI: A Multi-Engine Consensus Molecular Docking Pipeline

PoseAI is a modular computational framework designed to execute and harmonize ligand-binding simulations across multiple docking scoring functions. By integrating **Smina** (empirical), **Gnina** (deep learning-based), and **LeDock** (physics-based), the pipeline identifies high-confidence binding modes through unsupervised machine learning and symmetry-corrected spatial clustering.

The system has been optimized for execution within **Google Colab** to guarantee reproducibility, abstracting away complex local dependency management and system architectures by automatically provisioning a standardized Ubuntu environment.

## System Architecture

The pipeline is composed of distinct functional modules located in the `src/` directory:

*   **`preprocessor.py`**: Automates structural retrieval, solvent/ion stripping, and engine-specific formatting (PDBQT, mol2) via Open Babel.
*   **`site_finder.py`**: Identifies binding site centroids using co-crystallized ligand coordinates or dynamic pocket detection.
*   **`docking.py`**: Orchestrates multi-threaded subprocess execution for Gnina, Smina, and LeDock with dynamic configuration generation.
*   **`consensus.py`**: Analyzes generated poses, calculates symmetry-corrected Root-Mean-Square Deviation (RMSD), executes spatial clustering, and generates consensus confidence scores.
*   **`visualizer.py`**: Generates interactive 3D native-overlay visualizations of the consensus clusters using py3Dmol.
*   **`config.py`**: Manages global pipeline variables and dynamic engine tuning.

## Project Development Roadmap & Methodology

**Phase 1: Data Integration and Feature Engineering**
*   **Standardization:** Unified extraction of binding affinities and atomic Cartesian Coordinates (x, y, z) from heterogeneous engine trajectory files.
*   **Structure Preparation:** Automated generation of required topology and partial charge formats for multi-engine compatibility.

**Phase 2: Unsupervised Learning and Consensus Clustering**
*   **Distance Metrics:** Integration of robust, symmetry-aware RMSD algorithms to accurately calculate spatial divergence between poses, accounting for ligand graph symmetry.
*   **Clustering Implementation:** Identification of spatially dense binding hotspots to evaluate cross-engine convergence.
*   **Confidence Scoring:** Assignment of a Consensus Confidence Score based on the percentage of unique engines contributing to the top cluster.

**Phase 3: Scientific Validation and Visualization**
*   **Redocking Benchmarks:** Validation of the pipeline using standard PDBbind/CASF targets (e.g., Human Abl Kinase, PDB 1IEP). Success is defined as a primary consensus cluster achieving an RMSD ≤ 2.0 Å relative to the crystal structure.
*   **Visualization Deployment:** Automated generation of 3D HTML reports for immediate structural analysis of consensus results.

## Reproducibility Instructions (For Evaluation & Grading)

To ensure complete reproducibility of the pipeline and batch validation results, please follow these steps:

### 1. Cloud Storage Setup
1. Log into your Google account and open Google Drive.
2. Create a root directory named `PoseAI` in your Drive (`MyDrive/PoseAI`).
3. Upload the `src/` directory (containing the Python modules) to `MyDrive/PoseAI/src/`.
4. Upload the `dataset/` directory to `MyDrive/PoseAI/dataset/`.

### 2. Expected Dataset Structure
For the batch validation to function properly, the `dataset/` directory must follow the standard PDBbind/CASF naming conventions. Each target requires its own subfolder containing the strictly named protein and ligand files:

```text
PoseAI/
├── dataset/
│   ├── 1hsg/
│   │   ├── 1hsg_protein.pdb
│   │   └── 1hsg_ligand.mol2
│   ├── 1iep/
│   │   ├── 1iep_protein.pdb
│   │   └── 1iep_ligand.mol2
```

### 3. Execution Environment
1. Open the `PoseAI.ipynb` Notebook in Google Colab.
2. Go to **Runtime > Change runtime type** and ensure you are using a standard Python 3 compute instance.
3. Execute the Notebook sequentially from **Step 1 to Step 5**.

### 4. Automated Dependency Management
The pipeline is designed to act as an autonomous controller. **Step 1** of the notebook will automatically:
*   Mount your Google Drive to access the `PoseAI` directory.
*   Synchronize local modules.
*   Install required Python libraries (`rdkit`, `py3Dmol`, `hdbscan`, etc.).
*   Add 32-bit architecture to the Colab Ubuntu instance (required for legacy LeDock binaries).
*   Download and configure Smina, Gnina, and LeDock executables directly into the Colab environment.
*   Generate an `environment.yml` file for local Conda environment replication.

### 5. Validation Protocol
Executing **Step 5: Local Dataset Batch Validation** will automatically benchmark the pipeline against all targets in the `dataset/` folder. It evaluates structural integrity, performs ensemble docking, clusters the results, and outputs a `batch_summary.csv` alongside 3D `.html` visualization files in a newly generated `batch_results/` directory.

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

**Structural Data Sources**
*   **PDBbind:** Liu, Z., et al. (2017). Forging the Basis for Developing Protein-Ligand Interaction Scoring Functions. *Accounts of Chemical Research*, 50(2): 302-309.
*   **RCSB Protein Data Bank:** Berman, H. M., et al. (2000). The Protein Data Bank. *Nucleic Acids Research*, 28(1), 235–242.
