# PoseAI: A Multi-Engine Consensus Molecular Docking Pipeline

PoseAI is a modular computational framework designed to execute and harmonize ligand-binding simulations across multiple docking scoring functions. By integrating **Smina** (empirical), **LeDock** (physics-based), and **AutoDock 4** (semi-empirical), the pipeline identifies high-confidence binding modes through unsupervised machine learning and symmetry-corrected spatial clustering.

The system is optimized for **Apple Silicon (M2)** architectures, utilizing a Rosetta 2 Intel-bridge to ensure compatibility with the specific architectural requirements of legacy bioinformatics binaries.

---

## **System Architecture**

The pipeline is composed of four primary functional modules:

* **Structural Preprocessor (`src/preprocessor.py`):** Automates PDB retrieval, solvent/ion stripping, and engine-specific formatting (`PDBQT`, `pro.pdb`).
* **Pocket Detector (`src/detector.py`):** Identifies binding site centroids using co-crystallized ligand coordinates via Biopython and NumPy.
* **Docking Manager (`src/docking.py`):** Orchestrates multi-threaded subprocess execution for Smina and LeDock with dynamic configuration generation.
* **Consensus Parser (`src/parser.py`):** Normalizes binding affinities and extracts Cartesian coordinates for downstream geometric analysis.

---

## **Project Development Roadmap**

### **Phase 1: Data Integration and Feature Engineering (Week 1)**
* **Standardization:** Develop a unified extraction layer for Gibbs Free Energy ($\Delta G$) and atomic Cartesian Coordinates ($x, y, z$) from heterogeneous trajectory files.
* **Normalization:** Implement Z-score statistical scaling to allow comparative analysis of disparate energy scales across engines.
* **Feature Construction:** Transform 3D molecular conformations into fixed-dimension tensors for clustering input.

### **Phase 2: Unsupervised Learning and Consensus Clustering (Week 2)**
* **Distance Metrics:** Integrate symmetry-corrected Root-Mean-Square Deviation (RMSD) algorithms to account for ligand graph symmetry.
* **Clustering Implementation:** Utilize Density-Based Spatial Clustering of Applications with Noise (DBSCAN) to identify spatially dense binding hotspots.
* **Confidence Scoring:** Assign a **Consensus Ratio** based on the percentage of unique engines contributing to each cluster.

### **Phase 3: Scientific Validation and Visualization (Week 3)**
* **Redocking Benchmarks:** Validate the pipeline using Human Abl Kinase (**PDB 1IEP**). Success is defined as a primary consensus cluster achieving an $RMSD \leq 2.0$ Å relative to the crystal structure.
* **Visualization Deployment:** Automate the generation of PyMOL (`.pml`) session files for 3D structural analysis of consensus results.

---

## **Validation Protocol**

The project utilizes a three-tier Quality Assurance (QA) system:

1.  **Structural Audit:** Automated checksums and file-size verification to ensure simulation output integrity.
2.  **Geometric Accuracy:** The **Redocking Test**—verifying the pipeline's ability to recover the native crystal structure pose of a known ligand.
3.  **Consensus Convergence:** A metric-based scoring system where **High Confidence** is defined as $\geq 66\%$ agreement among the integrated docking engines.

## **Installation and Setup**

### **1. Environment Configuration**
The project requires a specific `osx-64` architecture to maintain compatibility with Intel-based bioinformatics binaries on Apple Silicon. Users should utilize Mamba or Conda to create the environment from the provided specification:

**Create the environment:**
`mamba env create -f environment_intel.yml`

**Activate the specialized Intel-bridge environment:**
`mamba activate poseai_intel`

### **2. Binary Dependencies**
The following executables must be present in the `bin/` directory. Note that these binaries are architecture-specific and are excluded from version control:
* **lepro / ledock:** Required for LeDock-specific protein protonation and simulation.
* **smina:** Integrated via the Conda environment, but can be manually linked if necessary.

### **3. Execution**
To execute the standard validation pipeline and benchmark the system against the Human Abl Kinase (PDB 1IEP) target, run the central controller:

`python3 main.py`

### **4. Directory Structure**
For the pipeline to function correctly, the following structure must be maintained:
* **bin/**: Contains Mach-O 64-bit x86_64 executables.
* **src/**: Contains the Python source modules for the pipeline.
* **data/**: Automated output directory for PDB structures and docking trajectories.

---

## **References and Citations**

The PoseAI framework integrates several peer-reviewed docking engines and bioinformatics libraries. Please cite the following primary literature when utilizing this pipeline for research or analysis:

### **Docking Engines**

* **Smina:** Koes, D. R., Baumgartner, M. P., & Camacho, C. J. (2013). Lessons learned from optimizing docking scoring functions. *Journal of Chemical Information and Modeling*, 53(8), 1893–1904. https://doi.org/10.1021/ci300604z
* **LeDock:** Zhao, H., & Caflisch, A. (2013). Molecular docking by simulated annealing and minimization. *European Journal of Medicinal Chemistry*, 61, 155–172. https://doi.org/10.1016/j.ejmech.2013.01.057
* **AutoDock Vina (Smina Core):** Trott, O., & Olson, A. J. (2010). AutoDock Vina: Improving the speed and accuracy of docking with a new scoring function, efficient optimization, and multithreading. *Journal of Computational Chemistry*, 31(2), 455–461. https://doi.org/10.1002/jcc.21334

### **Software and Libraries**

* **Biopython:** Cock, P. J., et al. (2009). Biopython: Freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics*, 25(11), 1422–1423. https://doi.org/10.1093/bioinformatics/btp163
* **Open Babel:** O'Boyle, N. M., et al. (2011). Open Babel: An open chemical toolbox. *Journal of Cheminformatics*, 3(1), 33. https://doi.org/10.1186/1758-2946-3-33
* **NumPy:** Harris, C. R., et al. (2020). Array programming with NumPy. *Nature*, 585(7825), 357–362. https://doi.org/10.1038/s41586-020-2649-2

### **Structural Data Sources**

* **RCSB Protein Data Bank:** Berman, H. M., et al. (2000). The Protein Data Bank. *Nucleic Acids Research*, 28(1), 235–242. https://doi.org/10.1093/nar/28.1.235