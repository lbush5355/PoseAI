import os
import numpy as np
from Bio.PDB import PDBParser

class PocketDetector:
    def __init__(self, pdb_id, work_dir="data"):
        self.pdb_id = pdb_id.lower()
        self.raw_pdb = os.path.join(work_dir, self.pdb_id, f"{self.pdb_id}.pdb")
        self.parser = PDBParser(QUIET=True)

    def detect_from_ligand(self):
        """Extracts the centroid of the co-crystallized ligand (STI)."""
        print(f"Action: Detecting ligand centroid in {self.pdb_id.upper()}...")
        
        if not os.path.exists(self.raw_pdb):
            print("Error: Raw PDB not found for detection.")
            return None

        structure = self.parser.get_structure(self.pdb_id, self.raw_pdb)
        coords = []
        
        # Common non-drug ligands to ignore
        ignore_list = ["HOH", "CL", "NA", "SO4", "PO4", "EDO", "DMS", "GOL"]

        for model in structure:
            for chain in model:
                for residue in chain:
                    # H_ indicates a hetero-atom (ligand)
                    if residue.id[0].startswith("H_") and residue.resname not in ignore_list:
                        print(f"Status: Found potential ligand: {residue.resname}")
                        for atom in residue:
                            coords.append(atom.get_coord())
        
        if coords:
            centroid = np.mean(coords, axis=0)
            # Round to 3 decimal places for docking config files
            return np.round(centroid, 3)
        
        print("Status: No suitable co-crystallized ligand found.")
        return None