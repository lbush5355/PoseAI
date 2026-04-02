import os
import subprocess
import urllib.request

class PoseAIPreprocessor:
    def __init__(self, pdb_id, work_dir="data"):
        self.pdb_id = pdb_id.lower()
        self.output_path = os.path.join(work_dir, self.pdb_id)
        self.raw_pdb = os.path.join(self.output_path, f"{self.pdb_id}.pdb")
        self.clean_pdb = os.path.join(self.output_path, f"{self.pdb_id}_clean.pdb")
        self.pdbqt = f"{self.clean_pdb}qt"
        self.pro_pdb = os.path.join(self.output_path, "pro.pdb")
        
        if not os.path.exists(self.output_path):
            os.makedirs(self.output_path)

    def fetch_pdb(self):
        if not os.path.exists(self.raw_pdb):
            url = f"https://files.rcsb.org/download/{self.pdb_id.upper()}.pdb"
            urllib.request.urlretrieve(url, self.raw_pdb)

    def clean_structure(self):
        with open(self.raw_pdb, "r") as f_in, open(self.clean_pdb, "w") as f_out:
            for line in f_in:
                if line.startswith("ATOM") or line.startswith("TER"):
                    f_out.write(line)
            f_out.write("END\n")

    def prep_autodock(self):
        # Using quotes for space-safe paths
        cmd = f'prepare_receptor4.py -r "{self.clean_pdb}" -o "{self.pdbqt}" -A hydrogens -U nwaters'
        subprocess.run(cmd, shell=True, capture_output=True)

    def prep_ledock(self):
        bin_path = os.path.abspath("bin/lepro")
        subprocess.run(f'"{bin_path}" "{os.path.basename(self.clean_pdb)}"', 
                       shell=True, cwd=self.output_path, capture_output=True)

    def validate_outputs(self):
        """Strict check to ensure all required files exist and are not empty."""
        requirements = {
            "Clean PDB": self.clean_pdb,
            "AutoDock PDBQT": self.pdbqt,
            "LeDock pro.pdb": self.pro_pdb
        }
        success = True
        print("\n--- Phase 1: Preprocessing Validation ---")
        for name, path in requirements.items():
            if os.path.exists(path) and os.path.getsize(path) > 0:
                print(f"PASS: {name} found ({os.path.getsize(path) // 1024} KB)")
            else:
                print(f"FAIL: {name} is missing or empty.")
                success = False
        return success

    def run_all(self):
        self.fetch_pdb()
        self.clean_structure()
        self.prep_autodock()
        self.prep_ledock()
        return self.validate_outputs()
