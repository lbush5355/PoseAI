import os
import subprocess
import urllib.request

class PoseAIPreprocessor:
    def __init__(self, pdb_id, work_dir="data"):
        self.pdb_id = pdb_id.lower()
        self.output_path = os.path.join(work_dir, self.pdb_id)
        self.raw_pdb = os.path.join(self.output_path, f"{self.pdb_id}.pdb")
        self.clean_pdb = os.path.join(self.output_path, f"{self.pdb_id}_clean.pdb")
        
        if not os.path.exists(self.output_path):
            os.makedirs(self.output_path)

    def fetch_pdb(self):
        if os.path.exists(self.raw_pdb):
            print(f"Status: {self.pdb_id.upper()}.pdb already exists.")
        else:
            url = f"https://files.rcsb.org/download/{self.pdb_id.upper()}.pdb"
            print(f"Action: Downloading {self.pdb_id.upper()}...")
            urllib.request.urlretrieve(url, self.raw_pdb)

    def clean_structure(self):
        print("Action: Cleaning PDB structure...")
        with open(self.raw_pdb, "r") as f_in, open(self.clean_pdb, "w") as f_out:
            for line in f_in:
                if line.startswith("ATOM") or line.startswith("TER"):
                    f_out.write(line)
            f_out.write("END\n")

    def prep_autodock(self):
        print("Action: Generating PDBQT for AutoDock 4...")
        output_pdbqt = f"{self.clean_pdb}qt"
        # Use quotes around paths to handle spaces in folder names
        cmd = f'prepare_receptor4.py -r "{self.clean_pdb}" -o "{output_pdbqt}" -A hydrogens -U nwaters'
        try:
            subprocess.run(cmd, shell=True, check=True, capture_output=True)
        except Exception as e:
            print(f"AutoDock Prep Warning: {e}")

    def prep_ledock(self):
        print("Action: Processing for LeDock...")
        bin_path = os.path.abspath("bin/lepro")
        # Quote the binary path and the input file path
        cmd = f'"{bin_path}" "{os.path.basename(self.clean_pdb)}"'
        try:
            subprocess.run(cmd, shell=True, check=True, cwd=self.output_path)
            print(f"Status: lepro generated pro.pdb and dock.in")
        except Exception as e:
            print(f"LeDock Prep Error: {e}")

    def run_all(self):
        self.fetch_pdb()
        self.clean_structure()
        self.prep_autodock()
        self.prep_ledock()
        print(f"Result: Preprocessing for {self.pdb_id} complete.")
