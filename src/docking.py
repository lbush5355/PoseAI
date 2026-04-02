import os
import subprocess

class DockingManager:
    def __init__(self, pdb_id, center, box_size=20.0, work_dir="data"):
        self.pdb_id = pdb_id.lower()
        self.center = center
        self.box_size = box_size
        self.work_dir = os.path.abspath(os.path.join(work_dir, self.pdb_id))
        self.bin_dir = os.path.abspath("bin")

    def _write_ledock_config(self):
        """Updates the dock.in file with Cartesian boundaries."""
        conf_path = os.path.join(self.work_dir, "dock.in")
        half = self.box_size / 2
        
        xmin, xmax = self.center[0] - half, self.center[0] + half
        ymin, ymax = self.center[1] - half, self.center[1] + half
        zmin, zmax = self.center[2] - half, self.center[2] + half

        if not os.path.exists(conf_path):
            print("Error: dock.in template missing.")
            return

        with open(conf_path, 'r') as f:
            lines = f.readlines()

        new_lines = []
        for line in lines:
            if line.startswith("-x"):
                new_lines.append(f"-x {xmin:.3f} {xmax:.3f}\n")
            elif line.startswith("-y"):
                new_lines.append(f"-y {ymin:.3f} {ymax:.3f}\n")
            elif line.startswith("-z"):
                new_lines.append(f"-z {zmin:.3f} {zmax:.3f}\n")
            else:
                new_lines.append(line)

        with open(conf_path, 'w') as f:
            f.writelines(new_lines)

    def run_smina(self, ligand_file):
        print("Action: Running Smina simulation...")
        receptor = os.path.join(self.work_dir, f"{self.pdb_id}_clean.pdbqt")
        output = os.path.join(self.work_dir, "out_smina.pdbqt")
        
        cmd = (
            f'smina -r "{receptor}" -l "{ligand_file}" '
            f'--center_x {self.center[0]} --center_y {self.center[1]} --center_z {self.center[2]} '
            f'--size_x {self.box_size} --size_y {self.box_size} --size_z {self.box_size} '
            f'--out "{output}" --exhaustiveness 8 --cpu 4'
        )
        
        # REMOVED capture_output=True to allow direct terminal feedback
        # This prevents the pipe buffer from hanging on M2
        subprocess.run(cmd, shell=True, check=True) 
        return os.path.exists(output)

    def run_ledock(self):
        print("Action: Running LeDock simulation...")
        self._write_ledock_config()
        ledock_bin = os.path.join(self.bin_dir, "ledock")
        
        try:
            subprocess.run(f'"{ledock_bin}" dock.in', shell=True, check=True, cwd=self.work_dir, capture_output=True)
            # LeDock generates output as {protein_name}.dok. For us, usually pro.pdb.dok
            return os.path.exists(os.path.join(self.work_dir, "pro.pdb.dok"))
        except Exception as e:
            print(f"LeDock Runtime Error: {e}")
            return False

    def run_all(self, ligand_file):
        s_ok = self.run_smina(ligand_file)
        l_ok = self.run_ledock()
        
        print("\n--- Phase 2: Docking Execution Summary ---")
        print(f"Smina Success: {s_ok}")
        print(f"LeDock Success: {l_ok}")
        return s_ok and l_ok
