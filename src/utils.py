import os
import subprocess

def convert_ligand(input_file, output_format="pdbqt"):
    """
    Uses Open Babel to efficiently convert any ligand format (sdf, mol2) 
    into a docking-ready format.
    """
    base = os.path.splitext(input_file)[0]
    output_file = f"{base}.{output_format}"
    
    # Efficient shell command for Open Babel
    # -p 7.4 adds hydrogens at physiological pH
    cmd = f'obabel "{input_file}" -O "{output_file}" -p 7.4 --gen3d'
    
    try:
        subprocess.run(cmd, shell=True, check=True, capture_output=True)
        return output_file
    except Exception as e:
        print(f"Conversion Error: {e}")
        return None