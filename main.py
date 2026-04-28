import os
from src.preprocessor import PoseAIPreprocessor
from src.detector import PocketDetector
from src.docking import DockingManager

def main():
    target = "1iep"
    print(f"--- PoseAI Pipeline: {target.upper()} ---")
    
    # 1. Preprocessing
    prep = PoseAIPreprocessor(target)
    if not prep.run_all():
        return

    # 2. Detection
    detector = PocketDetector(target)
    center = detector.detect_from_ligand()
    
    if center is not None:
        # For testing, we use the raw PDB as the ligand source
        # Smina is smart enough to extract HETATMs if pointed at the raw PDB
        ligand_ref = os.path.abspath(os.path.join("data", target, f"{target}.pdb"))
        
        # 3. Docking
        docker = DockingManager(target, center)
        if docker.run_all(ligand_ref):
            print("\n--- Pipeline Complete: Results ready for Consensus Scoring ---")
        else:
            print("\n--- Pipeline Error: One or more engines failed ---")
    else:
        print("Error: No binding site detected.")

if __name__ == "__main__":
    main()
