from src.preprocessor import PoseAIPreprocessor
from src.detector import PocketDetector

def test_pipeline():
    target = "1iep"
    print(f"--- PoseAI Pipeline: {target.upper()} ---")
    
    # 1. Preprocess
    prep = PoseAIPreprocessor(target)
    prep.run_all()
    
    # 2. Detect Coordinates
    detector = PocketDetector(target)
    center = detector.detect_from_ligand()
    
    if center is not None:
        print(f"Result: Binding Site Centroid at X:{center[0]} Y:{center[1]} Z:{center[2]}")
    else:
        print("Result: Automated detection failed. Fallback required.")

if __name__ == "__main__":
    test_pipeline()