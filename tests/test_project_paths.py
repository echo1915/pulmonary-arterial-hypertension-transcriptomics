from pathlib import Path
import sys


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from project_paths import DATA_ROOT, OUTPUT_ROOT, PROJECT_DATA_ROOT, PROJECT_SLUG


def main():
    assert PROJECT_SLUG == "project-001-pulmonary-arterial-hypertension-transcriptomics"
    assert DATA_ROOT == PROJECT_DATA_ROOT / "raw" / "current-data-snapshot"
    assert OUTPUT_ROOT == PROJECT_DATA_ROOT / "outputs" / "current-results"
    assert DATA_ROOT.exists(), f"Missing migrated data root: {DATA_ROOT}"
    print(f"project_data_root={PROJECT_DATA_ROOT}")
    print(f"data_root={DATA_ROOT}")
    print(f"output_root={OUTPUT_ROOT}")


if __name__ == "__main__":
    main()

