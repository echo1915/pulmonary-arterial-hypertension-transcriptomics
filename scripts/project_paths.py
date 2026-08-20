"""Centralized project paths for reproducible local runs."""

from __future__ import annotations

import os
from pathlib import Path


PROJECT_SLUG = "project-001-pulmonary-arterial-hypertension-transcriptomics"
PROJECT_ROOT = Path(__file__).resolve().parents[1]
WORKDATA_ROOT = Path(os.environ.get("WORKDATA_ROOT", r"G:\workdata"))
PROJECT_DATA_ROOT = Path(
    os.environ.get("PAH_DATA_ROOT", WORKDATA_ROOT / "projects" / PROJECT_SLUG)
)
DATA_ROOT = PROJECT_DATA_ROOT / "raw" / "current-data-snapshot"
OUTPUT_ROOT = PROJECT_DATA_ROOT / "outputs" / "current-results"
INTERIM_ROOT = PROJECT_DATA_ROOT / "interim"


def ensure_output_dirs() -> None:
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    INTERIM_ROOT.mkdir(parents=True, exist_ok=True)



