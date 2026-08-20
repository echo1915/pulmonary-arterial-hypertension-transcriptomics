"""Build the separately organized Figure 5 supplementary package (review-v2)."""

from __future__ import annotations

import shutil
import sys
from pathlib import Path

_REPO = Path(__file__).resolve().parents[1]
_LOCAL_DEPS = _REPO / ".tmp" / "figure5-python-deps"
if _LOCAL_DEPS.exists():
    sys.path.insert(0, str(_LOCAL_DEPS))

import matplotlib as mpl

sys.path.insert(0, str(Path(__file__).resolve().parent))
import plot_figure5_program_review as v1  # noqa: E402
from project_paths import PROJECT_DATA_ROOT  # noqa: E402


OUT = PROJECT_DATA_ROOT / "outputs" / "figure5-program-review-v2" / "supplementary"

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "font.size": 7,
})


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    v1.OUT = OUT
    frames = v1.load_inputs()
    v1.build_supplement(frames)
    v1.build_robustness_heatmaps()

    # Rename copies into an ordered manuscript-facing supplement sequence while
    # retaining the exact reproducible v1-derived render in the same directory.
    mapping = {
        "Supplementary_Figure_PDE8B_program_support_review_v1":
            "Supplementary_Figure_5_1_PDE8B_bridge_and_annotation_QC_review_v2",
        "Supplementary_Figure_PDE8B_program_robustness_review_v1":
            "Supplementary_Figure_5_2_PDE8B_gene_level_consistency_review_v2",
    }
    for source_stem, target_stem in mapping.items():
        for suffix in (".svg", ".pdf", ".png", ".tiff"):
            source = OUT / f"{source_stem}{suffix}"
            target = OUT / f"{target_stem}{suffix}"
            shutil.copy2(source, target)
            source.unlink()

    manifest = """# Figure 5 supplementary package — review-v2

## Supplementary Figure 5.1 — PDE8B bridge and annotation QC

- a, Directional relationship between PAH differential effects and PDE8B meta-correlations among the mapped robust-signature genes.
- b, Explicit overlap tests showing no gene-set membership enrichment; the stable mechanism module overlaps by only six genes.
- c, Descriptive biological-theme convergence, not evidence that the PDE8B program is inherited from the 188-gene signature.
- d–e, GSE293580 whole-atlas and mural-cell integration/annotation QC. These are cell-level QC views and are not disease-effect tests.

## Supplementary Figure 5.2 — PDE8B gene-level consistency

- a, Patient-level SMC pseudobulk heatmap for selected representative PDE8B-associated genes.
- b, The same selected genes' correlations with PDE8B across the four primary non-overlapping cohorts.
- This figure is gene-specific supporting evidence. Main Figure 5e instead summarizes all 201 core genes and their leave-one-cohort-out sign stability.

## Excluded from the package

- The decorative 12-gene network.
- Cell-state trends already assigned to tentative Figure 4.
- The SMC–pericyte ligand–receptor proxy that did not pass FDR and cannot support a communication claim.
- Historical files and source data remain preserved elsewhere; exclusion is compositional, not deletion of project evidence.

## Claim boundaries

- PDE8B remains an associated mechanism candidate, not a causal regulator or validated therapeutic target.
- Bridge results support directional concordance only; they do not establish gene-set inheritance.
- UMAP panels are annotation/integration QC only.
"""
    (OUT / "Figure_5_supplementary_review_v2_manifest.md").write_text(manifest, encoding="utf-8")


if __name__ == "__main__":
    main()
    for path in sorted(OUT.iterdir()):
        print(path)
