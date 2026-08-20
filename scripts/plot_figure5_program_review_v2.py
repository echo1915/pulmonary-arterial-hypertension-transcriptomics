"""Build the six-panel Figure 5 review-v2 candidate.

This revision preserves review-v1 panels a-d and adds two non-causal support
panels: cross-cohort/core-module stability (genes are the descriptive unit) and
patient-level PDE8B localization in GSE210248 (patients are the display unit).
"""

from __future__ import annotations

import sys
from pathlib import Path

# The project uses an embeddable Python runtime that ignores PYTHONPATH.
# Permit a repository-local, untracked dependency bundle without modifying it.
_LOCAL_DEPS = Path(__file__).resolve().parents[1] / ".tmp" / "figure5-python-deps"
if _LOCAL_DEPS.exists():
    sys.path.insert(0, str(_LOCAL_DEPS))

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import plot_figure5_program_review as v1  # noqa: E402
from project_paths import PROJECT_DATA_ROOT  # noqa: E402


ROOT = PROJECT_DATA_ROOT
O = ROOT / "outputs"
OUT = O / "figure5-program-review-v2"
COEX = O / "pde8b-coexpression"

DONOR = "#0072B2"
PAH = "#D55E00"
INK = "#263746"
MUTED = "#68757F"
NEUTRAL = "#D9DDE0"
GRID = "#E2E6E9"
POSITIVE = "#D58A42"
NEGATIVE = "#5B7FA3"

mpl.rcParams.update(v1.mpl.rcParams)
mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "font.size": 7,
})


def save_bundle(fig: plt.Figure, stem: Path) -> None:
    """Explicit publication exports (kept local so static QA can audit them)."""
    stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".tiff"), dpi=600, bbox_inches="tight", facecolor="white",
                pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)


def plot_core_stability(ax: plt.Axes, co: pd.DataFrame, loco: pd.DataFrame) -> pd.DataFrame:
    """Show cohort distributions and LOCO ranges without treating genes as subjects."""
    core = co.loc[co["core_module"].astype(bool)].copy()
    cohorts = ["GSE15197", "GSE113439", "GSE254617", "GSE208592"]
    cols = [f"r_{x}" for x in cohorts]
    long = core[["symbol", "adjusted_meta_r", *cols]].melt(
        id_vars=["symbol", "adjusted_meta_r"], value_vars=cols,
        var_name="cohort", value_name="r"
    )
    long["cohort"] = long["cohort"].str.replace("r_", "", regex=False)
    long["sign"] = np.where(long["adjusted_meta_r"].ge(0), "Positive core", "Negative core")

    positions = np.arange(len(cohorts))
    for i, cohort in enumerate(cohorts):
        q = long[long.cohort.eq(cohort)]
        for sign, color, dx in [("Negative core", NEGATIVE, -0.08), ("Positive core", POSITIVE, 0.08)]:
            vals_all = q.loc[q.sign.eq(sign), "r"]
            missing_n = int(vals_all.isna().sum())
            if missing_n:
                raise ValueError(f"Unexpected missing cohort correlations: {cohort}, {sign}, n={missing_n}")
            vals = vals_all.to_numpy()
            bp = ax.boxplot(vals, positions=[i + dx], widths=.13, patch_artist=True,
                            showfliers=False, whis=(10, 90), manage_ticks=False)
            bp["boxes"][0].set(facecolor=color, edgecolor=color, alpha=.55, linewidth=.65)
            for item in bp["whiskers"] + bp["caps"]:
                item.set(color=color, linewidth=.65)
            bp["medians"][0].set(color=INK, linewidth=.9)
    ax.axhline(0, color=MUTED, lw=.65)
    ax.set_xticks(positions, cohorts, rotation=27, ha="right")
    ax.set_ylabel("Gene-wise correlation with PDE8B")
    ax.set_title("Core-program stability across cohorts", loc="left")
    ax.grid(axis="y", color=GRID, lw=.45)

    loco_core = loco[loco.symbol.isin(core.symbol)].copy()
    stable = loco_core.groupby("symbol", as_index=False).agg(
        n_omissions=("omitted_cohort", "nunique"),
        min_loco_r=("meta_r", "min"), max_loco_r=("meta_r", "max")
    )
    stable = stable.merge(core[["symbol", "adjusted_meta_r"]], on="symbol", how="left")
    if len(stable) != len(core) or not stable.n_omissions.eq(4).all():
        raise ValueError("LOCO stability rows do not cover every core gene across all four omissions")
    return long


def plot_celltype_localization(ax: plt.Axes, pseudo: pd.DataFrame) -> pd.DataFrame:
    z = pseudo[(pseudo.gene.eq("PDE8B")) & (pseudo.n_cells.ge(5))].copy()
    summary = z.groupby("cell_type").agg(
        mean_expr=("log_cp10k", "mean"), eligible_patients=("sample", "nunique")
    ).sort_values("mean_expr", ascending=True)
    order = summary.index.tolist()
    ymap = {cell: i for i, cell in enumerate(order)}
    for disease, color, marker in [("Donor", DONOR, "o"), ("PAH", PAH, "o")]:
        q = z[z.disease.eq(disease)].copy()
        q["y"] = q.cell_type.map(ymap)
        ax.scatter(q.log_cp10k, q.y, s=21, color=color, marker=marker,
                   edgecolor="white", linewidth=.5, alpha=.9, label=disease, zorder=3)
    for i, cell in enumerate(order):
        q = z[z.cell_type.eq(cell)]
        ax.hlines(i, q.log_cp10k.min(), q.log_cp10k.max(), color=NEUTRAL, lw=.7, zorder=1)
    ax.set_yticks(range(len(order)), order)
    ax.set_xlabel("PDE8B expression (log CP10K)")
    ax.set_title("Patient-level cell-type localization", loc="left")
    ax.grid(axis="x", color=GRID, lw=.45)
    ax.legend(loc="lower right", fontsize=5.4, ncol=2, handletextpad=.25, columnspacing=.7)
    return z


def build_main() -> None:
    frames = v1.load_inputs()
    co = frames["coexpression"]
    pathways = v1.prepare_selected_pathways(frames["pathways"])
    classic = frames["classic"]
    values = frames["gse_values"]
    tests = frames["gse_tests"]
    loco = pd.read_csv(COEX / "pde8b_module_leave_one_cohort_out.csv")
    pseudo = pd.read_csv(O / "current-results" / "pah_audit" / "scrna_gse210248_candidate_pseudobulk.csv")

    fig = plt.figure(figsize=(7.35, 8.35))
    gs = fig.add_gridspec(3, 2, left=.085, right=.98, top=.94, bottom=.08,
                          width_ratios=[1.28, 1.0], height_ratios=[1.05, .98, 1.05],
                          wspace=.40, hspace=.58)
    axes = [fig.add_subplot(gs[i, j]) for i in range(3) for j in range(2)]
    a, b, c, d, e, f = axes

    # Reproduce the author's reviewed unified-v3 panel placement in physical
    # points.  The review PDF uses a wider right column and progressively
    # smaller vertical offsets across the three rows.
    panel_offsets = {
        a: (1.0, -24.0),
        b: (46.0, -24.0),
        c: (16.0, -13.0),
        d: (13.0, -13.0),
        e: (16.0, -8.0),
        f: (13.0, -8.0),
    }
    for ax, (dx_points, dy_points) in panel_offsets.items():
        pos = ax.get_position()
        ax.set_position([
            pos.x0 + dx_points / (72 * fig.get_figwidth()),
            pos.y0 + dy_points / (72 * fig.get_figheight()),
            pos.width,
            pos.height,
        ])
    v1.plot_coexpression(a, co); v1.panel_label(a, "a")
    v1.plot_pathways(b, pathways); v1.panel_label(b, "b", x_fig=.559)
    v1.plot_classic_connections(c, classic); v1.panel_label(c, "c")
    v1.plot_gse293580(d, values, tests); v1.panel_label(d, "d", x_fig=.559)
    stability = plot_core_stability(e, co, loco); v1.panel_label(e, "e")
    localization = plot_celltype_localization(f, pseudo); v1.panel_label(f, "f", x_fig=.559)

    stem = OUT / "Figure_5_PDE8B_program_main_unified_v3"
    save_bundle(fig, stem)

    OUT.mkdir(parents=True, exist_ok=True)
    stability.to_csv(OUT / "Figure_5_panel_e_core_stability_source_data.csv", index=False)
    localization.to_csv(OUT / "Figure_5_panel_f_celltype_localization_source_data.csv", index=False)
    pathways.to_csv(OUT / "Figure_5_panel_b_selected_pathways_source_data.csv", index=False)
    classic.to_csv(OUT / "Figure_5_panel_c_PAH_gene_connections_source_data.csv", index=False)
    values[(values.compartment == "strict_SMC") & values.n_cells.ge(5)].to_csv(
        OUT / "Figure_5_panel_d_GSE293580_patient_values_source_data.csv", index=False)

    notes = """# Figure 5 review-v2 QA notes

- Core conclusion remains associative: PDE8B covaries with a reproducible vascular/cyclic-nucleotide program and receives descriptive patient-level SMC support.
- Panels a-d are inherited from review-v1 without changing their statistical meaning.
- Panel e summarizes the distribution of core-gene correlations in each cohort and leave-one-cohort-out sign stability. Genes are explicitly descriptive units, not biological replicates; no inferential P value is computed from genes.
- Panel f uses GSE210248 patient-by-cell-type pseudobulk values. Every displayed point is a patient; patient-cell-type combinations with fewer than 5 cells are excluded. The panel is localization evidence, not a disease-effect test.
- Figure 3 vocabulary is retained: Donor/SMC2 blue #0072B2, PAH/SMC1 orange-red #D55E00; circles denote Donor and IPAH/PAH, triangles denote SSc-PAH where present.
- GSE293580 strict-SMC remains descriptive (Donor n=2, IPAH n=3, SSc-PAH n=4); exact P/q values and individual patient points remain visible.
- No Figure 3 or Figure 4 files were modified. Historical review-v1 outputs were preserved.
"""
    (OUT / "Figure_5_review_v2_QA_notes.md").write_text(notes, encoding="utf-8")


if __name__ == "__main__":
    build_main()
    for path in sorted(OUT.iterdir()):
        print(path)
