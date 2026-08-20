"""Build the revised Figure 5 main and supplementary review figures.

Figure contract
---------------
Core conclusion:
    PDE8B reproducibly covaries with a vascular/cyclic-nucleotide program across
    four lung cohorts and receives descriptive patient-level support in the small
    independent GSE293580 SMC dataset.
Archetype:
    Asymmetric quantitative grid with a coexpression hero panel.
Statistical unit:
    Subjects within cohorts for coexpression; cohorts for random-effects synthesis;
    patients for GSE293580 replication; genes/pathway sets for enrichment/bridge panels.
Claim boundary:
    Association and directional replication only; no causal regulation or cell-cell
    communication inference.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
import sys

_LOCAL_DEPS = Path(__file__).resolve().parents[1] / ".tmp" / "figure5-python-deps"
if _LOCAL_DEPS.exists():
    sys.path.insert(0, str(_LOCAL_DEPS))

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import Normalize

sys.path.insert(0, str(Path(__file__).resolve().parent))
from project_paths import PROJECT_DATA_ROOT  # noqa: E402


ROOT = PROJECT_DATA_ROOT
O = ROOT / "outputs"
OUT = O / "figure5-program-review-v1"
COEX = O / "pde8b-coexpression"

# Inherit the confirmed Figure 3 visual vocabulary.
DONOR = "#0072B2"
PAH = "#D55E00"
SMC1 = "#D55E00"
SMC2 = "#0072B2"
INK = "#263746"
MUTED = "#68757F"
NEUTRAL = "#D9DDE0"
GRID = "#E2E6E9"
NEGATIVE = "#5B7FA3"
POSITIVE = "#D58A42"

mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "font.size": 7,
        "axes.titlesize": 8,
        "axes.titleweight": "semibold",
        "axes.labelsize": 7,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.75,
        "legend.frameon": False,
        "xtick.labelsize": 6.5,
        "ytick.labelsize": 6.5,
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
    }
)


def panel_label(ax, label: str, x: float = -0.10, y: float = 1.10,
                x_fig: float | None = None) -> None:
    fig = ax.figure
    y_fig = fig.transFigure.inverted().transform(ax.transAxes.transform((0, y)))[1]
    if x_fig is None:
        x_fig = .015 if ax.get_position().x0 < .50 else .55
    fig.text(x_fig, y_fig, label.upper(), fontsize=9, fontweight=700,
             va="bottom", ha="left", color=INK)


def shift_axes(ax, dx_points: float) -> None:
    """Reproduce the author's horizontal panel adjustment in physical points."""
    pos = ax.get_position()
    dx = dx_points / (72 * ax.figure.get_figwidth())
    ax.set_position([pos.x0 + dx, pos.y0, pos.width, pos.height])


def fixed_jitter(n: int) -> np.ndarray:
    return np.linspace(-0.055, 0.055, n) if n > 1 else np.zeros(1)


def save_bundle(fig, stem: Path) -> None:
    stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".tiff"), dpi=600, bbox_inches="tight", facecolor="white",
                pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)


def load_inputs() -> dict[str, pd.DataFrame | dict]:
    frames: dict[str, pd.DataFrame | dict] = {
        "coexpression": pd.read_csv(COEX / "pde8b_cross_cohort_coexpression_meta.csv"),
        "pathways": pd.read_csv(COEX / "enrichment" / "pde8b_core_pathway_enrichment.csv"),
        "classic": pd.read_csv(COEX / "state-deconvolution" / "pde8b_classic_pah_gene_connections.csv"),
        "gse_values": pd.read_csv(ROOT / "processed" / "gse293580_scanpy" / "gse293580_pde8b_patient_pseudobulk.csv"),
        "gse_tests": pd.read_csv(O / "gse293580-reanalysis" / "gse293580_pde8b_patient_level_tests.csv"),
        "bridge": pd.read_csv(O / "signature-pde8b-bridge" / "robust188_pde8b_bridge_gene_table.csv"),
        "overlap": pd.read_csv(O / "signature-pde8b-bridge" / "robust188_pde8b_overlap_tests.csv"),
        "theme": pd.read_csv(O / "pathway-bridge" / "pathway_theme_bridge_matrix.csv"),
        "atlas_umap": pd.read_csv(ROOT / "processed" / "gse293580_scanpy" / "gse293580_harmony_umap_clusters.csv.gz"),
        "mural_umap": pd.read_csv(ROOT / "processed" / "gse293580_scanpy" / "gse293580_mural_harmony_umap.csv.gz"),
    }
    with (O / "signature-pde8b-bridge" / "robust188_pde8b_bridge_summary.json").open(encoding="utf-8") as handle:
        frames["bridge_summary"] = json.load(handle)
    return frames


def prepare_selected_pathways(pathways: pd.DataFrame) -> pd.DataFrame:
    selected = [
        ("core_positive", "cilium organization", "Cilium organization"),
        ("core_positive", "muscle system process", "Muscle system process"),
        ("core_positive", "actin filament-based process", "Actin filament process"),
        ("core_positive", "cGMP-PKG signaling pathway", "cGMP–PKG signaling"),
        ("core_positive", "smooth muscle contraction", "Smooth-muscle contraction"),
        ("core_positive", "cellular response to calcium ion", "Calcium-ion response"),
        ("core_negative", "Translation", "Translation"),
        ("core_negative", "Cellular responses to stress", "Stress response"),
    ]
    rows = []
    for query, term, label in selected:
        hit = pathways[(pathways.query_set == query) & (pathways.term_name == term)]
        if len(hit) != 1:
            raise ValueError(f"Expected one pathway row for {query}: {term}; found {len(hit)}")
        row = hit.iloc[0].copy()
        row["display"] = label
        rows.append(row)
    return pd.DataFrame(rows)


def plot_coexpression(ax, co: pd.DataFrame) -> None:
    co = co.copy()
    co["mlog"] = -np.log10(co.adjusted_fdr.clip(lower=1e-300))
    core = co.core_module.astype(bool)
    ax.scatter(co.loc[~core, "adjusted_meta_r"], co.loc[~core, "mlog"], s=3.3,
               color=NEUTRAL, alpha=.45, linewidths=0, rasterized=True)
    neg = core & co.adjusted_meta_r.lt(0)
    pos = core & co.adjusted_meta_r.gt(0)
    ax.scatter(co.loc[neg, "adjusted_meta_r"], co.loc[neg, "mlog"], s=7,
               color=NEGATIVE, alpha=.78, linewidths=0, rasterized=True,
               label=f"Negative core (n={neg.sum()})")
    ax.scatter(co.loc[pos, "adjusted_meta_r"], co.loc[pos, "mlog"], s=7,
               color=POSITIVE, alpha=.78, linewidths=0, rasterized=True,
               label=f"Positive core (n={pos.sum()})")
    ax.axvline(0, color=MUTED, lw=.7)
    ax.set_xlabel("Meta-correlation with PDE8B")
    ax.set_ylabel("−log$_{10}$(FDR)")
    ax.set_title("Stable cross-cohort coexpression", loc="left")
    ax.legend(loc="upper left", fontsize=5.8, markerscale=1.5, handletextpad=.25)
    labels = ["PTPRD", "TRPC3", "PRKG1", "PDE5A", "IRAK1"]
    offsets = {"PTPRD": (4, 3), "TRPC3": (4, -7), "PRKG1": (4, 3), "PDE5A": (4, -7), "IRAK1": (-29, 3)}
    for gene in labels:
        q = co.loc[co.symbol.eq(gene)]
        if q.empty:
            continue
        r = q.iloc[0]
        ax.annotate(gene, (r.adjusted_meta_r, r.mlog), xytext=offsets[gene],
                    textcoords="offset points", fontsize=5.6, color=INK)


def plot_pathways(ax, selected: pd.DataFrame) -> None:
    z = selected.copy()
    z["score"] = -np.log10(z.fdr)
    z = z.sort_values(["query_set", "score"], ascending=[True, True])
    y = np.arange(len(z))
    colors = np.where(z.query_set.eq("core_positive"), POSITIVE, NEGATIVE)
    sizes = 24 + 120 * z.precision.astype(float)
    ax.hlines(y, 0, z.score, color=colors, lw=1.4, alpha=.70)
    ax.scatter(z.score, y, s=sizes, c=colors, edgecolor="white", linewidth=.55, zorder=3)
    ax.axvline(-np.log10(.05), color=MUTED, lw=.7, ls="--")
    ax.set_yticks(y, z.display)
    ax.set_xlabel("−log$_{10}$(FDR)")
    ax.set_title("Selected core-module pathways", loc="left")
    handles = [
        ax.scatter([], [], s=36, color=POSITIVE, edgecolor="white", label="Positive core"),
        ax.scatter([], [], s=36, color=NEGATIVE, edgecolor="white", label="Negative core"),
    ]
    ax.legend(handles=handles, loc="lower right", fontsize=5.7, handletextpad=.25)


def plot_classic_connections(ax, classic: pd.DataFrame) -> None:
    genes = ["PRKG1", "PDE5A", "TRPC3", "PLCB4", "ROCK2", "EDNRA", "BMPR2", "SMAD9", "CAV1", "KDR", "SOX17", "ACVRL1"]
    z = classic.set_index("symbol").reindex(genes)
    missing = z.adjusted_meta_r.isna()
    if missing.any():
        raise ValueError(f"Missing predefined PAH-gene estimates: {z.index[missing].tolist()}")
    z = z.sort_values("adjusted_meta_r")
    y = np.arange(len(z))
    x = z.adjusted_meta_r.to_numpy()
    lo = z.adjusted_ci_low_r.to_numpy()
    hi = z.adjusted_ci_high_r.to_numpy()
    colors = np.where(z.adjusted_fdr.lt(.05), PAH, MUTED)
    ax.axvline(0, color=MUTED, lw=.7)
    for xx, yy, low, high, color in zip(x, y, lo, hi, colors):
        ax.errorbar(xx, yy, xerr=[[xx-low], [high-xx]], fmt="none", ecolor=color,
                    elinewidth=1.0, capsize=1.8, alpha=.85)
    ax.scatter(x, y, s=27, c=colors, marker="o", edgecolor="white", linewidth=.55, zorder=3)
    ax.set_yticks(y, z.index)
    ax.set_xlabel("Meta-correlation with PDE8B (95% CI)")
    ax.set_title("Predefined PAH and cyclic-nucleotide genes", loc="left")


def plot_gse293580(ax, values: pd.DataFrame, tests: pd.DataFrame) -> None:
    # Include every annotated donor in this descriptive strict-SMC display.
    # A donor with no retained strict-SMC cells is shown at zero detected PDE8B.
    z = values[values.compartment == "strict_SMC"].copy()
    z["PDE8B_log1p_CPM"] = z["PDE8B_log1p_CPM"].fillna(0.0)
    groups = ["Donor", "IPAH", "SSc-PAH"]
    colors = [DONOR, PAH, PAH]
    markers = ["o", "o", "^"]
    for xpos, (group, color, marker) in enumerate(zip(groups, colors, markers)):
        q = z[z.condition.eq(group)].sort_values("sample")
        ax.scatter(xpos + fixed_jitter(len(q)), q.PDE8B_log1p_CPM, s=39,
                   color=color, marker=marker, edgecolor="white", linewidth=.65,
                   alpha=.92, zorder=3)
        ax.hlines(q.PDE8B_log1p_CPM.mean(), xpos-.18, xpos+.18, color=INK, lw=1.35)
    group_counts = z.groupby("condition")["sample"].nunique()
    ax.set_xticks(
        range(3),
        [f"{group}\n(n={int(group_counts.get(group, 0))})" for group in groups],
    )
    ax.set_ylabel("PDE8B expression\n(log$_e$[CPM+1])")
    ax.set_title("Independent GSE293580 strict-SMC replication", loc="left")
    ax.grid(axis="y", color=GRID, lw=.55)



def build_main(frames: dict[str, pd.DataFrame | dict]) -> None:
    co = frames["coexpression"]
    selected = prepare_selected_pathways(frames["pathways"])
    classic = frames["classic"]
    values = frames["gse_values"]
    tests = frames["gse_tests"]
    assert isinstance(co, pd.DataFrame) and isinstance(selected, pd.DataFrame)
    assert isinstance(classic, pd.DataFrame) and isinstance(values, pd.DataFrame) and isinstance(tests, pd.DataFrame)

    fig = plt.figure(figsize=(7.35, 6.55))
    gs = fig.add_gridspec(2, 2, left=.085, right=.98, top=.93, bottom=.095,
                          width_ratios=[1.30, 1.0], height_ratios=[1.08, 1.0],
                          wspace=.39, hspace=.46)
    a = fig.add_subplot(gs[0, 0])
    b = fig.add_subplot(gs[0, 1]); shift_axes(b, 26)
    c = fig.add_subplot(gs[1, 0]); shift_axes(c, 9)
    d = fig.add_subplot(gs[1, 1]); shift_axes(d, -6)
    plot_coexpression(a, co); panel_label(a, "a")
    plot_pathways(b, selected); panel_label(b, "b", x_fig=.535)
    plot_classic_connections(c, classic); panel_label(c, "c")
    plot_gse293580(d, values, tests); panel_label(d, "d", x_fig=.535)
    # Preserve the vertical trim of the author's edited PDF after removing the
    # visible page title and footer sentence.
    fig.add_artist(mpl.lines.Line2D([.5, .5], [.025, .983], transform=fig.transFigure,
                                   color="none", linewidth=.01))
    save_bundle(fig, OUT / "Figure_5_PDE8B_program_main_review_v1")

    selected.to_csv(OUT / "Figure_5_panel_b_selected_pathways_source_data.csv", index=False)
    classic.to_csv(OUT / "Figure_5_panel_c_PAH_gene_connections_source_data.csv", index=False)
    panel_d_values = values[values.compartment == "strict_SMC"].copy()
    panel_d_values["PDE8B_log1p_CPM"] = panel_d_values["PDE8B_log1p_CPM"].fillna(0.0)
    panel_d_values.to_csv(
        OUT / "Figure_5_panel_d_GSE293580_patient_values_source_data.csv", index=False
    )


def plot_bridge(ax, bridge: pd.DataFrame, summary: dict) -> None:
    eligible = bridge[["meta_g", "adjusted_meta_r"]].notna().all(axis=1)
    excluded_n = int((~eligible).sum())
    expected_rows = int(summary["robust_in_coexpression_universe"])
    # The exported table is already restricted to mapped genes. PDE8B self-correlation
    # and DMRTC2 are documented only in the summary, so enforce the mapped row count.
    if excluded_n != 0 or len(bridge) != expected_rows:
        raise ValueError(f"Expected {expected_rows} complete bridge rows; found {len(bridge)} with {excluded_n} incomplete")
    z = bridge.loc[eligible].copy()
    stable = z.in_stable_mechanism.astype(bool)
    ax.scatter(z.loc[~stable, "meta_g"], z.loc[~stable, "adjusted_meta_r"], s=10,
               color=NEUTRAL, alpha=.70, linewidths=0)
    ax.scatter(z.loc[stable, "meta_g"], z.loc[stable, "adjusted_meta_r"], s=22,
               color=POSITIVE, edgecolor="white", linewidth=.45)
    ax.axhline(0, color=MUTED, lw=.65); ax.axvline(0, color=MUTED, lw=.65)
    ax.set_xlabel("PAH differential effect (Hedges g)")
    ax.set_ylabel("Meta-correlation with PDE8B")
    ax.set_title("Disease effects and coexpression directions", loc="left")
    rho = summary["spearman_meta_g_vs_adjusted_meta_r"]["rho"]
    p = summary["spearman_meta_g_vs_adjusted_meta_r"]["p"]
    frac = summary["direction_concordance"]["fraction"]
    ax.text(.03, .97, f"Spearman ρ={rho:.2f}; P={p:.1e}\nDirection concordance={frac:.1%}",
            transform=ax.transAxes, va="top", fontsize=5.8)


def plot_overlap(ax, overlap: pd.DataFrame) -> None:
    z = overlap.copy()
    labels = z.module
    observed = z.overlap_n.astype(float)
    expected = z.expected_overlap.astype(float)
    y = np.arange(len(z))[::-1]
    ax.barh(y, expected, color=NEUTRAL, height=.56, label="Expected")
    ax.scatter(observed, y, s=35, color=POSITIVE, edgecolor="white", linewidth=.5, label="Observed")
    ax.set_yticks(y, labels)
    ax.set_xlabel("Overlapping genes")
    ax.set_title("No gene-set membership enrichment", loc="left")
    for yy, row in zip(y, z.itertuples()):
        ax.text(max(row.overlap_n, row.expected_overlap)+.12, yy,
                f"P={row.p_enrichment:.3g}", va="center", fontsize=5.5)
    ax.legend(loc="lower right", fontsize=5.7)


def plot_theme_bridge(ax, theme: pd.DataFrame) -> None:
    pivot = theme.pivot(index="theme", columns="query_set", values="best_fdr")
    order = ["Cilium/projection", "ECM/adhesion", "Muscle/cytoskeleton", "Vascular signaling", "Stress/protein"]
    columns = ["robust_positive", "core_positive", "robust_negative", "core_negative"]
    pivot = pivot.reindex(index=order, columns=columns)
    matrix = -np.log10(pivot.astype(float))
    masked = np.ma.masked_invalid(matrix.to_numpy())
    im = ax.imshow(masked, cmap="cividis", norm=Normalize(vmin=0, vmax=max(4, np.nanmax(matrix.to_numpy()))), aspect="auto")
    ax.set_yticks(range(len(order)), [x.replace("/", " / ") for x in order])
    ax.set_xticks(range(len(columns)), ["Robust +", "PDE8B +", "Robust −", "PDE8B −"], rotation=25, ha="right")
    ax.set_title("Descriptive biological-theme convergence", loc="left")
    for i in range(masked.shape[0]):
        for j in range(masked.shape[1]):
            text = "–" if masked.mask[i, j] else f"{masked[i, j]:.1f}"
            ax.text(j, i, text, ha="center", va="center", fontsize=5.6,
                    color="white" if (not masked.mask[i, j] and masked[i, j] > 2.2) else INK)
    return im


def plot_umap(ax, data: pd.DataFrame, field: str, palette: dict[str, str], title: str) -> None:
    for category in palette:
        z = data[data[field].astype(str).eq(str(category))]
        ax.scatter(z.UMAP1, z.UMAP2, s=1.25, color=palette[category], alpha=.60,
                   linewidths=0, rasterized=True, label=category)
    ax.set(xticks=[], yticks=[], xlabel="UMAP 1", ylabel="UMAP 2")
    ax.set_title(title, loc="left")
    ax.legend(loc="upper left", fontsize=5.2, markerscale=2.6, handletextpad=.2,
              labelspacing=.15, ncol=1)


def build_supplement(frames: dict[str, pd.DataFrame | dict]) -> None:
    bridge = frames["bridge"]; overlap = frames["overlap"]; theme = frames["theme"]
    summary = frames["bridge_summary"]; atlas = frames["atlas_umap"]; mural = frames["mural_umap"]
    assert isinstance(bridge, pd.DataFrame) and isinstance(overlap, pd.DataFrame)
    assert isinstance(theme, pd.DataFrame) and isinstance(summary, dict)
    assert isinstance(atlas, pd.DataFrame) and isinstance(mural, pd.DataFrame)

    fig = plt.figure(figsize=(7.35, 8.35))
    gs = fig.add_gridspec(3, 2, left=.09, right=.95, top=.94, bottom=.075,
                          height_ratios=[1.0, .95, 1.0], wspace=.42, hspace=.58)
    a = fig.add_subplot(gs[0, 0]); shift_axes(a, -11)
    b = fig.add_subplot(gs[0, 1]); shift_axes(b, 2)
    c = fig.add_subplot(gs[1, :]); shift_axes(c, 4)
    d = fig.add_subplot(gs[2, 0]); shift_axes(d, -17)
    e = fig.add_subplot(gs[2, 1]); shift_axes(e, -27)
    plot_bridge(a, bridge, summary); panel_label(a, "a", y=1.106, x_fig=-.038)
    plot_overlap(b, overlap); panel_label(b, "b", y=1.106, x_fig=.45)
    im = plot_theme_bridge(c, theme); panel_label(c, "c", y=1.106, x_fig=-.038)
    cb = fig.colorbar(im, ax=c, fraction=.018, pad=.018); cb.set_label("−log$_{10}$(FDR)", fontsize=6)
    cb_pos = cb.ax.get_position()
    cb.ax.set_position([
        cb_pos.x0 + 1 / (72 * fig.get_figwidth()),
        cb_pos.y0 - 3 / (72 * fig.get_figheight()),
        cb_pos.width,
        cb_pos.height,
    ])
    condition_palette = {"Donor": DONOR, "IPAH": PAH, "SSc-PAH": "#CC79A7"}
    plot_umap(d, atlas, "condition", condition_palette, "GSE293580 whole-atlas integration QC")
    panel_label(d, "d", x_fig=-.038)
    mural_palette = {str(x): c for x, c in zip(range(6), ["#6A51A3", "#9E9AC8", "#2B8CBE", SMC1, SMC2, "#8C6D5A"])}
    plot_umap(e, mural, "mural_leiden_0.4", mural_palette, "GSE293580 mural-cell refinement QC")
    panel_label(e, "e", x_fig=.45)
    save_bundle(fig, OUT / "Supplementary_Figure_PDE8B_program_support_review_v1")


def build_robustness_heatmaps() -> None:
    source = O / "pde8b-smc-integrated-review-v1"
    patient = pd.read_csv(source / "patient_smc_heatmap_source_data.csv", index_col=0)
    cohort = pd.read_csv(source / "independent_cohort_heatmap_source_data.csv", index_col=0)
    cmap = mpl.colors.LinearSegmentedColormap.from_list("f3_diverging", [DONOR, "#F7F7F7", PAH])

    fig, axes = plt.subplots(1, 2, figsize=(7.35, 4.15), gridspec_kw={"width_ratios":[1.45, 1.0]})
    fig.subplots_adjust(left=.10, right=.96, top=.82, bottom=.22, wspace=.42)
    a, b = axes
    im = a.imshow(patient, aspect="auto", cmap=cmap, vmin=-2, vmax=2, interpolation="nearest")
    a.set_yticks(range(len(patient)), patient.index, fontsize=5.7)
    a.set_xticks(range(len(patient.columns)), [x.replace(" | ", "\n") for x in patient.columns],
                 rotation=90, fontsize=5.0)
    a.set_title("Patient-level SMC pseudobulk consistency", loc="left")
    panel_label(a, "a")
    cb = fig.colorbar(im, ax=a, fraction=.035, pad=.025); cb.set_label("Gene-wise z score", fontsize=5.5)

    im2 = b.imshow(cohort, aspect="auto", cmap=cmap, vmin=-.75, vmax=.75, interpolation="nearest")
    b.set_yticks(range(len(cohort)), cohort.index, fontsize=5.7)
    b.set_xticks(range(len(cohort.columns)), cohort.columns, rotation=35, ha="right", fontsize=5.4)
    for i in range(len(cohort)):
        for j in range(len(cohort.columns)):
            val = cohort.iloc[i, j]
            b.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=5.0,
                   color="white" if abs(val) > .43 else INK)
    b.set_title("Four-cohort correlation consistency", loc="left")
    panel_label(b, "b")
    cb2 = fig.colorbar(im2, ax=b, fraction=.045, pad=.025); cb2.set_label("PDE8B correlation", fontsize=5.5)
    save_bundle(fig, OUT / "Supplementary_Figure_PDE8B_program_robustness_review_v1")


def write_qa_notes() -> None:
    notes = """# Figure 5 review-v1 QA notes

- Core conclusion: PDE8B reproducibly covaries with vascular/cyclic-nucleotide genes across four lung cohorts and has descriptive support in independent GSE293580 SMCs.
- Main figure statistical units: subjects within cohorts, cohort-level random-effects meta-analysis, and patient-level GSE293580 pseudobulk.
- Enrichment and bridge panels use genes/pathway sets rather than subjects.
- Removed from the revised composition: prior SMC-state composite already represented by tentative Figure 4; decorative 12-gene network; cell-level state trends duplicated in Figure 4; ligand-receptor proxies with no FDR-significant primary family.
- Historical source figures and raw/processed data were preserved; removal means exclusion from the revised composition, not filesystem deletion.
- GSE293580 strict-SMC replication is explicitly descriptive and includes all annotated subjects; individual patient points are shown.
- PDE8B is presented as an associated mechanism candidate, not a causal regulator or validated therapeutic target.
- Visual vocabulary inherits Figure 3: Donor blue, PAH orange-red, circular patient marks with a triangle for SSc-PAH subtype, white background, cividis for continuous heatmaps, and restrained neutral/core colors.
"""
    (OUT / "Figure_5_review_v1_QA_notes.md").write_text(notes, encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    frames = load_inputs()
    mode = sys.argv[1] if len(sys.argv) > 1 else "all"
    if mode not in {"all", "main", "supplement"}:
        raise ValueError("Mode must be one of: all, main, supplement")
    if mode in {"all", "main"}:
        build_main(frames)
    if mode in {"all", "supplement"}:
        build_supplement(frames)
        build_robustness_heatmaps()
    write_qa_notes()
    for path in sorted(OUT.iterdir()):
        print(path)


if __name__ == "__main__":
    main()
