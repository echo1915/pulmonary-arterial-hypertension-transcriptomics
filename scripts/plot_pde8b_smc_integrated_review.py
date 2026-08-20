from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize, TwoSlopeNorm
import numpy as np
import pandas as pd

from project_paths import PROJECT_DATA_ROOT


ROOT = PROJECT_DATA_ROOT / "outputs"
COEX = ROOT / "pde8b-coexpression"
TREND = ROOT / "smc-gene-state-trends-review-v1" / "Figure_SMC_gene_state_trends_source_data.csv"
OUT = ROOT / "pde8b-smc-integrated-review-v1"
COHORTS = ["GSE15197", "GSE113439", "GSE254617", "GSE208592"]
SMC_COLORS = {"SMC 1": "#E86B6B", "SMC 2": "#42B883"}
BLUE, PURPLE, RED, INK = "#4169A1", "#7652A6", "#F21A24", "#263746"
CMAP = LinearSegmentedColormap.from_list("bwr_soft", [BLUE, "#F7F7F7", RED])

mpl.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7, "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 0.75, "pdf.fonttype": 42, "svg.fonttype": "none",
})


def panel(ax, label: str) -> None:
    ax.text(-0.08, 1.03, label, transform=ax.transAxes, fontsize=11,
            fontweight="bold", va="bottom", ha="left")


def smooth(x: np.ndarray, y: np.ndarray, bins: int = 24) -> tuple[np.ndarray, np.ndarray]:
    edges = np.unique(np.quantile(x, np.linspace(0, 1, bins + 1)))
    xx, yy = [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (x >= lo) & (x <= hi if hi == edges[-1] else x < hi)
        if m.sum() >= 10:
            xx.append(np.median(x[m])); yy.append(np.mean(y[m]))
    return np.asarray(xx), np.asarray(yy)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    core = pd.read_csv(COEX / "pde8b_high_confidence_core_module.csv")
    pb = pd.read_csv(COEX / "single-cell-module" / "gse210248_pde8b_core_gene_pseudobulk.csv")
    trend = pd.read_csv(TREND)

    smc_pb = pb[pb["cell_type"].isin(["SMC 1", "SMC 2"])].copy()
    detect = (smc_pb.groupby("gene")
              .agg(mean_expr=("log_cp10k", "mean"), variance=("log_cp10k", "var"),
                   detected=("detection_fraction", lambda v: (v > 0).mean()))
              .reset_index())
    candidates = core.merge(detect, left_on="symbol", right_on="gene", how="inner")
    candidates = candidates[(candidates["detected"] >= 0.35) & (candidates["variance"] > 0)]
    candidates["rank_score"] = candidates["adjusted_meta_r"].abs() * np.sqrt(candidates["variance"])
    pos = candidates[candidates["adjusted_meta_r"] > 0].nlargest(6, "rank_score")
    neg = candidates[candidates["adjusted_meta_r"] < 0].nlargest(6, "rank_score")
    selected = pd.concat([pos, neg]).drop_duplicates("symbol").head(12).copy()
    genes = selected["symbol"].tolist()

    # Patient-by-state heatmap: gene-wise z scores across 12 biological pseudobulks.
    hp = smc_pb[smc_pb["gene"].isin(genes)].copy()
    hp["column"] = hp["sample"] + " | " + hp["cell_type"].str.replace(" ", "")
    mat = hp.pivot_table(index="gene", columns="column", values="log_cp10k")
    sample_order = [f"{s} | {ct.replace(' ', '')}" for s in
                    ["Donor_1", "Donor_2", "Donor_3", "PAH_1", "PAH_2", "PAH_3"]
                    for ct in ["SMC 1", "SMC 2"]]
    mat = mat.reindex(index=genes, columns=sample_order).fillna(0)
    mat_z = mat.apply(lambda r: (r-r.mean())/r.std(ddof=1) if r.std(ddof=1)>0 else 0, axis=1)

    # Independent cohort heatmap: adjusted PDE8B-gene correlations.
    corr_cols = [f"r_{c}" for c in COHORTS]
    cm = selected.set_index("symbol")[corr_cols].copy()
    cm.columns = COHORTS

    fig = plt.figure(figsize=(11.2, 8.2))
    gs = fig.add_gridspec(3, 6, height_ratios=[1.22, 1.05, 1.0],
                          left=.065, right=.96, top=.96, bottom=.08,
                          hspace=.54, wspace=.72)

    # a: restrained star network, co-expression only.
    ax = fig.add_subplot(gs[0, :3]); panel(ax, "a")
    theta = np.linspace(0, 2*np.pi, len(selected), endpoint=False) + np.pi/12
    radius = 1.0
    norm = TwoSlopeNorm(vmin=-.7, vcenter=0, vmax=.7)
    for angle, (_, row) in zip(theta, selected.iterrows()):
        x, y = radius*np.cos(angle), .78*radius*np.sin(angle)
        r = row["adjusted_meta_r"]
        ax.plot([0, x], [0, y], color=CMAP(norm(r)), lw=.7+2.1*abs(r), alpha=.72, zorder=1)
        ax.scatter(x, y, s=42+90*abs(r), color=CMAP(norm(r)), edgecolor="white", lw=.7, zorder=3)
        ax.text(x*1.14, y*1.14, row["symbol"], fontsize=6.5, ha="center", va="center")
    ax.scatter(0, 0, s=255, color=PURPLE, edgecolor="white", lw=1.2, zorder=4)
    ax.text(0, 0, "PDE8B", color="white", fontsize=7.5, fontweight="bold", ha="center", va="center")
    ax.set_xlim(-1.42, 1.42); ax.set_ylim(-1.12, 1.12); ax.axis("off")
    ax.text(.5, -.05, "Stable cross-cohort co-expression", transform=ax.transAxes,
            ha="center", color="#59666F", fontsize=6.8)

    # b: patient-level SMC heatmap.
    ax = fig.add_subplot(gs[0, 3:]); panel(ax, "b")
    im = ax.imshow(mat_z, aspect="auto", cmap=CMAP, norm=Normalize(-2, 2), interpolation="nearest")
    ax.set_yticks(np.arange(len(genes)), genes, fontsize=6.3)
    ax.set_xticks(np.arange(len(sample_order)), [x.replace(" | ", "\n") for x in sample_order],
                  rotation=90, fontsize=5.6)
    for x in [1.5, 3.5, 5.5, 7.5, 9.5]: ax.axvline(x, color="white", lw=1.5)
    ax.set_xlabel("Patient-level SMC pseudobulk")
    cb = fig.colorbar(im, ax=ax, fraction=.035, pad=.025); cb.set_label("Gene-wise z score", fontsize=6.5)

    # c: independent bulk cohorts.
    ax = fig.add_subplot(gs[1, 1:5]); panel(ax, "c")
    im2 = ax.imshow(cm, aspect="auto", cmap=CMAP, norm=Normalize(-.75, .75), interpolation="nearest")
    ax.set_yticks(np.arange(len(genes)), genes, fontsize=6.4)
    ax.set_xticks(np.arange(len(COHORTS)), COHORTS, fontsize=6.5)
    for i in range(len(genes)):
        for j in range(len(COHORTS)):
            ax.text(j, i, f"{cm.iloc[i,j]:.2f}", ha="center", va="center", fontsize=5.5,
                    color="white" if abs(cm.iloc[i,j]) > .42 else INK)
    cb2 = fig.colorbar(im2, ax=ax, fraction=.025, pad=.025); cb2.set_label("PDE8B correlation", fontsize=6.5)

    # d-f: single-cell visualization along a defined state score.
    for label, col, gene in zip(["d", "e", "f"], range(0, 6, 2), ["MYH11", "FN1", "PDE8B"]):
        ax = fig.add_subplot(gs[2, col:col+2]); panel(ax, label)
        for ct in ["SMC 1", "SMC 2"]:
            q = trend[trend["cell_type"].eq(ct)]
            ax.scatter(q["state_score"], q[gene], s=4.2, alpha=.22, color=SMC_COLORS[ct],
                       linewidths=0, rasterized=True, label=ct)
        q = trend.sort_values("state_score")
        xx, yy = smooth(q["state_score"].to_numpy(), q[gene].to_numpy())
        ax.plot(xx, yy, color="#D64A61", lw=1.8)
        ax.set_title(gene, fontstyle="italic", pad=6)
        ax.set_xlabel("SMC state score\nECM − contractile")
        ax.set_ylabel("Expression (log CP10K)")
    handles, labs = ax.get_legend_handles_labels()
    ax.legend(handles, labs, frameon=False, loc="upper left", bbox_to_anchor=(1.01, 1.0), fontsize=6.5)

    stem = OUT / "Figure_PDE8B_SMC_integrated_review_v1"
    fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".tiff"), dpi=600, bbox_inches="tight", facecolor="white",
                pil_kwargs={"compression": "tiff_lzw"})
    selected.to_csv(OUT / "selected_high_confidence_genes.csv", index=False)
    mat_z.to_csv(OUT / "patient_smc_heatmap_source_data.csv")
    cm.to_csv(OUT / "independent_cohort_heatmap_source_data.csv")
    print(stem.with_suffix(".png"))


if __name__ == "__main__":
    main()
