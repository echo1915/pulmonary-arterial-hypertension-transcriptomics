"""Build a portable, non-destructive source-data and research-record package.

The package captures the exact figures used by ``build_manuscript_docx.py`` and
the smallest practical set of source tables, statistical outputs, scripts, and
raw-data locators needed to audit them. Public GEO matrices are not duplicated;
their locations are recorded because they remain immutable under ``raw/``.
"""

from __future__ import annotations

import csv
import hashlib
import json
import shutil
import sys
from collections import Counter, defaultdict
from datetime import date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.project_paths import PROJECT_DATA_ROOT, PROJECT_ROOT


STAMP = "2026-08-18"
PACKAGE_ROOT = PROJECT_DATA_ROOT / "outputs" / f"manuscript-research-record-package-{STAMP}"
OUTPUTS = PROJECT_DATA_ROOT / "outputs"

FIGURES = {
    "Figure_1": OUTPUTS / "figure1-publication-v4" / "Figure_1_overall_project_workflow_v5_draft",
    "Figure_2": OUTPUTS / "manuscript-figures-v3" / "Figure_2_gene_set_and_independent_validation_v6_draft",
    "Figure_3": OUTPUTS / "figure3-localization-review-v1" / "Figure_3_cell_atlas_markers_PDE8B_SMC_review_v3",
    "Figure_4": OUTPUTS / "figure4-smc-main-supp-review-v1" / "Figure_4_SMC_state_continuum_main_review_v1",
    "Figure_5": OUTPUTS / "figure5-program-review-v2" / "Figure_5_PDE8B_program_main_unified_v3",
    "Figure_6": OUTPUTS / "figure6-cyclic-nucleotide-program" / "Figure_6_PDE8B_cyclic_nucleotide_program",
    "Supplementary_Figure_1": OUTPUTS / "candidate-comparison" / "Supplementary_Figure_candidate_prioritization_bridge",
    "Supplementary_Figure_2": OUTPUTS / "manuscript-figures-v3" / "Supplementary_Figure_PDE8B_LOCO_v4_draft",
    "Supplementary_Figure_3": OUTPUTS / "figure4-smc-main-supp-review-v1" / "Supplementary_Figure_SMC_state_support_review_v1",
    "Supplementary_Figure_4": OUTPUTS / "figure5-program-review-v1" / "Supplementary_Figure_PDE8B_program_support_review_v1",
    "Supplementary_Figure_5": OUTPUTS / "figure5-program-review-v1" / "Supplementary_Figure_PDE8B_program_robustness_review_v1",
}

SOURCE_TABLES = {
    "Figure_1": [
        OUTPUTS / "figure1-publication-v4" / "Table_1_bulk_lung_cohort_overview.csv",
        OUTPUTS / "current-results" / "pah_audit" / "core_sample_manifest.csv",
    ],
    "Figure_2": [
        OUTPUTS / "current-results" / "pah_audit" / "lung_mechanism_random_effects_meta.csv",
        OUTPUTS / "current-results" / "pah_audit" / "gse53408_raw_cel_candidate_validation.csv",
        OUTPUTS / "current-results" / "pah_audit" / "pilot_all_common_gene_effects.csv",
    ],
    "Figure_3": [
        OUTPUTS / "current-results" / "pah_audit" / "scrna_gse210248_candidate_pseudobulk.csv",
        OUTPUTS / "current-results" / "pah_audit" / "scrna_gse210248_candidate_celltype_effects.csv",
        OUTPUTS / "scrna-patient-pseudobulk-audit" / "gse210248_pde8b_smc_patient_values.csv",
    ],
    "Figure_4": [
        OUTPUTS / "pde8b-coexpression" / "state-deconvolution" / "gse210248_smc_marker_pseudobulk.csv",
        OUTPUTS / "pde8b-coexpression" / "state-deconvolution" / "gse210248_smc_state_marker_effects.csv",
        OUTPUTS / "pde8b-coexpression" / "state-deconvolution" / "gse210248_celltype_proportions.csv",
        OUTPUTS / "pde8b-coexpression" / "state-deconvolution" / "gse210248_celltype_proportion_effects.csv",
    ],
    "Figure_5": [
        OUTPUTS / "figure5-program-review-v2" / "Figure_5_panel_b_selected_pathways_source_data.csv",
        OUTPUTS / "figure5-program-review-v2" / "Figure_5_panel_c_PAH_gene_connections_source_data.csv",
        OUTPUTS / "figure5-program-review-v2" / "Figure_5_panel_d_GSE293580_patient_values_source_data.csv",
        OUTPUTS / "figure5-program-review-v2" / "Figure_5_panel_e_core_stability_source_data.csv",
        OUTPUTS / "figure5-program-review-v2" / "Figure_5_panel_f_celltype_localization_source_data.csv",
        OUTPUTS / "pde8b-coexpression" / "pde8b_high_confidence_core_module.csv",
    ],
    "Figure_6": [
        OUTPUTS / "figure6-cyclic-nucleotide-program" / "Figure_6_panels_a-c_source_data.csv",
        OUTPUTS / "figure6-cyclic-nucleotide-program" / "Figure_6_panel_d_source_data.csv",
    ],
    "Supplementary_Figure_1": [
        OUTPUTS / "candidate-comparison" / "three_candidate_evidence_matrix.csv",
        OUTPUTS / "candidate-comparison" / "three_candidate_comparison_summary.json",
    ],
    "Supplementary_Figure_2": [
        OUTPUTS / "bulk-meta-robustness" / "lung_meta_loco_summary.csv",
        OUTPUTS / "bulk-meta-robustness" / "lung_meta_threshold_grid.csv",
    ],
    "Supplementary_Figure_3": [
        OUTPUTS / "pde8b-coexpression" / "state-deconvolution" / "gse210248_smc_marker_pseudobulk.csv",
        OUTPUTS / "pde8b-coexpression" / "state-deconvolution" / "gse210248_smc_state_marker_effects.csv",
    ],
    "Supplementary_Figure_4": [
        OUTPUTS / "pde8b-coexpression" / "state-deconvolution" / "pde8b_classic_pah_gene_connections.csv",
        OUTPUTS / "pde8b-coexpression" / "single-cell-module" / "gse210248_pde8b_module_scores.csv",
        OUTPUTS / "pde8b-coexpression" / "single-cell-module" / "gse293580_pde8b_module_scores.csv",
    ],
    "Supplementary_Figure_5": [
        OUTPUTS / "pde8b-coexpression" / "pde8b_core_loco_stability_summary.csv",
        OUTPUTS / "pde8b-coexpression" / "pde8b_module_leave_one_cohort_out.csv",
    ],
}

STATISTICAL_FILES = [
    OUTPUTS / "current-results" / "pah_audit" / "formal_random_effects_meta.csv",
    OUTPUTS / "current-results" / "pah_audit" / "formal_meta_summary.json",
    OUTPUTS / "current-results" / "pah_audit" / "lung_mechanism_meta_summary.json",
    OUTPUTS / "bulk-meta-robustness" / "bulk_meta_robustness_summary.json",
    OUTPUTS / "pde8b-coexpression" / "pde8b_cross_cohort_coexpression_meta.csv",
    OUTPUTS / "pde8b-coexpression" / "pde8b_per_cohort_correlations.csv",
    OUTPUTS / "gse293580-reanalysis" / "gse293580_pde8b_patient_level_tests.csv",
    OUTPUTS / "gse293580-reanalysis" / "gse293580_pde8b_cellstate_paired_tests.csv",
    OUTPUTS / "gse293580-reanalysis" / "gse293580_pde8b_disease_by_cellstate_interaction.csv",
    OUTPUTS / "gse293580-reanalysis" / "gse293580_smc_program_group_tests.csv",
    OUTPUTS / "r-reanalysis-2026-08-16" / "sessionInfo.txt",
    OUTPUTS / "r-reanalysis-2026-08-16" / "tables" / "gse210248_candidate_effects_R.csv",
    OUTPUTS / "r-reanalysis-2026-08-16" / "tables" / "gse293580_strict_SMC_tests_R.csv",
]

SCRIPTS = [
    "scripts/pah/lung_mechanism_meta.py",
    "scripts/pah/formal_meta_and_confounding.py",
    "scripts/pah/reprocess_gse53408_raw.R",
    "scripts/pah/scrna_candidate_localization.py",
    "scripts/pah/gse293580_pde8b_patient_level.py",
    "scripts/pah/pde8b_cross_cohort_coexpression.py",
    "scripts/pah/pde8b_pathway_enrichment.py",
    "scripts/pah/reanalyse_scrna_with_pahscrna.R",
    "scripts/make_figure1_publication_v3.py",
    "scripts/make_manuscript_figure2_v3.py",
    "scripts/assemble_figure3_localization_review.py",
    "scripts/assemble_figure4_smc_main_supp_review.py",
    "scripts/plot_figure5_program_review_v2.py",
    "scripts/plot_figure5_supplement_review_v2.py",
    "scripts/plot_figure6_cyclic_nucleotide_program.py",
    "scripts/compare_three_pah_candidates.py",
    "scripts/project_paths.py",
]

ACCESSIONS = ["GSE15197", "GSE113439", "GSE254617", "GSE208592", "GSE53408", "GSE210248", "GSE293580"]


def copy_file(src: Path, dst: Path) -> None:
    if not src.exists():
        raise FileNotFoundError(src)
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def build_n_audit() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = [
        {"figure": "Figure 1", "panel": "workflow/Table 1", "n_definition": "independent subjects", "groups": "discovery PAH=129; control=94; validation PAH=12; control=11", "status": "verified from manuscript and cohort table", "raw_image_rule": "not applicable; public transcriptomic subjects, not tissue images"},
        {"figure": "Figure 2", "panel": "a-c", "n_definition": "subjects within 4 discovery cohorts; cohorts are meta-analytic units", "groups": "GSE15197 18/13; GSE113439 14/11; GSE254617 82/52; GSE208592 15/18", "status": "verified", "raw_image_rule": "one source-data record per subject/effect; no microscopy image"},
        {"figure": "Figure 2", "panel": "d", "n_definition": "subjects", "groups": "GSE53408 PAH=12; control=11", "status": "verified", "raw_image_rule": "one expression value per subject"},
        {"figure": "Figure 3", "panel": "a-c", "n_definition": "cells shown descriptively", "groups": "GSE210248 donor patients=3; PAH patients=3", "status": "verified; cells are not biological replicates", "raw_image_rule": "UMAP/dot plot is computed from count matrix"},
        {"figure": "Figure 3", "panel": "d-e", "n_definition": "patient pseudobulk", "groups": "donor=3; PAH=3 per eligible SMC state", "status": "verified", "raw_image_rule": "one plotted value per patient"},
        {"figure": "Figure 4", "panel": "b/e", "n_definition": "patient pseudobulk / paired patient-state contrasts", "groups": "donor=3; PAH=3 where eligible", "status": "verified at subject level", "raw_image_rule": "cell-level curves descriptive; patient rows carry inference"},
        {"figure": "Figure 5", "panel": "a-c/e", "n_definition": "subjects within cohort; genes are features", "groups": "4 independent bulk cohorts", "status": "verified", "raw_image_rule": "genes/pathways are not biological n"},
        {"figure": "Figure 5", "panel": "d", "n_definition": "patient pseudobulk", "groups": "strict SMC donor=2; IPAH=3; SSc-PAH=4", "status": "verified; donor comparison descriptive", "raw_image_rule": "one plotted value per patient"},
        {"figure": "Figure 5", "panel": "f", "n_definition": "patient-cell-type pseudobulk", "groups": "GSE210248 donor=3; PAH=3 subject pool", "status": "verified", "raw_image_rule": "one value per eligible patient-cell-type combination"},
        {"figure": "Figure 6", "panel": "a-d", "n_definition": "subjects within 4 cohorts; genes/pathways are features", "groups": "4 discovery cohorts", "status": "verified", "raw_image_rule": "network/pathway nodes are not biological n"},
    ]
    return rows


def build_subject_traceability() -> list[dict[str, object]]:
    trace: list[dict[str, object]] = []
    sources = [
        ("Figure_3_d-e", OUTPUTS / "scrna-patient-pseudobulk-audit" / "gse210248_pde8b_smc_patient_values.csv", "sample", "disease", "cell_type", "n_cells"),
        ("Figure_5_d", OUTPUTS / "figure5-program-review-v2" / "Figure_5_panel_d_GSE293580_patient_values_source_data.csv", "sample", "condition", "compartment", "n_cells"),
    ]
    for figure_panel, path, sample_col, group_col, stratum_col, count_col in sources:
        for row in csv_rows(path):
            trace.append({
                "figure_panel": figure_panel,
                "dataset": "GSE210248" if "210248" in path.name else "GSE293580",
                "subject_id": row.get(sample_col, ""),
                "group": row.get(group_col, ""),
                "stratum": row.get(stratum_col, ""),
                "cells_contributing": row.get(count_col, ""),
                "source_file": str(path),
                "biological_n": 1,
                "raw_image_available": "no; public count matrix, not microscopy",
            })
    return trace


def raw_locator_rows() -> list[dict[str, object]]:
    raw_root = PROJECT_DATA_ROOT / "raw"
    found: list[dict[str, object]] = []
    for path in raw_root.rglob("*"):
        if not path.is_file():
            continue
        upper = path.name.upper()
        hits = [acc for acc in ACCESSIONS if acc in upper or acc in str(path.parent).upper()]
        for acc in hits:
            found.append({
                "accession": acc,
                "absolute_path": str(path),
                "relative_to_data_root": str(path.relative_to(PROJECT_DATA_ROOT)),
                "bytes": path.stat().st_size,
                "raw_read_only": "yes",
                "duplicated_in_package": "no",
            })
    return found


def research_record_text() -> str:
    return f"""# PAH–PDE8B 转录组稿件科研记录（可誊写版）

日期：{STAMP}

项目编号：project-001

标准项目名：project-001-pulmonary-arterial-hypertension-transcriptomics

代码目录：{PROJECT_ROOT}

外置数据目录：{PROJECT_DATA_ROOT}

## 一、研究目的

整合多个独立人肺转录组队列，筛选在肺动脉高压（PAH）中方向一致、效应量稳定且异质性较低的基因；随后利用独立 bulk 队列及两个单细胞 RNA 测序队列验证候选基因，并判断 PDE8B 是否定位于血管平滑肌细胞及其是否与环核苷酸相关重塑程序共同变化。

## 二、研究设计

主发现队列为 GSE15197、GSE113439、GSE254617 和 GSE208592，共 223 名独立受试者（PAH 129，control 94）。GSE53408 为独立 raw-CEL 验证队列，共 23 名受试者（PAH 12，control 11）。GSE210248 用于肺动脉单细胞定位，共 3 名供体和 3 名 PAH 患者；GSE293580 用于 IPAH 与 SSc-PAH 的独立病因学复现。GSE117261 与 GSE254617 中 PHBI 受试者存在重叠，未作为独立主队列重复计数。

## 三、数据处理与统计单位

各 bulk 队列分别预处理，不跨平台直接拼接矩阵。每个基因在每个队列内计算 PAH 与对照间 Hedges' g，再进行随机效应 Meta 分析。受试者是 bulk 与单细胞伪批量分析的生物学重复，队列是 Meta 分析的贡献单位。单个细胞、基因和通路条目均不作为独立生物学重复。多重检验采用 Benjamini–Hochberg FDR 校正；检验均为双侧，除非另有说明。

单细胞比较先按患者和细胞类型聚合为伪批量值。GSE210248 中每个患者–细胞类型组合少于 5 个细胞时不进入展示比较。UMAP、细胞级表达图和平滑曲线仅用于描述，不用于把细胞当作 n 进行推断。

## 四、主要分析与结果

1. 四队列 Meta 分析在共同基因中筛得 188 个稳健基因，标准为四队列方向一致、FDR<0.05、|合并 Hedges' g|≥0.50 且 I²<50%。
2. PDE8B 在四个发现队列中均升高，合并 Hedges' g=0.744，FDR=8.81×10⁻⁵，I²=0%。留一队列分析的合并 g 范围为 0.680–1.009，方向始终为正。
3. 独立 GSE53408 raw-CEL 验证中，PDE8B 仍升高，Hedges' g=1.077，FDR=0.018。PIEZO2 和 SLC16A12 也通过 bulk 验证，但跨细胞层和病因层证据不如 PDE8B 连续。
4. GSE210248 的患者级伪批量结果显示，PDE8B 在 SMC1 与 SMC2 中 PAH–donor 效应均为正（g 约 0.73 和 0.74）。该数据只有 3 名 donor 与 3 名 PAH，因此用于定位和方向支持，不作为强因果证据。
5. PAH SMC 中收缩程序降低、ECM/重塑程序升高，PDE8B 同时在 SMC1 和 SMC2 中升高。该结果支持 PDE8B 与 SMC 状态重塑相关，但不能证明单个细胞发生 SMC1→SMC2 谱系转换，也不能证明 PDE8B 驱动转换。
6. PDE8B 跨四个 bulk 队列的组别校正共表达分析得到 543 个稳定机制模块基因，其中 201 个达到高置信核心阈值。PDE5A、PRKG1、PDE3B 等与 PDE8B 稳定正相关；cGMP–PKG 通路富集 FDR=0.0275。
7. PRKACA 与 PDE8B 无稳定相关，CREB1 的合并相关未达显著。因此结果应表述为“选择性的环核苷酸相关转录程序”，不能表述为经典 cAMP–PKA–CREB 全轴普遍激活。

## 五、图件 n 值及“原始图片”说明

本研究为公开转录组数据再分析，没有每名受试者对应的组织切片、Western blot 或免疫荧光原始图片。稿件中 n 主要指独立受试者；每个 n 对应源表达矩阵中的一名受试者及一个患者级数据点，而不是一张实验图片。细胞数只表示该患者伪批量值由多少细胞贡献，不等于生物学 n。逐图 n 值见 `04_n值与逐样本追溯/n_value_audit.csv`，逐患者来源见 `subject_level_traceability.csv`。

## 六、质量控制与敏感性分析

已避免将 GSE117261 与 GSE254617 重复计数；GSE254617 技术列先聚合到唯一 PHBI 受试者；GSE53408 因公开处理矩阵经过过滤，仅作为 raw-CEL 重处理后的次级验证；Meta 结果执行留一队列分析；单细胞推断使用患者级伪批量，避免细胞级伪重复；GSE293580 strict-SMC 中 donor 仅 n=2，因此相应比较只作描述。

## 七、当前结论

PDE8B 是 PAH 肺组织中跨队列稳定升高、定位于血管 SMC 并与 SMC 重塑及选择性环核苷酸程序相关的机制候选基因。现有证据属于关联、定位和跨队列复现，不足以证明 PDE8B 是 PAH 的致病因子、直接调控 cGMP，或已构成可用治疗靶点。

## 八、后续实验建议

在原代人肺动脉平滑肌细胞中进行 PDE8B 敲低/抑制与过表达，检测区室化 cAMP、cGMP、PKA/CREB、PKG、增殖、迁移及 ACTA2/MYH11、COL3A1/FN1 等标志；在人 PAH 组织中进行 PDE8B 与收缩型及 ECM 型 SMC 标志的共定位；实验按 IPAH 与 SSc-PAH 分层，并以独立供体作为 n。

## 九、归档说明

本归档包不覆盖原始数据和既有结果。`raw/` 保持只读，公共 GEO 大型矩阵不重复复制，仅在 `01_公共原始数据定位` 中记录位置。归档内所有复制文件均列入 SHA-256 清单，以便后续核验。
"""


def main() -> None:
    if PACKAGE_ROOT.exists():
        raise FileExistsError(f"Refusing to overwrite existing package: {PACKAGE_ROOT}")
    PACKAGE_ROOT.mkdir(parents=True)

    copy_file(PROJECT_ROOT / "manuscript" / "PAH_PDE8B_manuscript_v0.1.md", PACKAGE_ROOT / "00_稿件快照" / "PAH_PDE8B_manuscript_v0.1.md")
    for docx in sorted((PROJECT_ROOT / "manuscript" / "word").glob("PAH_PDE8B_manuscript*.docx"), key=lambda p: p.stat().st_mtime, reverse=True)[:2]:
        copy_file(docx, PACKAGE_ROOT / "00_稿件快照" / docx.name)

    for label, stem in FIGURES.items():
        for suffix in (".png", ".pdf", ".svg", ".tiff"):
            src = stem.with_suffix(suffix)
            if src.exists():
                copy_file(src, PACKAGE_ROOT / "02_稿件图片" / label / src.name)

    for label, paths in SOURCE_TABLES.items():
        for src in paths:
            copy_file(src, PACKAGE_ROOT / "03_逐图源数据_Source_Data" / label / src.name)

    for src in STATISTICAL_FILES:
        copy_file(src, PACKAGE_ROOT / "05_统计学分析" / "results" / src.name)
    for relative in SCRIPTS:
        src = PROJECT_ROOT / relative
        copy_file(src, PACKAGE_ROOT / "05_统计学分析" / "scripts" / relative.replace("/", "__"))

    raw_rows = raw_locator_rows()
    write_csv(PACKAGE_ROOT / "01_公共原始数据定位" / "raw_data_locator.csv", raw_rows, ["accession", "absolute_path", "relative_to_data_root", "bytes", "raw_read_only", "duplicated_in_package"])
    manifests = PROJECT_DATA_ROOT / "manifests"
    if manifests.exists():
        for src in manifests.glob("*"):
            if src.is_file():
                copy_file(src, PACKAGE_ROOT / "01_公共原始数据定位" / "existing_manifests" / src.name)

    n_rows = build_n_audit()
    write_csv(PACKAGE_ROOT / "04_n值与逐样本追溯" / "n_value_audit.csv", n_rows, ["figure", "panel", "n_definition", "groups", "status", "raw_image_rule"])
    trace_rows = build_subject_traceability()
    write_csv(PACKAGE_ROOT / "04_n值与逐样本追溯" / "subject_level_traceability.csv", trace_rows, ["figure_panel", "dataset", "subject_id", "group", "stratum", "cells_contributing", "source_file", "biological_n", "raw_image_available"])

    record = research_record_text()
    (PACKAGE_ROOT / "06_科研记录").mkdir(parents=True, exist_ok=True)
    (PACKAGE_ROOT / "06_科研记录" / "科研记录_可誊写版.md").write_text(record, encoding="utf-8")

    readme = f"""# 稿件原始数据、图片与统计归档包

- 创建日期：{STAMP}
- 项目编号：project-001
- 标准项目名：project-001-pulmonary-arterial-hypertension-transcriptomics
- 代码目录：{PROJECT_ROOT}
- 外置数据目录：{PROJECT_DATA_ROOT}
- 当前稿件图：6 张主图、5 张补充图

推荐阅读顺序：先读 `06_科研记录/科研记录_可誊写版.md`，再查 `04_n值与逐样本追溯/n_value_audit.csv`。大型 GEO 原始矩阵未重复复制；其只读位置记录在 `01_公共原始数据定位/raw_data_locator.csv`。本项目不存在逐受试者实验显微原图，n 对应独立受试者或队列，详细解释见科研记录。
"""
    (PACKAGE_ROOT / "README.md").write_text(readme, encoding="utf-8")

    manifest_rows = []
    for path in sorted(PACKAGE_ROOT.rglob("*")):
        if path.is_file() and path.name != "package_sha256_manifest.csv":
            manifest_rows.append({"relative_path": str(path.relative_to(PACKAGE_ROOT)), "bytes": path.stat().st_size, "sha256": sha256(path)})
    write_csv(PACKAGE_ROOT / "package_sha256_manifest.csv", manifest_rows, ["relative_path", "bytes", "sha256"])

    summary = {
        "package_root": str(PACKAGE_ROOT),
        "files": len(manifest_rows) + 1,
        "bytes_excluding_manifest": sum(int(row["bytes"]) for row in manifest_rows),
        "figures": len(FIGURES),
        "source_tables": sum(len(v) for v in SOURCE_TABLES.values()),
        "raw_locator_rows": len(raw_rows),
        "subject_traceability_rows": len(trace_rows),
        "created": str(date.today()),
    }
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
