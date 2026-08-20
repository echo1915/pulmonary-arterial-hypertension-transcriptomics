"""Build the canonical artifact for the integrated PDE8B-SMC evidence report."""

from __future__ import annotations

import json
import sqlite3
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SCRIPT_DIR))
from project_paths import PROJECT_DATA_ROOT, PROJECT_ROOT  # noqa: E402


OUTPUTS = PROJECT_DATA_ROOT / "outputs"
REPORT_DIR = PROJECT_ROOT / "reports" / "integrated-evidence-report"


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def row(df: pd.DataFrame, **filters: object) -> pd.Series:
    keep = pd.Series(True, index=df.index)
    for column, value in filters.items():
        keep &= df[column] == value
    selected = df.loc[keep]
    if len(selected) != 1:
        raise ValueError(f"expected one row for {filters}, found {len(selected)}")
    return selected.iloc[0]


def source(source_id: str, label: str, path: str, description: str) -> dict:
    return {
        "id": source_id,
        "label": label,
        "path": path,
        "query": {
            "engine": "local_reproducible_analysis",
            "language": "python",
            "description": description,
            "tables_used": [path],
            "filters": ["Patient/cohort-level statistical units; saved reviewed outputs only"],
        },
    }


def sql_literal(value: object) -> str:
    if value is None:
        return "NULL"
    if isinstance(value, bool):
        return "1" if value else "0"
    if isinstance(value, (int, float)):
        return str(value)
    return "'" + str(value).replace("'", "''") + "'"


def materialize_rows(
    rows: list[dict[str, object]], table_name: str, order_by: str
) -> tuple[list[dict[str, object]], str]:
    columns = list(rows[0])
    values = ",\n    ".join(
        "(" + ", ".join(sql_literal(item[column]) for column in columns) + ")"
        for item in rows
    )
    query = (
        f"WITH {table_name}({', '.join(columns)}) AS (\n"
        f"  VALUES\n    {values}\n"
        ")\n"
        f"SELECT {', '.join(columns)} FROM {table_name} ORDER BY {order_by};"
    )
    connection = sqlite3.connect(":memory:")
    connection.row_factory = sqlite3.Row
    materialized = [dict(item) for item in connection.execute(query).fetchall()]
    connection.close()
    return materialized, query


def main() -> None:
    coexpr = read_json(OUTPUTS / "pde8b-coexpression" / "pde8b_coexpression_summary.json")
    sensitivity = read_json(
        OUTPUTS / "gse293580-reanalysis" / "sensitivity" / "gse293580_sensitivity_summary.json"
    )
    annotation = read_json(
        OUTPUTS / "gse293580-reanalysis" / "sensitivity" / "gse293580_annotation_boundary_summary.json"
    )
    state_tests = pd.read_csv(
        OUTPUTS / "gse293580-reanalysis" / "gse293580_smc_program_group_tests.csv"
    )
    correlations = pd.read_csv(
        OUTPUTS / "gse293580-reanalysis" / "gse293580_pde8b_smc_program_correlations.csv"
    )
    interaction = pd.read_csv(
        OUTPUTS / "gse293580-reanalysis" / "gse293580_pde8b_disease_by_cellstate_interaction.csv"
    )
    gse210_state = pd.read_csv(
        OUTPUTS / "smc-composition-within-state" / "gse210248_smc_composition_within_state_effects.csv"
    )
    lr_tests = pd.read_csv(
        OUTPUTS / "gse293580-reanalysis" / "ligand-receptor-proxy" / "gse293580_mural_lr_patient_tests.csv"
    )
    lr_cross = pd.read_csv(
        OUTPUTS
        / "gse210248-lr-smc-component-replication"
        / "gse210248_gse293580_lr_smc_component_concordance.csv"
    )

    contractile_293 = row(
        state_tests, contrast="PAH_combined_vs_Donor", module="contractile"
    )
    matrix_293 = row(
        state_tests, contrast="PAH_combined_vs_Donor", module="matrix_remodeling"
    )
    contractile_210_smc1 = row(
        gse210_state, layer="Within-state", cell_type="SMC 1", feature="Contractile"
    )
    contractile_210_smc2 = row(
        gse210_state, layer="Within-state", cell_type="SMC 2", feature="Contractile"
    )
    matrix_210_smc1 = row(
        gse210_state, layer="Within-state", cell_type="SMC 1", feature="Synthetic/ECM"
    )
    matrix_210_smc2 = row(
        gse210_state, layer="Within-state", cell_type="SMC 2", feature="Synthetic/ECM"
    )
    interaction_primary = row(interaction, contrast="PAH_combined_vs_Donor")
    ssc_ecm_forward = row(
        lr_tests,
        smc_boundary="core_SMC",
        direction="SMC_to_pericyte",
        axis="ECM_integrin",
        contrast="SSc-PAH_vs_Donor",
    )
    ssc_ecm_reverse = row(
        lr_tests,
        smc_boundary="core_SMC",
        direction="pericyte_to_SMC",
        axis="ECM_integrin",
        contrast="SSc-PAH_vs_Donor",
    )
    ligand_concordance = float(lr_cross.loc[lr_cross["smc_role"] == "ligand", "direction_concordant"].mean())
    receptor_concordance = float(lr_cross.loc[lr_cross["smc_role"] == "receptor", "direction_concordant"].mean())
    weak_corr = correlations.loc[correlations["scope"] == "PAH_only"]

    evidence = [
        {
            "grade_rank": 1,
            "grade": "A",
            "claim": "PDE8B-centered lung co-expression program is reproducible",
            "evidence_scope": "Four independent human lung bulk cohorts",
            "key_result": f"{coexpr['high_confidence_core_abs_r_ge_0.40']} high-confidence genes; {coexpr['stable_mechanism_module']} stable genes",
            "robustness": "Adjusted within-cohort correlations plus random-effects meta-analysis and leave-one-cohort-out stability",
            "failure_or_boundary": "Co-expression cannot identify PDE8B as the upstream regulator",
            "allowed_language": "Robust PDE8B-centered molecular association",
        },
        {
            "grade_rank": 2,
            "grade": "B",
            "claim": "PDE8B elevation localizes to SMC-associated states",
            "evidence_scope": "GSE210248 SMC1/SMC2 and GSE293580 mural-cell reanalysis",
            "key_result": f"GSE293580 strict-SMC g={sensitivity['primary_combined_PAH_SMC_PDE8B_hedges_g']:.2f}; core-SMC g={annotation['core_SMC_hedges_g']:.2f}",
            "robustness": "Positive at 3, 5, and 10-cell thresholds; 9/9 leave-one-patient directions at primary threshold",
            "failure_or_boundary": "Only two eligible GSE293580 donors; single-cell detection is sparse",
            "allowed_language": "SMC-associated PDE8B elevation",
        },
        {
            "grade_rank": 2,
            "grade": "B",
            "claim": "PAH SMCs shift away from contraction toward matrix remodeling",
            "evidence_scope": "Directionally replicated in GSE210248 and GSE293580",
            "key_result": (
                f"Contractile g: {contractile_210_smc1.hedges_g:.2f}/{contractile_210_smc2.hedges_g:.2f} in GSE210248, "
                f"{contractile_293.hedges_g:.2f} in GSE293580; matrix g={matrix_293.hedges_g:.2f}"
            ),
            "robustness": (
                f"GSE210248 ECM g={matrix_210_smc1.hedges_g:.2f}/{matrix_210_smc2.hedges_g:.2f}; "
                "GSE293580 direction stable at thresholds 5 and 10"
            ),
            "failure_or_boundary": "Module definitions and score scales differ; patient groups are small",
            "allowed_language": "Directionally replicated SMC remodeling state",
        },
        {
            "grade_rank": 2,
            "grade": "B",
            "claim": "Several SMC ligand programs replicate across cohorts",
            "evidence_scope": "GSE293580 core SMC versus GSE210248 SMC1/SMC2",
            "key_result": f"Ligand direction concordance {ligand_concordance:.1%}; receptor concordance {receptor_concordance:.1%}",
            "robustness": "PDGF, TGF-beta, JAG, ECM, and THBS1 ligand directions agree in both SMC subtypes",
            "failure_or_boundary": "GSE210248 has no pericyte annotation; this is component—not communication—replication",
            "allowed_language": "Cross-cohort SMC ligand-program concordance",
        },
        {
            "grade_rank": 3,
            "grade": "C",
            "claim": "SSc-PAH may preferentially engage ECM-integrin exchange",
            "evidence_scope": "GSE293580 only, four eligible SSc-PAH versus two donors",
            "key_result": f"SMC→pericyte g={ssc_ecm_forward.hedges_g:.2f}; pericyte→SMC g={ssc_ecm_reverse.hedges_g:.2f}; exact P=0.067 each",
            "robustness": "Both directions retain sign in 6/6 leave-one-patient iterations",
            "failure_or_boundary": "Not seen in IPAH; expression proxy lacks spatial/protein confirmation",
            "allowed_language": "SSc-PAH-prioritized ECM-integrin hypothesis",
        },
        {
            "grade_rank": 3,
            "grade": "C",
            "claim": "Disease-specific SMC-pericyte PDE8B localization shift is not established",
            "evidence_scope": "GSE293580 paired SMC and pericyte compartments",
            "key_result": f"Disease-by-cell-state g={interaction_primary.hedges_g:.2f}, exact P={interaction_primary.exact_permutation_p:.3f}",
            "robustness": "Effect weakens as the minimum-cell threshold rises",
            "failure_or_boundary": "No formal evidence for a stable interaction",
            "allowed_language": "Exploratory localization contrast only",
        },
        {
            "grade_rank": 3,
            "grade": "C",
            "claim": "Within-PAH PDE8B-state coupling is weak and patient-sensitive",
            "evidence_scope": "Seven GSE293580 PAH patients with eligible strict SMCs",
            "key_result": f"PAH-only |rho| range {weak_corr.spearman_rho.abs().min():.2f}–{weak_corr.spearman_rho.abs().max():.2f}; all q≥{weak_corr.BH_q_within_scope.min():.3f}",
            "robustness": "Leave-one-patient ranges cross zero for every module",
            "failure_or_boundary": "Does not support stable patient-level coupling",
            "allowed_language": "PDE8B elevation accompanies, but is not tightly coupled to, state scores",
        },
        {
            "grade_rank": 3,
            "grade": "C",
            "claim": "Ligand-receptor compatibility axes are candidate-ranking signals",
            "evidence_scope": "GSE293580 expression-based SMC-pericyte proxy",
            "key_result": "TGF-beta, ECM-integrin, PDGF, CXCL12-CXCR4, and JAG-NOTCH3 show directional effects",
            "robustness": "All 14 axis directions preserved after adding SMC_like",
            "failure_or_boundary": "All proxy-score q≥0.417; component directions can oppose the composite score",
            "allowed_language": "Expression compatibility proxy, not activated communication",
        },
        {
            "grade_rank": 4,
            "grade": "D",
            "claim": "PDE8B is a causal therapeutic target in PAH",
            "evidence_scope": "Not tested by the present observational transcriptomic analyses",
            "key_result": "No perturbation, spatial confirmation, protein activity, hemodynamic association, or outcome validation",
            "robustness": "Not applicable",
            "failure_or_boundary": "Causal claim is unsupported",
            "allowed_language": "Mechanism candidate requiring perturbational validation",
        },
    ]

    priorities = [
        {
            "priority": 1,
            "target": "PDE8B × SMC-state multiplex localization",
            "why_now": "Directly tests whether PDE8B-positive SMCs lose MYH11/ACTA2 and gain COL3A1/FN1",
            "preferred_validation": "Human lung or pulmonary-artery multiplex IF/RNAscope, stratified by IPAH and SSc-PAH",
            "decision_rule": "Advance only if patient-level colocalization reproduces the predicted state shift",
        },
        {
            "priority": 2,
            "target": "SSc-PAH ECM-integrin axis",
            "why_now": "Largest disease-subtype-specific communication proxy with bidirectional leave-one-out stability",
            "preferred_validation": "Spatial ligand/receptor colocalization plus integrin-pathway protein/activity readout",
            "decision_rule": "Treat as SSc-specific unless independently reproduced in IPAH",
        },
        {
            "priority": 3,
            "target": "TGF-beta SMC ligand program",
            "why_now": "Positive SMC ligand direction in both GSE293580 and both GSE210248 SMC subtypes",
            "preferred_validation": "Patient-stratified TGFB ligand and phospho-SMAD measurements",
            "decision_rule": "Require ligand and downstream activity to agree; expression alone is insufficient",
        },
        {
            "priority": 4,
            "target": "Larger patient-level replication cohort",
            "why_now": "Current single-cell comparisons are limited by two or three donors",
            "preferred_validation": "Independent IPAH and SSc-PAH cohort with subject-level SMC/pericyte summaries",
            "decision_rule": "Predefine cell thresholds and disease strata before testing",
        },
        {
            "priority": 5,
            "target": "PDE8B perturbation",
            "why_now": "Needed to move from association to mechanism",
            "preferred_validation": "SMC-focused knockdown/inhibition with cAMP, contractility, proliferation, and ECM readouts",
            "decision_rule": "Deprioritize if perturbation does not shift both cyclic-nucleotide and remodeling phenotypes",
        },
    ]

    cohorts = [
        {
            "evidence_layer": "Bulk lung co-expression",
            "cohorts": "GSE15197, GSE113439, GSE254617, GSE208592",
            "replicate_unit": "cohort subjects",
            "comparison": "Within-cohort residual rank correlation, random-effects meta-analysis",
            "principal_limit": "Association rather than cell localization or causality",
        },
        {
            "evidence_layer": "GSE210248 pulmonary-artery scRNA-seq",
            "cohorts": "3 PAH, 3 donors",
            "replicate_unit": "patient",
            "comparison": "SMC1/SMC2 patient pseudobulk",
            "principal_limit": "No independent pericyte annotation",
        },
        {
            "evidence_layer": "GSE293580 whole-lung scRNA-seq",
            "cohorts": "IPAH, SSc-PAH, donors",
            "replicate_unit": "patient/sample",
            "comparison": "Locked mural annotations; ≥5 cells per primary compartment",
            "principal_limit": "Only two eligible donors for most SMC/pericyte contrasts",
        },
        {
            "evidence_layer": "Causal/spatial/protein validation",
            "cohorts": "Not available",
            "replicate_unit": "Not applicable",
            "comparison": "Not tested",
            "principal_limit": "Blocks therapeutic-target and activated-communication claims",
        },
    ]
    evidence, evidence_sql = materialize_rows(evidence, "integrated_evidence", "grade_rank, claim")
    priorities, priority_sql = materialize_rows(priorities, "validation_priorities", "priority")
    cohorts, cohort_sql = materialize_rows(cohorts, "cohort_context", "evidence_layer")
    grade_count_values = [
        (
            grade,
            rank,
            sum(item["grade"] == grade for item in evidence),
            definition,
            primary_use,
        )
        for grade, rank, definition, primary_use in (
            ("A", 1, "Independent-cohort replication plus robustness", "Manuscript anchor"),
            ("B", 2, "Directional cross-dataset replication with small-n limits", "Supported context"),
            ("C", 3, "Single-cohort, heterogeneous, or non-robust evidence", "Hypothesis generation"),
            ("D", 4, "Untested causal or validation layer", "Explicit evidence gap"),
        )
    ]
    value_sql = ",\n    ".join(
        "(" + ", ".join(
            [
                "'" + str(value).replace("'", "''") + "'" if isinstance(value, str) else str(value)
                for value in values
            ]
        ) + ")"
        for values in grade_count_values
    )
    grade_count_sql = (
        "WITH evidence_grades(grade, grade_rank, claim_count, definition, primary_use) AS (\n"
        f"  VALUES\n    {value_sql}\n"
        ")\n"
        "SELECT grade, grade_rank, claim_count, definition, primary_use\n"
        "FROM evidence_grades\n"
        "ORDER BY grade_rank;"
    )
    connection = sqlite3.connect(":memory:")
    connection.row_factory = sqlite3.Row
    grade_counts = [dict(item) for item in connection.execute(grade_count_sql).fetchall()]
    connection.close()

    generated = datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")
    sources = [
        source(
            "integrated_evidence",
            "Integrated evidence build",
            "scripts/pah/build_integrated_evidence_report.py",
            "Combines reviewed saved outputs into the evidence matrix and validation priorities.",
        ),
        source(
            "coexpression_outputs",
            "PDE8B bulk co-expression outputs",
            "external-data-root/outputs/pde8b-coexpression/",
            "Four-cohort co-expression meta-analysis and leave-one-cohort-out stability outputs.",
        ),
        source(
            "single_cell_outputs",
            "Patient-level single-cell outputs",
            "external-data-root/outputs/gse293580-reanalysis/",
            "GSE293580 patient pseudobulk, state, sensitivity, and ligand-receptor proxy outputs.",
        ),
        source(
            "crosscohort_outputs",
            "GSE210248 cross-cohort outputs",
            "external-data-root/outputs/gse210248-lr-smc-component-replication/",
            "GSE210248 patient SMC component effects and cross-cohort direction comparisons.",
        ),
    ]
    sources.append({
        "id": "grade_distribution_sql",
        "label": "Evidence grade distribution query",
        "path": "reports/integrated-evidence-report/grade_distribution.sql",
        "query": {
            "engine": "sqlite",
            "language": "sql",
            "sql": grade_count_sql,
            "description": "Counts the nine reviewed claims by prespecified evidence grade.",
            "tables_used": ["evidence_grades"],
            "filters": ["All nine prespecified integrated evidence claims"],
            "metric_definitions": ["claim_count = number of evidence-matrix claims assigned to each grade"],
            "executed_at": generated,
        },
    })
    for source_id, label, filename, query, description, table_name in (
        (
            "evidence_matrix_sql",
            "Evidence matrix query",
            "evidence_matrix.sql",
            evidence_sql,
            "Materializes the nine reviewed evidence claims and their inference boundaries.",
            "integrated_evidence",
        ),
        (
            "cohort_context_sql",
            "Cohort context query",
            "cohort_context.sql",
            cohort_sql,
            "Materializes cohort, replicate-unit, comparison, and limitation definitions.",
            "cohort_context",
        ),
        (
            "validation_priorities_sql",
            "Validation priority query",
            "validation_priorities.sql",
            priority_sql,
            "Materializes evidence-backed validation priorities and decision rules.",
            "validation_priorities",
        ),
    ):
        sources.append({
            "id": source_id,
            "label": label,
            "path": f"reports/integrated-evidence-report/{filename}",
            "query": {
                "engine": "sqlite",
                "language": "sql",
                "sql": query,
                "description": description,
                "tables_used": [table_name],
                "filters": ["Reviewed integrated evidence snapshot"],
                "executed_at": generated,
            },
        })

    artifact = {
        "surface": "report",
        "manifest": {
            "version": 1,
            "surface": "report",
            "title": "PAH PDE8B–SMC evidence integration",
            "description": "Technical synthesis of bulk, single-cell, sensitivity, and ligand-receptor evidence.",
            "generatedAt": generated,
            "cards": [],
            "charts": [
                {
                    "id": "grade_distribution",
                    "title": "Claims by evidence grade",
                    "subtitle": "Nine prespecified claims; lower letter grades require more cautious language",
                    "type": "bar",
                    "dataset": "evidence_grade_counts",
                    "sourceId": "grade_distribution_sql",
                    "valueFormat": "number",
                    "encodings": {
                        "x": {"field": "grade", "type": "nominal", "label": "Evidence grade"},
                        "y": {"field": "claim_count", "type": "quantitative", "label": "Claims"},
                        "tooltip": [
                            {"field": "definition", "type": "nominal", "label": "Definition"},
                            {"field": "primary_use", "type": "nominal", "label": "Primary use"},
                        ],
                    },
                }
            ],
            "tables": [
                {
                    "id": "evidence_matrix",
                    "title": "Evidence matrix",
                    "subtitle": "Claims ranked by replication, robustness, and inference boundary",
                    "dataset": "evidence_matrix",
                    "sourceId": "evidence_matrix_sql",
                    "defaultSort": {"field": "grade_rank", "direction": "asc"},
                    "columns": [
                        {"field": "grade_rank", "label": "Rank", "type": "number"},
                        {"field": "grade", "label": "Grade", "type": "text"},
                        {"field": "claim", "label": "Claim", "type": "text"},
                        {"field": "evidence_scope", "label": "Evidence scope", "type": "text"},
                        {"field": "key_result", "label": "Key result", "type": "text"},
                        {"field": "robustness", "label": "Robustness", "type": "text"},
                        {"field": "failure_or_boundary", "label": "Failure / boundary", "type": "text"},
                        {"field": "allowed_language", "label": "Allowed language", "type": "text"},
                    ],
                },
                {
                    "id": "cohort_context",
                    "title": "Cohorts and comparison units",
                    "subtitle": "Biological replicate, comparison basis, and principal limitation by evidence layer",
                    "dataset": "cohort_context",
                    "sourceId": "cohort_context_sql",
                    "defaultSort": {"field": "evidence_layer", "direction": "asc"},
                    "columns": [
                        {"field": "evidence_layer", "label": "Evidence layer", "type": "text"},
                        {"field": "cohorts", "label": "Cohort / groups", "type": "text"},
                        {"field": "replicate_unit", "label": "Replicate unit", "type": "text"},
                        {"field": "comparison", "label": "Comparison", "type": "text"},
                        {"field": "principal_limit", "label": "Principal limitation", "type": "text"},
                    ],
                },
                {
                    "id": "validation_priorities",
                    "title": "Validation priorities",
                    "subtitle": "Experiments ordered by expected information gain and current evidence gaps",
                    "dataset": "validation_priorities",
                    "sourceId": "validation_priorities_sql",
                    "defaultSort": {"field": "priority", "direction": "asc"},
                    "columns": [
                        {"field": "priority", "label": "Priority", "type": "number"},
                        {"field": "target", "label": "Target", "type": "text"},
                        {"field": "why_now", "label": "Why now", "type": "text"},
                        {"field": "preferred_validation", "label": "Preferred validation", "type": "text"},
                        {"field": "decision_rule", "label": "Decision rule", "type": "text"},
                    ],
                },
            ],
            "sources": [{"id": item["id"], "label": item["label"], "path": item["path"]} for item in sources],
            "blocks": [
                {"id": "title", "type": "markdown", "body": "# PAH PDE8B–SMC evidence integration"},
                {
                    "id": "technical_summary",
                    "type": "markdown",
                    "body": (
                        "## Technical summary\n\n"
                        "The strongest result is a reproducible **PDE8B-centered lung association** that localizes to an SMC-remodeling context. "
                        "PDE8B elevation, loss of contractile programs, and gain of matrix programs replicate directionally across datasets, but patient-level PDE8B-to-state coupling is weak. "
                        "Ligand programs reproduce better than receptor programs; a complete SMC–pericyte communication architecture is not established. "
                        "PDE8B should remain a mechanism candidate—not a validated causal target."
                    ),
                },
                {
                    "id": "evidence_heading",
                    "type": "markdown",
                    "body": (
                        "## Replicated association is stronger than communication or causality\n\n"
                        "Grades A–B identify results that can anchor the manuscript. Grade C results should be framed as hypothesis-generating, and Grade D marks claims blocked by missing evidence. "
                        "The table preserves negative results and inference boundaries alongside positive findings."
                    ),
                },
                {"id": "grade_chart", "type": "chart", "chartId": "grade_distribution"},
                {"id": "evidence_table", "type": "table", "tableId": "evidence_matrix"},
                {
                    "id": "scope_heading",
                    "type": "markdown",
                    "body": (
                        "## Patient and cohort units define what can be inferred\n\n"
                        "Bulk cohorts establish reproducibility, while single-cell cohorts establish localization and state context. "
                        "GSE210248 cannot validate SMC–pericyte communication because it lacks a pericyte annotation; GSE293580 communication scores remain expression proxies with only two eligible donors."
                    ),
                },
                {"id": "cohort_table", "type": "table", "tableId": "cohort_context"},
                {
                    "id": "methodology",
                    "type": "markdown",
                    "body": (
                        "## Evidence grading favors replication and sensitivity over isolated P values\n\n"
                        "Grade A requires independent-cohort replication plus prespecified robustness checks. Grade B requires directional cross-dataset replication with acknowledged small-n limits. "
                        "Grade C denotes single-cohort, heterogeneous, or non-robust evidence. Grade D denotes an untested causal layer. "
                        "Exact permutation tests, Hedges g, threshold grids, annotation-boundary checks, and leave-one-patient/cohort analyses were retained. "
                        "Tables replace a summary chart because effect scales, cohort designs, and inference layers are not quantitatively interchangeable."
                    ),
                },
                {
                    "id": "limitations",
                    "type": "markdown",
                    "body": (
                        "## Robust directions coexist with decisive evidence gaps\n\n"
                        "The main limitations are two to three donors in single-cell comparisons, sparse expression, nonidentical SMC definitions, absent pericyte annotation in GSE210248, and lack of spatial, protein, perturbational, hemodynamic, or outcome validation. "
                        "Notably, within-PAH PDE8B–program correlations cross zero under leave-one-patient analysis, receptor-side communication directions are heterogeneous, and no ligand–receptor family survives multiplicity correction."
                    ),
                },
                {
                    "id": "recommendations",
                    "type": "markdown",
                    "body": (
                        "## Validation should first resolve localization, disease subtype, and pathway activity\n\n"
                        "The highest-information experiments directly test the remaining links in the argument: PDE8B-positive SMC state, SSc-PAH ECM–integrin localization, and TGF-beta pathway activity. "
                        "Perturbation should follow localization rather than precede it."
                    ),
                },
                {"id": "priority_table", "type": "table", "tableId": "validation_priorities"},
                {
                    "id": "further_questions",
                    "type": "markdown",
                    "body": (
                        "## Questions that could change the conclusion\n\n"
                        "- Does PDE8B colocalize with loss of mature SMC markers and gain of ECM markers in the same human cells?\n"
                        "- Is ECM–integrin exchange truly enriched in SSc-PAH after spatial and protein validation?\n"
                        "- Does PDE8B perturbation alter cAMP signaling and remodeling phenotypes in the same direction?\n"
                        "- Do larger, disease-stratified patient cohorts reproduce the receptor-side heterogeneity?"
                    ),
                },
            ],
        },
        "snapshot": {
            "version": 1,
            "generatedAt": generated,
            "status": "ready",
            "datasets": {
                "evidence_matrix": evidence,
                "evidence_grade_counts": grade_counts,
                "cohort_context": cohorts,
                "validation_priorities": priorities,
            },
            "accessIssues": [],
        },
        "sources": sources,
    }

    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    (REPORT_DIR / "grade_distribution.sql").write_text(grade_count_sql + "\n", encoding="utf-8")
    (REPORT_DIR / "evidence_matrix.sql").write_text(evidence_sql + "\n", encoding="utf-8")
    (REPORT_DIR / "cohort_context.sql").write_text(cohort_sql + "\n", encoding="utf-8")
    (REPORT_DIR / "validation_priorities.sql").write_text(priority_sql + "\n", encoding="utf-8")
    artifact_path = REPORT_DIR / "artifact.json"
    artifact_path.write_text(json.dumps(artifact, indent=2, ensure_ascii=False), encoding="utf-8")
    print(artifact_path)
    print(f"evidence_rows={len(evidence)} priority_rows={len(priorities)} cohort_rows={len(cohorts)}")


if __name__ == "__main__":
    main()
