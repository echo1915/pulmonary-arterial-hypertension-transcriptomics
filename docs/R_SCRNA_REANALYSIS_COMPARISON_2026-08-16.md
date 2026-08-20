# R/pahscrna 单细胞重分析比较（2026-08-16）

## 结论

使用当前 `pahscrna` R 包对既往单细胞结果进行独立复算后，主要结论没有改变。GSE210248 的 69 个候选基因、7 个保留细胞类型共 483 个基因—细胞类型比较中，效应方向一致率为 100%；新旧 `delta_log_cp10k` 的最大绝对差为 `4.44e-16`，Hedges' g 的 Spearman 相关为 0.9994。

PDE8B 在 GSE210248 的 SMC 1 和 SMC 2 中仍呈 PAH 高于供体的描述性效应，Hedges' g 分别为 0.73 和 0.74。但 PDE8B 原始计数较稀疏，没有通过 edgeR `filterByExpr`，因此不能把该方向性结果表述为经过原始计数差异表达模型验证的显著结果。

GSE293580 的严格 SMC 患者级验证仍为低样本量证据：IPAH 对供体的 Hedges' g 为 11.06、精确置换 p = 0.10；SSc-PAH 对供体的 Hedges' g 为 1.37、精确置换 p = 0.267。供体仅 2 例，因此只作为定位/方向性支持，不作为独立确证。

## 方法和复现边界

- GSE210248 从原始 10x 稀疏计数矩阵读入，在 R 中按“受试者 × 细胞类型”聚合；最低 20 个细胞，且每组至少 3 名受试者。
- 所有统计比较以受试者为生物学重复；细胞没有被当作独立样本。
- 描述性结果使用 `log(1 + CP10K)`、Hedges' g 和精确置换检验；候选集内 p 值采用 BH 校正。
- edgeR 使用受试者级伪 bulk 原始计数、TMM 标准化和 quasi-likelihood 模型；PDE8B 因低表达过滤未进入 SMC 模型。
- GSE293580 的作者 RDS 当前不可读，H5AD 的 nullable-string-array 也不受现有 R 读取器支持。因此该队列从既有、已锁定的患者级伪 bulk CSV 重新执行统计检验，而不是在 R 中从细胞矩阵重新聚合。这是本次复算的主要技术限制。

## 图形优化

新版图采用 183 × 128 mm 双栏布局，统一 Arial、受试者散点和均值横线，同时直接标明 483 个比较的方向一致性、PDE8B 的描述性 SMC 效应及 edgeR 过滤限制。提供 SVG、PDF、600 dpi TIFF 和 300 dpi PNG。

图中：

- a：旧分析与 R/pahscrna 复算的效应一致性；
- b：PDE8B 在 SMC 1/2 的受试者级描述性 Hedges' g；
- c：GSE210248 中每名受试者的 PDE8B 表达；
- d：GSE293580 严格 SMC 患者级低样本量验证。

## 输出

- 脚本：`scripts/pah/reanalyse_scrna_with_pahscrna.R`
- 汇总表：`outputs/r-reanalysis-2026-08-16/tables/reanalysis_summary.csv`
- 新旧比较：`outputs/r-reanalysis-2026-08-16/tables/gse210248_old_vs_R_comparison.csv`
- GSE210248 候选结果：`outputs/r-reanalysis-2026-08-16/tables/gse210248_candidate_effects_R.csv`
- GSE293580 检验：`outputs/r-reanalysis-2026-08-16/tables/gse293580_strict_SMC_tests_R.csv`
- 图件：`outputs/r-reanalysis-2026-08-16/figures/Figure_R_reanalysis_PDE8B_SMC_comparison.*`
- R 环境记录：`outputs/r-reanalysis-2026-08-16/sessionInfo.txt`

上述 `outputs/` 均位于外置数据目录，不在代码仓库内。

## 项目标识

- 项目编号：001
- 标准项目名：`project-001-pulmonary-arterial-hypertension-transcriptomics`
- 代码目录：`C:\workspace\project-001-pulmonary-arterial-hypertension-transcriptomics`
- 外置数据目录：`G:\workdata\projects\project-001-pulmonary-arterial-hypertension-transcriptomics`
