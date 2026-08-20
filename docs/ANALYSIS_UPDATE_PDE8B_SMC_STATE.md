# PDE8B–SMC 状态转换阶段分析

## 项目信息

- 项目编号：`project-001`
- 项目名称：`project-001-pulmonary-arterial-hypertension-transcriptomics`
- 代码目录：`C:\workspace\project-001-pulmonary-arterial-hypertension-transcriptomics`
- 数据目录：`G:\workdata\projects\project-001-pulmonary-arterial-hypertension-transcriptomics`

## 已验证结果

1. 四个肺组织队列的 PDE8B 共表达核心包含 201 个基因。离一队列后，201/201 方向不变且均为 P<0.05，198/201 的最弱相关仍为 |r|≥0.30。
2. PDE8B 与 TRPC3、PRKG1、PDE5A、PLCB4、ROCK2、EDNRA、BMPR2 和 KCNK3 呈正相关。其中 BMPR2 的校正后 meta-r=0.411、FDR=6.96×10⁻⁶、I²=16.1%；KCNK3 的 meta-r=0.415，但异质性较高（I²=57.8%）。
3. GSE210248 中，正相关核心和收缩/cGMP 子模块在细胞类型间主要定位于 SMC1/SMC2。
4. 在 PAH 与 donor 的样本级比较中，PDE8B 在 SMC1/SMC2 升高（Hedges g=0.73/0.74），但经典收缩标志物模块下降（g=-0.83/-1.92）。
5. SMC2 中合成/ECM 模块升高（g=0.63），COL3A1、FN1、VIM 增强，而 ACTA2、TAGLN、MYH11、CNN1、MYL9 和 TPM2 下降。
6. PAH 中 SMC1 比例下降（g=-1.16），SMC2 比例上升（g=0.60）。SMC2 的完整模块评分与亚群比例呈负相关（Spearman rho=-0.83；n=6，仅作探索性证据）。
7. 独立 GSE293580 中，IPAH 的 PDE8B（g=1.33）、正核心模块（g=1.52）和收缩/cGMP 子模块（g=1.21）升高；SSc-PAH 中 PDE8B 仍升高（g=0.72），但完整模块效应较弱。

## 当前最合理模型

PDE8B 不是“所有收缩基因同步升高”的简单标志物。数据更支持：PDE8B 定位于肺血管平滑肌相关的环核苷酸/收缩分子背景，但在 PAH 中可伴随 SMC 从成熟收缩型向重塑或合成型状态转换而代偿性升高。该关系在 IPAH 与 SSc-PAH 之间可能不同。

因此，现阶段应使用“SMC 状态标志物”“表型重编程相关 PDE8B 上调”或“病因依赖性解耦”等表述，避免宣称 PDE8B 已被证明是共表达模块的因果上游调控因子。

## 初稿建议结果小标题

1. Cross-cohort meta-analysis identifies a robust PDE8B-centered PAH lung signature
2. The PDE8B co-expression core is enriched for cyclic-nucleotide and smooth-muscle programs
3. PDE8B localizes to vascular smooth muscle but is uncoupled from the mature contractile state in PAH
4. PDE8B upregulation accompanies a contractile-to-remodeling transition in PAH smooth muscle
5. PDE8B–module coupling differs between IPAH and SSc-PAH

## 最小后续验证

- 计算 PDE8B 与收缩/ECM 模块在各样本 SMC1/SMC2 内的相关性；由于 n=6，仅作为方向性证据。
- 若增加组织实验，优先 PDE8B+ACTA2/MYH11 和 PDE8B+COL3A1/FN1 双染，而不是扩大无针对性的分子实验。
- 若有人群样本，优先按 IPAH、SSc-PAH 分层，结合血流动力学和预后；不建议把所有 PAH 亚型直接合并。
- 初稿阶段保留轻量 CSV、SVG/PDF；最终图版确认后再统一生成投稿级高分辨率图片。
