# 变更记录

格式遵循 Keep a Changelog；版本号遵循语义化版本。

## [0.9.0] - 2026-09-12

### Added
- 论文全文（`manuscript/sections/`）：摘要、问题重述（含三处不确定性的解释）、问题分析（TikZ 技术路线图）、编号假设与符号表、四个问题的模型建立/离散格式/算法/结果表/结果分析（含离散守恒与离散极值原理两个命题及证明）、模型检验（解析级数解、网格与时间收敛、Richardson 误差界、守恒与定性性质、独立二维轴对称交叉检验与端部效应、结果文件契约、灵敏度与不确定性预算、备选解释）、模型评价与推广、AI 使用详情与复现附录；参考文献 23 条（GB/T 7714）。
- `derived` 阶段与 MDR-0009：论文引用的无量纲数、时间尺度、轴向 Fourier 数、扩散系数量级比、烘干时间的 Richardson 外推误差界、温度平衡时刻、求解器统计与收缩几何比（39 个数值宏）。
- `manuscript/main.tex`：载入 gbt7714（缺失时回退 natbib+unsrtnat）、`\sci` 与 `\smalltable` 宏、"补充表格"附录（大表移出正文）。

### Fixed
- 参考文献包未载入导致 `\cite` 展开为完整作者串、产生 170 pt 溢出行；符号表（longtable）无题注却占用表号 1。

### Verified（run 20260912-104036-184d092）
- `paper`：165 页，正文 30 页，摘要 1 页，0 处未定义引用，0 处溢出行，26 个图文件、15 张表。
- `qa`：16/16 通过（含 result1–4.xlsx 契约）；`lint` 0 项；`test` 48 项通过。

## [Unreleased]

### Added
- 问题 A 科学阶段：`q1`–`q4`（径向有限体积 + 自适应 BDF，Landau 坐标处理收缩边界）、`convergence`（网格/时间收敛）、`verification`（贝塞尔级数解对照、独立二维轴对称程序交叉检验与端部效应）、`sensitivity`（弹性与备选解释）；`results`（result1–4.xlsx）、`figures`、`tables`。
- 模块：`pipelines/a/physics.py`（附录物性、烘房条件、半径）、`solver.py`、`verify.py`（与求解器分离的校验器）、`common.py`、`science.py`、`outputs.py`、`figures.py`。
- MDR-0001 … MDR-0007；`docs/DATA_NOTES.md`；`docs/RESULTS.md`；`tests/unit/test_drying.py`（物理、求解器、校验器的单元/性质测试）。
- `configs/default.toml` 新增 `[a.*]` 配置段；`paper.sources` 与 `package.stages` 登记新增阶段。
- `extension` 阶段（MDR-0008，模型评价）：问题 3 解的蒸发潜热一致性后验诊断、烘房湿空气解读（Magnus 饱和蒸气压、ASHRAE 湿球温度）、湿球能量上限下的恒速段变体与纯能量下界、平衡含水率驱动力 C_eq 的界定分析；`tab_extension`、`fig_extension`；求解器增加可选的蒸发通量上限（默认关闭，守恒不变量保持）；`tests/unit/test_extension.py`。
- `docs/mdr/0008-*.md`；`[a.extension]` 配置段；`package.stages` 登记 `extension`。
- `convergence` 时间无关性研究增加限制最大步长（Δt ≤ 3600/900/300 s）的行（键 TimeStep{Coarse,Fine}{DiffMoist,DryDiffSec}），`tab_time_independence` 相应扩展。

### Verified（run 20260912-104036-184d092，代码 184d092）
- 17 个阶段（ingest … test, paper, qa）全部 completed；qa 16/16 通过（含 result1–4.xlsx 契约）；lint 0 项；45 项测试通过。
- 独立校验：q1–q4、verification（贝塞尔级数解 1×10⁻⁵、二维交叉检验）、extension 的校验报告全部 `passed=true`。
- 身份信息观察名单命中 5 处，均为附录程序清单中 `configs/default.toml` 打印的观察名单模式本身（"大学、学院、赛区、队号、指导教师"），不含任何身份信息；复核结论：无需处理。

## [0.1.0] - 2026-09-12

### Added
- Forge 云端运行时（阶段/运行/清单模型、数据与结果契约、事件日志、确定性打包）。
- 通用流水线：ingest、validate、fmt、lint、test、paper、qa、package、release；全部在 Modal 执行。
- CUMCM 格式论文骨架（XeLaTeX + BibTeX gbt7714），自动生成数字宏、支撑材料清单、复现记录与完整程序附录。
- 工程规范、建模与写作标准、运行手册、ADR/MDR、AI 使用记录、Unlicense。

### Verified
- 三题仓库全链路冒烟通过：25 项单元/性质/集成测试，13 项论文与结果 QA 检查，发布散列核对一致。
