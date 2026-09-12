# 变更记录

格式遵循 Keep a Changelog；版本号遵循语义化版本。

## [1.0.0] - 2026-09-12

本版本处理三位审稿人（rigor / compliance / writing）的 43 条意见（17 major、26 minor）：采纳 41 条、部分采纳 1 条、不采纳 1 条，逐条处理记录见 `docs/REVIEW_RESPONSE.md`。**四个问题的主结果未变**（t_end = 57.48 h 与 51.09 h，表 2–表 7 的全部数值不变）；变化的是检验的深度、表述的准确性与交付的完整性。

### Added
- **Duhamel 解析基准（MDR-0012）**：`verify.duhamel_field` 给出分段线性烘房温度下 Robin 圆柱的精确级数解，作为问题 1 温度场的**真解**对照（表格时刻误差 2.0e-6 °C、观测阶 1.97；result1.xlsx 温度表全部 37 800 个单元误差 2.0e-6 °C）。表 14 增列第三组。
- **完整输出网格的收敛研究（MDR-0011）**：`convergence._output_grid_study` 在 k=6、7 上重解问题 1、2 的 1 s 输出网格，按实测三网格阶作 Richardson 外推，逐格评估交付数组的误差并统计超过半个末位（5e-5）的单元；新增附录表 `tab_outputgrid`。
- **不受参照偏置的三网格观测阶（MDR-0011）**：`verify.convergence_table` 与 `verify.triplet_orders`；表 13、表 15 增列，表注写明自收敛估计量的偏置因子。
- **守恒型水分方程（MDR-0010）**：`ProblemSpec.mass_form`（drybasis | conservative）与 `physics.storage_coefficient`、`verify.conservative_moisture_balance`；作为第 (d) 处解释性不确定性进入 §1.3 与表 11（−18.71%）。
- **问题 2/3 的区间一致性校验**：`verification.horizon_consistency`（181 个共有 60 s 整点上差异为 0）。
- **补充交付 `result2_全过程.xlsx`**：同一模板、温度与水分浓度两张表、覆盖 0–t_end、60 s 分辨率（MDR-0005 修订；全过程每 1 s 的版本实测 57.3 MB，超出 20 MB 预算）。
- **能量上限的局部密度换算**：`FluxCap.local_rho_s`，潜热影响改以区间报告（+17.3% 至 +26.3%）。
- **数据自洽性与烘房统计**：`derived` 登记题给 R(t) 与 ρ(C) 在干物质守恒下的 13.4% 差、由 ρ(C) 反推的末半径 1.287 cm、末 1 h 窗口的噪声统计与 14400 s 处的 0.17 K 跳跃、附录 3 的初始传质 Biot 数 2.84、收缩幅度 s 的标度预测 −1.34。
- **结论节**：四问答案汇总表与三条跨问题结论；§5–§7、§9 补节首概述。
- MDR-0010、0011、0012；`docs/REVIEW_RESPONSE.md`。

### Changed
- **结果文件末行**：`solve_until_dry` 增加 `display_decimals`，result3/result4 写到四位小数下全场均 ≤0.1499 的首个 60 s 整点（行数 3449→3452、3066→3067），消除"显示为 0.1500"的误读。
- **论断改写**：§8.2 的"每一位小数都有意义"改为分层结论（温度全网格可靠；水分仅 42 与 29 个表面单元在最初数十秒内不可靠）；问题 3 的 t_end 改报为"57.48 h（生产网格）/ 57.47 h（网格收敛极限）"并给出 0.41 min 的界；Richardson 外推改用实测观测阶。
- 命题 4.2（离散极值原理）改用上 Dini 导数叙述，删去循环论证；假设 7 改为"径向均匀收缩（轴向长度不变）"。
- 附录 C 的复现命令由 `\nolinkurl` 改为 `\texttt`（`\nolinkurl` 会吞掉空格使命令不可复制）；`_reproduce_md` 输出真实阶段串与 tag。
- 图：`fig_properties` 纵轴限为可读范围、图例顺序统一、判据线加标注；`fig_sensitivity` 横轴收紧；技术路线图改为单一总线分四路。
- 表 9（时间积分无关性）与表 12（模型评价扩展）移入附录，为结论节腾出版面。

### Fixed
- §4.5 关于 r ≤ 1 cm 处水分仍为初值的表述与表 3 矛盾；§2.2 与 §5.2 跨附录混用 D 与 Biot 数；§7.2 关于 result4.xlsx 留空方式的表述；§9.2 "22 倍"的先行词；摘要中缺单位、能量不变量缺条件、"无插值"对问题 4 不成立；nfev 22510/22528 的不自洽；三张图与六个编号公式未被正文引用。

### Platform（`forge/`，最小改动）
- `forge/xlsx.py`：`SheetContract` 增加 `first_col_step`、`header_values`、`header_value_columns`，使 qa 真正检查时间列步长与首行刻度；`result4.xlsx` 补 `header_len=23`。
- `pipelines/common/stages.py`：复现记录表把"正在写该表的阶段"标注为 `completed*`；`REPRODUCE.md` 的复现命令由占位符改为真实阶段串与 tag。

## [0.9.0] - 2026-09-12

### Added
- 论文全文（`manuscript/sections/`）：摘要、问题重述（含三处不确定性的解释）、问题分析（TikZ 技术路线图）、编号假设与符号表、四个问题的模型建立/离散格式/算法/结果表/结果分析（含离散守恒与离散极值原理两个命题及证明）、模型检验（解析级数解、网格与时间收敛、Richardson 误差界、守恒与定性性质、独立二维轴对称交叉检验与端部效应、结果文件契约、灵敏度与不确定性预算、备选解释）、模型评价与推广、AI 使用详情与复现附录；参考文献 23 条（GB/T 7714）。
- `derived` 阶段与 MDR-0009：论文引用的无量纲数、时间尺度、轴向 Fourier 数、扩散系数量级比、烘干时间的 Richardson 外推误差界、温度平衡时刻、求解器统计与收缩几何比（39 个数值宏）。
- `manuscript/main.tex`：载入 gbt7714（缺失时回退 natbib+unsrtnat）、`\sci` 与 `\smalltable` 宏、"补充表格"附录（大表移出正文）。

### Fixed
- 参考文献包未载入导致 `\cite` 展开为完整作者串、产生 170 pt 溢出行；符号表（longtable）无题注却占用表号 1；`manuscript/ai_usage.tex` 独立编译时 `\ref{app:code}` 显示为 `??`（qa 的引用检查只覆盖主文档）。
- 结果文件末行的水分浓度按题目要求四舍五入到四位小数后显示为 `0.1500`（精确值略低于 0.15）；`derived` 阶段登记精确值，论文以“注”明示，避免误读。

### Verified（run 20260912-104036-184d092）
- `paper`：165 页，正文 30 页，摘要 1 页，0 处未定义引用，0 处溢出行，26 个图文件、15 张表。
- `qa`：16/16 通过（含 result1–4.xlsx 契约）；`lint` 0 项；`test` 49 项通过；`package` 257 个成员、漂移 0。
- 论文数字与交付工作簿的独立抽查 15/15 一致（q1/q2 的表格单元、q3/q4 的末行时刻与行数、`tab_q3_moist` 与 result3.xlsx 逐格一致）。

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
