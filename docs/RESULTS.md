# 结果说明（供论文撰写使用）

本文件汇总最终云端运行的全部科学结果：每个问题的方法、关键数值、表格/图文件、`ctx.number` 键名、验证结论、局限与备选解释。所有数值均来自 `numbers.json`（论文中用 `\val<Key>` 引用）与 `tables` 阶段生成的 `tables/*.tex`，禁止手抄。

- 最终 run：`20260912-042323-b639c9c`（代码 `b639c9c`）。
- 阶段：ingest, validate, q1, q2, q3, q4, convergence, verification, sensitivity, results, figures, tables, lint, test（全部 completed）。
- 生产网格：Δr = 0.1/2⁵ = 0.003125 cm（641 节点，`\valGridNodes`、`\valGridDrCm`）；时间积分 BDF，rtol=1e-9、atol=1e-11（`\valTimeRtol`、`\valTimeAtol`），在附件 1 的每个样本时刻与附件 2 的每个节点重启（MDR-0004）。

## 统一模型（MDR-0001 … 0004）

一维径向热质耦合方程（圆柱、无限长近似，中截面）：ρc_p ∂T/∂t = (1/r)∂_r(k r ∂_r T)，∂C/∂t = (1/r)∂_r(D r ∂_r C)；中心对称；表面 −k∂_r T = h(T−T_air(t))，−D∂_r C = h_m(C−C_air(t))；初值 28 °C、2.55 kg/kg。节点型有限体积（r 加权控制体，界面系数取界面状态法），自适应 BDF；输出时刻由稠密插值精确采样。h=25、h_m=8×10⁻⁷ 全程沿用（MDR-0002）；14400 s 后烘房条件取末 1 h 均值平台 T_air=`\valQthreeAirTempPlateau` °C、C_air=`\valQthreeAirHumPlateau`（MDR-0003）。

## 问题 1（q1）

- 方法一句话：附录 2 常物性 + D(C)，0–1800 s，每 1 s 输出，节点与 0.1 cm 网格重合。
- 关键数值：1800 s 时中心温度 `\valQoneCentreTempEnd` °C、表面温度 `\valQoneSurfaceTempEnd` °C（烘房 `\valQoneAirTempEnd` °C，表面滞后 `\valQoneSurfaceTempLag` K）；中心水分 `\valQoneCentreMoistEnd`、表面水分 `\valQoneSurfaceMoistEnd`、体积平均 `\valQoneMeanMoistEnd`，30 min 失水 `\valQoneMoistLossPct`%。Bi=`\valQoneBiotHeat`、Bi_m=`\valQoneBiotMass`、α=`\valQoneAlpha` m²/s、D(C₀)=`\valQoneDinit` m²/s。
- 表：`tables/tab_q1_temp.tex`（表 1）、`tables/tab_q1_moist.tex`（表 2）；CSV 同名。
- 图：`figures/fig_q1_profiles.pdf`（建议图注："问题 1 在表 1/2 各时刻的温度（a）与水分浓度（b）径向剖面；热扩散在 30 min 内已渗透到中心，而水分边界层仍限于表面约 0.5 cm。"）；`figures/fig_q1_history.pdf`（"问题 1 中心、表面与烘房温度（a）及中心、表面、平均水分浓度（b）时程。"）。
- numbers 键：QoneCentreTempEnd, QoneSurfaceTempEnd, QoneCentreMoistEnd, QoneSurfaceMoistEnd, QoneMeanMoistEnd, QoneMoistLossPct, QoneBiotHeat, QoneBiotMass, QoneAlpha, QoneDinit, QoneAirTempEnd, QoneSurfaceTempLag, GridNodes, GridDrCm, TimeRtol, TimeAtol。
- 结果文件：`results/result1.xlsx`（温度、水分浓度两表，t=1…1800 s × r=0…2 cm）。
- 验证：`q1/verification_q1.json`——水分与能量的离散守恒（积分器内累积边界通量，相对残差 ~1e-15）、极值原理、径向单调性，全部通过；网格收敛（表 `tab_convergence`）、贝塞尔级数解对照（表 `tab_analytic`）。

## 问题 2（q2）

- 方法一句话：附录 3 变物性（ρ、c_p、k 随 C；D 随 C、T），附件 1 逐点插值的烘房条件，0–3 h 每 1 s 输出（MDR-0005：result2 覆盖 3 h）。
- 关键数值：3 h 时中心/表面温度 `\valQtwoCentreTempThreeH` / `\valQtwoSurfaceTempThreeH` °C；中心/表面/平均水分 `\valQtwoCentreMoistThreeH` / `\valQtwoSurfaceMoistThreeH` / `\valQtwoMeanMoistThreeH`，3 h 失水 `\valQtwoMoistLossPct`%；初始 D=`\valQtwoDinit`，3 h 表面 D=`\valQtwoDsurfaceThreeH` m²/s，Bi_m(3 h)=`\valQtwoBiotMassThreeH`。
- 表：`tab_q2_temp.tex`（表 3）、`tab_q2_moist.tex`（表 4）。
- 图：`fig_q2_profiles.pdf`（"问题 2 各半小时时刻的温度（a）与水分（b）剖面"）、`fig_q2_history.pdf`（"问题 2 温度与水分时程；温度在约 2 h 内与烘房平衡"）。
- numbers 键：QtwoCentreTempThreeH, QtwoSurfaceTempThreeH, QtwoCentreMoistThreeH, QtwoSurfaceMoistThreeH, QtwoMeanMoistThreeH, QtwoMoistLossPct, QtwoDinit, QtwoDsurfaceThreeH, QtwoBiotMassThreeH, QtwoAlphaInit。
- 结果文件：`results/result2.xlsx`（t=1…10800 s）。
- 验证：`q2/verification_q2.json`（守恒、极值、单调性通过）；收敛表、贝塞尔对照同上。

## 问题 3（q3）

- 方法一句话：问题 2 模型延续至 max_r C < 0.15 kg/kg（事件求根精确判定，MDR-0005），每 60 s 输出。
- 关键数值：烘干时间 `\valQthreeDryHours` h（`\valQthreeDryHoursFour` h，`\valQthreeDrySeconds` s，`\valQthreeDryDays` 天）；60 s 分辨率 `\valQthreeDryHoursGrid` h；结束时表面/平均水分 `\valQthreeSurfaceMoistAtDry` / `\valQthreeMeanMoistAtDry`；平均含水率判据 `\valQthreeMeanCriterionHours` h；失水一半 `\valQthreeHalfLossHours` h；result3 行数 `\valQthreeResultRows`。
- 表：`tab_q3_moist.tex`（表 5，每 6 h + 烘干结束时间行）。
- 图：`fig_q3.pdf`（"问题 3：每 6 h 与结束时刻的水分剖面（a）及中心、表面、平均水分时程（b）；虚线为 0.15 kg/kg"）。
- numbers 键：QthreeDryHours, QthreeDryHoursFour, QthreeDryHoursGrid, QthreeDrySeconds, QthreeDryDays, QthreeSurfaceMoistAtDry, QthreeMeanMoistAtDry, QthreeMeanCriterionHours, QthreeHalfLossHours, QthreeResultRows, QthreeAirTempPlateau, QthreeAirHumPlateau。
- 结果文件：`results/result3.xlsx`（t=60…结束的 60 s 整点）。
- 验证：`q3/verification_q3.json`；烘干时间网格收敛（`tab_convergence_drytime`：与最细网格之差 `\valConvDryDiffQthreeMin` min）；时间积分无关性（`tab_time_independence`：`\valTimeStrictDryDiffSec` s、`\valTimeRadauDryDiffSec` s）；二维端部效应（`\valEndEffectDryDiffPct`%）。

## 问题 4（q4）

- 方法一句话：附录 4 物性 + 附件 2 半径 R(t)（PCHIP），干物质守恒的拉格朗日 Landau 坐标（MDR-0006）；备选欧拉写法与固定半径对照。
- 关键数值：烘干时间 `\valQfourDryHours` h（`\valQfourDryHoursFour`，60 s 分辨率 `\valQfourDryHoursGrid`）；欧拉备选 `\valQfourDryHoursAlt` h（相对 `\valQfourAltDiffPct`%）；半径固定 2 cm 时 `\valQfourDryHoursFixed` h（`\valQfourFixedDiffPct`%）；结束时半径 `\valQfourRadiusAtDryCm` cm；退化检验最大差 `\valQfourDegenerateMaxDiff`。
- 表：`tab_q4_moist.tex`（表 6，含"药材表面"列；超出当前半径处为 --）、`tab_q4_moist_alt.tex`（备选写法）。
- 图：`fig_q4.pdf`（"问题 4：物理坐标下每 6 h 的水分剖面与表面位置（a）；中心水分时程在主模型、欧拉备选与固定半径下的比较（b）"）。
- numbers 键：QfourDryHours, QfourDryHoursFour, QfourDryHoursGrid, QfourDrySeconds, QfourDryDays, QfourDryHoursAlt, QfourDryHoursFixed, QfourAltDiffPct, QfourFixedDiffPct, QfourRadiusAtDryCm, QfourSurfaceMoistAtDry, QfourMeanMoistAtDry, QfourResultRows, QfourDegenerateMaxDiff。
- 结果文件：`results/result4.xlsx`（每 60 s；列 0…2 cm，超出当前半径留空，末列"药材表面"）。
- 验证：`q4/verification_q4.json`（拉格朗日写法水分守恒相对残差 ~1e-15；欧拉写法通过退缩边界伪损失约 34% 的水分——见 `alternative_formulation.moisture_balance`）；收敛表；退化检验。

## 模型检验汇总（verification, convergence）

- 贝塞尔级数解对照（常物性、恒定环境）：温度 `\valAnalyticErrTemp`、水分 `\valAnalyticErrMoist`（生产网格）；观测阶 `\valAnalyticOrderTemp` / `\valAnalyticOrderMoist`；表 `tab_analytic.tex`，图 `fig_convergence.pdf`(b)。
- 网格收敛：表 `tab_convergence.tex`（观测阶 `\valConvOrderQoneTemp` 等）、`tab_convergence_drytime.tex`；图 `fig_convergence.pdf`(a)。
- 时间积分无关性：表 `tab_time_independence.tex`（`\valTimeStrictDiffMoist`、`\valTimeRadauDiffMoist`）。
- 守恒：图 `fig_conservation.pdf`；表 `tab_verification.tex`。
- 独立二维轴对称程序：端面绝热时与一维解之差 `\valCrossCheckDiffMoist` / `\valCrossCheckDiffTemp`（代码独立性）；端面 Robin 时中截面差 `\valEndEffectMidplaneDiff`，烘干时间差 `\valEndEffectDryDiffPct`%（`\valEndEffectDryHoursTwoD` vs `\valEndEffectDryHoursOneD` h）；图 `fig_endeffect.pdf`。

## 灵敏度与备选解释（sensitivity）

- 表 `tab_sensitivity.tex`：h、h_m、D₀、C_air、T 平台、收缩幅度 s 的 ×0.8/0.9/1.1/1.2 烘干时间与弹性（键 ElasH, ElasHm, ElasDzero, ElasCair, ElasTplateau, ElasShrink, ElasDzeroQfour）；图 `fig_sensitivity.pdf`。
- 表 `tab_alternatives.tex`：平台取末样本/名义值、界面系数算术/调和平均、欧拉写法、固定半径、平均含水率判据（键 Alt…Hours / Alt…Pct）。

## 局限

平衡含水率与等温吸附未建模（表面条件为题给线性形式）；蒸发潜热对能量方程的耦合未计入；各向同性均匀收缩假设；h、h_m 恒定；一维近似（端部效应已量化）；附件 1 之后的烘房条件为平台外推。

## 数据与图表文件清单

见 `figures/figure_index.json`（每张图用途）与 `tables/*.csv`。
