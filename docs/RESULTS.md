# 结果说明（供论文撰写使用）

本文件汇总最终云端运行的全部科学结果：每个问题的方法、关键数值、表格/图文件、`ctx.number` 键名、验证结论、局限与备选解释。所有数值均来自 `numbers.json`（论文中用 `\val<Key>` 引用）与 `tables` 阶段生成的 `tables/*.tex`，禁止手抄。

- 最终 run：`20260912-153103-eb858bf`（代码 `eb858bf` 及其后续文档/论文提交；v1.0.0 发布）。上一轮为 `20260912-104036-184d092`（v0.9.0）。
- 阶段：ingest, validate, q1, q2, q3, q4, convergence, verification, sensitivity, extension, derived, results, figures, tables, lint, test, paper, qa, package, release（全部 completed；qa 含 `result:*` 契约检查）。
- **本轮（v1.0.0，审稿意见处理）的变化**：主结果未变（t_end = 57.4805 h 与 51.0893 h，表 1–表 6 的全部数值不变）。新增/改变的是检验与交付：
  - `convergence`：新增结果文件完整输出网格（1 s × 0.1 cm）上的误差研究（k=6,7 两级加密 + Richardson 极限），并为所有收敛表增加不受参照偏置的**三网格观测阶**（MDR-0011）。
  - `verification`：新增问题 1 温度场的 **Duhamel 真解**对照（MDR-0012）与问题 2/3 的**区间一致性校验**（181 个共有 60 s 整点上逐位相同）。
  - `sensitivity`：新增守恒型水分方程备选 `alt:mass_conservative`（MDR-0010，−18.71%）。
  - `extension`：能量上限增加局部 ρ_s(C_s) 换算的变体（+26.3% vs 常数换算的 +17.3%）。
  - `derived`：题给 R(t) 与 ρ(C) 的干物质自洽性（差 13.4%）、烘房平台噪声统计、附录 3 的初始传质 Biot 数、收缩幅度 s 的标度预测、按实测观测阶的 Richardson 外推。
  - `results`：result3/result4 的末行延长到四位小数下全场 ≤0.1499 的首个 60 s 整点（行数 3452 / 3067）；新增补充交付 `result2_全过程.xlsx`。
  - 审稿意见的逐条处理见 `docs/REVIEW_RESPONSE.md`。
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

- 方法一句话：附录 3 变物性（ρ、c_p、k 随 C；D 随 C、T），附件 1 逐点插值的烘房条件，0–3 h 每 1 s 输出（MDR-0005：result2.xlsx 覆盖 3 h；全过程的温度与水分另由 `result2_全过程.xlsx` 以 60 s 分辨率交付，全过程每 1 s 的版本实测 57.3 MB，超出 20 MB 预算）。
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
- 时间积分无关性：表 `tab_time_independence.tex`——更严格容差与 Radau（`\valTimeStrictDiffMoist`、`\valTimeRadauDiffMoist`；烘干时间差 `\valTimeStrictDryDiffSec`、`\valTimeRadauDryDiffSec` s）以及限制最大步长 Δt ≤ 3600/900/300 s 的自适应 BDF（`\valTimeStepCoarseDiffMoist`、`\valTimeStepCoarseDryDiffSec` s；`\valTimeStepFineDiffMoist`、`\valTimeStepFineDryDiffSec` s）：差异均在积分容差量级，结果与步长序列无关。
- 守恒：图 `fig_conservation.pdf`；表 `tab_verification.tex`。
- 独立二维轴对称程序：端面绝热时与一维解之差 `\valCrossCheckDiffMoist` / `\valCrossCheckDiffTemp`（代码独立性）；端面 Robin 时中截面差 `\valEndEffectMidplaneDiff`，烘干时间差 `\valEndEffectDryDiffPct`%（`\valEndEffectDryHoursTwoD` vs `\valEndEffectDryHoursOneD` h）；图 `fig_endeffect.pdf`。

## 灵敏度与备选解释（sensitivity）

- 表 `tab_sensitivity.tex`：h、h_m、D₀、C_air、T 平台、收缩幅度 s 的 ×0.8/0.9/1.1/1.2 烘干时间与弹性（键 ElasH, ElasHm, ElasDzero, ElasCair, ElasTplateau, ElasShrink, ElasDzeroQfour）；图 `fig_sensitivity.pdf`。
- 表 `tab_alternatives.tex`：平台取末样本/名义值、界面系数算术/调和平均、欧拉写法、固定半径、平均含水率判据（键 Alt…Hours / Alt…Pct）。

## 模型评价扩展（extension，MDR-0008）——供"模型评价与推广"一节使用

主结果严格按题给参数与边界条件计算；本阶段量化题给模型省略的两项物理机制（问题 3 设定，固定半径、附录 3 物性），不改变问题 1–4 的答案。

- 潜热一致性诊断（后验，由问题 3 的解计算）：题给传质律隐含的蒸发质量通量 $j_w=\rho_{s0}h_m(C_s-C_{\mathrm{air}})$（干固体密度 $\rho_{s0}=\rho(C_0)/(1+C_0)$=`\valExtDrySolidDensity` kg/m³），其潜热负荷折算为"等效表面温降" $\Delta T_{\mathrm{lat}}=L_v j_w/h$：初始 `\valExtLatentDeltaTPeak` K，3 h 时 `\valExtLatentDeltaTThreeH` K，12 h 时 `\valExtLatentDeltaTTwelveH` K，24 h 时 `\valExtLatentDeltaTTwentyFourH` K；降到 1 K 以下的时刻 `\valExtLatentBelowOneKHours` h（0.1 K 以下 `\valExtLatentBelowTenthKHours` h）。全程蒸发水量 `\valExtWaterRemovedKg` kg/m（单位长度），潜热 `\valExtLatentEnergyMJ` MJ/m，而题给能量边界条件全程供入的对流热量仅 `\valExtConvEnergyKJ` kJ/m（= 显热储量），比值 `\valExtLatentToConvRatio`。结论：题给模型在前十余小时高估药材温度（因而高估 D），后期潜热可忽略；这是题给参数体系的固有简化（MDR-0008 备选 1 说明直接耦合不适定）。
- 烘房状态的湿空气解读（把"水分浓度"读作含湿量 kg/kg 干空气）：初始 28 °C、RH `\valExtRhInitPct`%（湿球 `\valExtWetBulbInit` °C）；恒温阶段 50 °C、RH `\valExtRhPlateauPct`%，湿球 `\valExtWetBulbPlateau` °C，湿球温差 `\valExtWetBulbDepression` K；对流供热可维持的蒸发通量上限 $j_{\max}=h(T_{\mathrm{air}}-T_{\mathrm{wb}})/L_v$ = `\valExtFluxCapPlateau` kg/(m²·s) = `\valExtFluxCapPlateauKgPerHour` kg/(m²·h)，而题给模型的初始蒸发通量 `\valExtModelFluxInit` kg/(m²·s) 为初始上限的 `\valExtFluxRatioInit` 倍。
- 情形 A（能量限制/恒速段）：蒸发通量取 $\min\{h_m(C_s-C_{\mathrm{air}}),\,j_{\max}(t)/\rho_{s0}\}$，烘干时间 `\valExtCappedHours` h（相对主结果 `\valExtCappedPct`%），上限在 `\valExtCappedReleaseHours` h 后不再起作用（此后与题给模型相同的降速段）；上限 ×0.5、×2 的结果 `\valExtCapOneHours` / `\valExtCapThreeHours` h（`\valExtCapOnePct`% / `\valExtCapThreePct`%）。纯能量下界（全程按湿球极限蒸发）`\valExtEnergyBoundHours` h（`\valExtEnergyBoundPct`%）——主结果远高于该下界，说明题给模型的烘干时长由内部扩散控制而非能量供给控制。
- 情形 B（平衡含水率驱动力 $C_s-C_{\mathrm{eq}}$）：$C_{\mathrm{eq}}$=`\valExtCeqOneValue`/`\valExtCeqTwoValue`/`\valExtCeqThreeValue` 时烘干时间 `\valExtCeqOneHours`/`\valExtCeqTwoHours`/`\valExtCeqThreeHours` h（`\valExtCeqOnePct`/`\valExtCeqTwoPct`/`\valExtCeqThreePct`%）；$C_{\mathrm{eq}}\to0.15$ 时烘干时间趋于无穷，说明 0.15 kg/kg 的终点判据必须与烘房湿度（等温吸附曲线）匹配。注意该单调性不是普遍的：驱动力越小表面越湿，而 $D\propto\exp(-a/C)$ 使湿表面的扩散系数大得多——在附录 2 的律（$a=0.89$）下 $C_{\mathrm{eq}}=0.10$ 反而比 0.05 略快（干表面形成低 $D$ 硬壳，"case hardening"；见 `tests/unit/test_extension.py` 的说明），附录 3（$a=0.45$）下则单调变慢；常数 $D$ 时可证单调（线性问题，阈值 $(0.15-C_{\mathrm{eq}})/(2.55-C_{\mathrm{eq}})$ 随 $C_{\mathrm{eq}}$ 减小）。
- 表：`tab_extension.tex`（题给模型、能量限制 ×0.5/1/2、C_eq 三档、纯能量下界）。图：`fig_extension.pdf`（建议图注："模型评价扩展：(a) 题给模型蒸发率的潜热等效温降 $L_v j_w/h$ 与模型实际温差 $T_{\mathrm{air}}-T_s$（对数坐标，点线为 1 K）；(b) 题给模型、能量限制与平衡含水率驱动力三种情形下的中心水分时程，点线为 0.15 kg/kg。"；说明：(a) 中模型温差在 2–4 h 间因烘房温度样本噪声（±0.3 °C）而在 0 附近改变符号，负值在对数坐标下被截断，故呈振荡带）。
- numbers 键：ExtDrySolidDensity, ExtRhInitPct, ExtRhPlateauPct, ExtWetBulbInit, ExtWetBulbPlateau, ExtWetBulbDepression, ExtLatentHeatMJ, ExtFluxCapPlateau, ExtFluxCapPlateauKgPerHour, ExtFluxCapInit, ExtModelFluxInit, ExtFluxRatioInit, ExtLatentDeltaTPeak, ExtLatentDeltaTThreeH, ExtLatentDeltaTTwelveH, ExtLatentDeltaTTwentyFourH, ExtLatentBelowOneKHours, ExtLatentBelowTenthKHours, ExtLatentEnergyMJ, ExtConvEnergyKJ, ExtLatentToConvRatio, ExtWaterRemovedKg, ExtEnergyBoundHours, ExtEnergyBoundPct, ExtCappedHours, ExtCappedPct, ExtCappedReleaseHours, ExtCap{One,Two,Three}{Factor,Hours,Pct,ReleaseHours}, ExtCeq{One,Two,Three}{Value,Hours,Pct}。
- 验证：`extension/verification_extension.json`——每个情形的守恒（含通量上限的损失律）、极值原理、单调性与判据一致性检查全部通过。
- 假设与局限：湿空气解读依赖"水分浓度=含湿量"的读法（DATA_NOTES）；$L_v$、Magnus 公式、ASHRAE 湿球方程为标准常数/公式；$\rho_{s0}$ 取初始值；情形 A 只限制质量通量、不改变能量方程（温度仍按题给模型），因此仍是界定性分析而非完整耦合模型。

## 派生量（derived，MDR-0009）

`derived` 阶段只读取 q1–q4 与 convergence 的已存产物，不做新的模拟，登记论文分析与检验部分引用的派生量：
长径比与未知量个数（`\valAspectRatio`、`\valSystemSize`）；问题 1 的热/质扩散时间尺度、Fourier 数、扩散长度与刚性比
（`\valQoneHeatTimeScaleS`、`\valQoneMassTimeScaleH`、`\valQoneFourierHeat`、`\valQoneFourierMass`、
`\valQoneDiffusionLengthCm`、`\valQoneStiffnessRatio`）；问题 2 中心温度与烘房相差 1 K / 0.5 K 的时刻
（`\valQtwoCentreWithinOneKelvinHours`、`\valQtwoCentreWithinHalfKelvinHours`）；一维化论证的轴向 Fourier 数
（`\valFoAxialHeat`、`\valFoAxialMassUpper`）；附录 3/4 扩散系数的量级与比值（`\valDmaxAppendixThree`、
`\valDtargetAppendixThree`、`\valDsurfaceAtDryQthree`、`\valDratioQthree`、`\valBiotMassAtDryQthree`、
`\valDmaxAppendixFour`、`\valDtargetAppendixFour`、`\valDratioThreeToFourInit`）；烘干时间的 Richardson 外推与
生产网格误差（`\valQthreeDryHoursExtrapolated`、`\valQthreeGridErrorMin`、`\valQfourDryHoursExtrapolated`、
`\valQfourGridErrorMin`）；求解器统计（`\valQthreeNfev`、`\valQthreeNlu`、`\valQthreeSegments`、`\valQfourNfev`）；
收缩几何比与固定/移动半径的时间比（`\valQfourAreaRatioPct`、`\valQfourRadiusRatioSquared`、
`\valQfourFixedToMovingRatio`）；数据规模（`\valChamberSamples`、`\valChamberEndHours`、`\valRadiusSamples`、
`\valRadiusEndHours`）。

## 论文（paper、qa）

最终排版同属 run `20260912-104036-184d092`：165 页（摘要 1 页 + 正文 30 页 + 附录），0 处未定义引用、0 处溢出行、
字体全嵌入；qa 16/16 通过（含四个结果工作簿的模板契约）。大表（网格收敛、解析解对照、烘干时间收敛、问题 4 备选写法）
置于附录 A"补充表格"，正文保留题目要求格式的表 2–表 7。

## 局限

平衡含水率与等温吸附未建模（表面条件为题给线性形式；`extension` 给出 C_eq 的界定分析）；蒸发潜热对能量方程的耦合未计入（`extension` 给出后验诊断与能量限制情形）；各向同性均匀收缩假设；h、h_m 恒定；一维近似（端部效应已量化）；附件 1 之后的烘房条件为平台外推。

## 数据与图表文件清单

见 `figures/figure_index.json`（每张图用途）与 `tables/*.csv`。
