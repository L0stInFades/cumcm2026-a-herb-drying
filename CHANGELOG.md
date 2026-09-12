# 变更记录

格式遵循 Keep a Changelog；版本号遵循语义化版本。

## [Unreleased]

### Added
- 问题 A 科学阶段：`q1`–`q4`（径向有限体积 + 自适应 BDF，Landau 坐标处理收缩边界）、`convergence`（网格/时间收敛）、`verification`（贝塞尔级数解对照、独立二维轴对称程序交叉检验与端部效应）、`sensitivity`（弹性与备选解释）；`results`（result1–4.xlsx）、`figures`、`tables`。
- 模块：`pipelines/a/physics.py`（附录物性、烘房条件、半径）、`solver.py`、`verify.py`（与求解器分离的校验器）、`common.py`、`science.py`、`outputs.py`、`figures.py`。
- MDR-0001 … MDR-0007；`docs/DATA_NOTES.md`；`docs/RESULTS.md`；`tests/unit/test_drying.py`（物理、求解器、校验器的单元/性质测试）。
- `configs/default.toml` 新增 `[a.*]` 配置段；`paper.sources` 与 `package.stages` 登记新增阶段。
- `extension` 阶段（MDR-0008，模型评价）：问题 3 解的蒸发潜热一致性后验诊断、烘房湿空气解读（Magnus 饱和蒸气压、ASHRAE 湿球温度）、湿球能量上限下的恒速段变体与纯能量下界、平衡含水率驱动力 C_eq 的界定分析；`tab_extension`、`fig_extension`；求解器增加可选的蒸发通量上限（默认关闭，守恒不变量保持）；`tests/unit/test_extension.py`。
- `docs/mdr/0008-*.md`；`[a.extension]` 配置段；`package.stages` 登记 `extension`。

## [0.1.0] - 2026-09-12

### Added
- Forge 云端运行时（阶段/运行/清单模型、数据与结果契约、事件日志、确定性打包）。
- 通用流水线：ingest、validate、fmt、lint、test、paper、qa、package、release；全部在 Modal 执行。
- CUMCM 格式论文骨架（XeLaTeX + BibTeX gbt7714），自动生成数字宏、支撑材料清单、复现记录与完整程序附录。
- 工程规范、建模与写作标准、运行手册、ADR/MDR、AI 使用记录、Unlicense。

### Verified
- 三题仓库全链路冒烟通过：25 项单元/性质/集成测试，13 项论文与结果 QA 检查，发布散列核对一致。
