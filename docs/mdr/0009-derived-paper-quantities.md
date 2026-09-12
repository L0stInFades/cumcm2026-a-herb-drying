# MDR-0009 论文引用的派生量：无量纲数、时间尺度、Richardson 误差界与求解器统计

日期：2026-09-12 · 状态：已采纳 · 关联问题：问题 1–4（论文写作）

## 问题
工程规范要求论文中的每个数字来自 `ctx.number` 登记或生成表格，禁止手工算术。论文的分析、讨论与检验部分需要一批
由主结果派生的量：Biot/Fourier 数、热与质扩散时间尺度、扩散长度、轴向 Fourier 数（一维化论证）、扩散系数在干燥
前后的量级比（刚性论证）、烘干时间的 Richardson 外推误差界、温度与烘房平衡的时刻、求解器统计（右端项调用、LU
分解与重启段数）以及收缩几何比。这些量此前没有登记。

## 备选方案
1. 在 q1–q4、convergence 阶段追加 `ctx.number` 并以 `--force` 重跑这些阶段：重算耗时约 3 min，且会改写已核验阶段
   的清单。
2. 新增独立的轻量阶段 `derived`（依赖 q1、q2、q3、q4、convergence 的已存产物），只读取 `summary.json`、
   `history.json`、`convergence.json` 并计算派生量；不做新的模拟。
3. 在论文中手写这些数字：违反规范，且不可追溯。

## 决策与理由
采用方案 2。定义如下（全部在 `pipelines/a/derived.py`）：
- 热/质扩散时间尺度 $R_0^2/\alpha$、$R_0^2/D_0$，Fourier 数 $\mathrm{Fo}=\alpha t/R_0^2$（$t=1800$ s），扩散长度
  $\sqrt{D_0 t}$，刚性比 $\alpha/D_0$（附录 2 物性，初始状态）。
- 一维化论证的轴向 Fourier 数：$\alpha_3(C_0)\,t_{\mathrm{end}}/(L/2)^2$ 与上界 $D_{\max}t_{\mathrm{end}}/(L/2)^2$，
  $D_{\max}=D_3(C_0,T_{\mathrm{plateau}})$ 为全过程扩散系数的最大值（保守估计）。
- 刚性论证：$D_3$ 在 $(C_0,T_{\mathrm{plateau}})$、$(0.15,T_{\mathrm{plateau}})$ 与结束时表面状态下的值及其比值。
- 烘干时间离散误差界：以观测阶 $p=2$（收敛表证实）对最细两级做 Richardson 外推，
  $t^*=t_6+(t_6-t_5)/3$，生产网格误差 $|t_5-t^*|$（分钟）。
- 温度平衡时刻：$|T(0,t)-T_{\mathrm{air}}(t)|$ 此后始终小于 1 K（及 0.5 K）的首个时刻（问题 2，3 h 内）。
- 求解器统计取自各阶段 `summary.json` 的 `stats`；收缩几何比 $(R_{\mathrm{end}}/R_0)^2$ 与
  $(R_0/R_{\mathrm{end}})^2$ 取附件 2 末值。

## 对结果的影响与验证
不改变任何主结果；`derived` 登记的键仅供论文正文引用。辅助函数（Fourier 数、Richardson 外推、"此后始终小于"判据）
由 `tests/unit/test_derived.py` 覆盖；阶段加入 `package.stages`。
