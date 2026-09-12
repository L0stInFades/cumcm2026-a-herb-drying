"""Paper figures (vector PDF + PNG preview) built only from stage outputs.

Style: house palette (Okabe-Ito, fixed order for categorical series), single-hue sequential ramps for
time-ordered profiles (light = early, dark = late), one measure per axis, thin marks, recessive grid.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from forge import plotting
from forge.context import StageContext
from forge.runner import stage
from pipelines.a.physics import APPENDIX2, APPENDIX3, APPENDIX4, C_TARGET, Radius

HOUR = 3600.0
INK = "#1F2933"
MUTED = "#6B7280"
CAT = plotting.OKABE_ITO


def _ramp(name: str, count: int) -> list[Any]:
    from matplotlib import colormaps

    cmap = colormaps[name]
    return [cmap(x) for x in np.linspace(0.35, 0.95, count)]


def _load(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def _finish(ctx: StageContext, fig: Any, name: str, index: list[dict[str, str]], purpose: str) -> None:
    info = plotting.save(fig, ctx.out("figures", f"{name}.pdf"))
    index.append({"file": f"{name}.pdf", "purpose": purpose, "bytes": str(info["bytes"])})
    ctx.log.info("figure", name=name, bytes=info["bytes"])


def _profiles(ax: Any, snaps: Any, field: str, unit: str, cmap: str, label_fmt: Any) -> None:
    t = snaps["t"]
    colors = _ramp(cmap, len(t))
    for i, ti in enumerate(t):
        r_cm = snaps["xi"] * snaps["radius"][i] * 100.0
        ax.plot(r_cm, snaps[field][i], color=colors[i], label=label_fmt(ti), lw=1.4)
    ax.set_xlabel("到药材中心的距离 $r$ / cm")
    ax.set_ylabel(unit)
    ax.legend(frameon=False, ncol=2, title="时间", title_fontsize=8)


def _history(ax: Any, hist: dict[str, Any], keys: list[tuple[str, str]], unit: str, scale: float = 1.0) -> None:
    t = np.asarray(hist["t"]) / scale
    for j, (key, label) in enumerate(keys):
        style = "--" if key.startswith("air") else "-"
        ax.plot(t, hist[key], style, color=CAT[j], label=label, lw=1.4)
    ax.set_ylabel(unit)
    ax.legend(frameon=False)


@stage(
    "figures",
    deps=("ingest", "q1", "q2", "q3", "q4", "convergence", "verification", "sensitivity", "extension"),
    description="Paper figures: data, properties, profiles/histories of q1-q4, convergence, verification, sensitivity",
)
def figures(ctx: StageContext) -> dict[str, Any]:
    plotting.setup()
    index: list[dict[str, str]] = []
    data = ctx.dep("ingest") / "data"
    a1 = pd.read_parquet(data / "附件1__Sheet1.parquet")
    a2 = pd.read_parquet(data / "附件2__Sheet1.parquet")

    # 1. chamber conditions and radius data ---------------------------------------------------------
    q3sum = _load(ctx.dep("q3") / "summary.json")
    fig, axes = plotting.new_figure(6.3, 5.4, nrows=3, sharex=False)
    t1 = a1["时间"].to_numpy(float) / HOUR
    axes[0].plot(t1, a1["温度"], color=CAT[1], lw=1.2, label="附件 1 样本")
    axes[0].axhline(
        q3sum["spec"]["chamber"]["temp_plateau"], color=MUTED, ls=":", lw=1, label="恒温阶段平台（末 1 h 均值）"
    )
    axes[0].set_ylabel("烘房温度 / °C")
    axes[0].legend(frameon=False, loc="lower right")
    axes[1].plot(t1, a1["水分浓度"], color=CAT[0], lw=1.2, label="附件 1 样本")
    axes[1].axhline(q3sum["spec"]["chamber"]["hum_plateau"], color=MUTED, ls=":", lw=1, label="平台（末 1 h 均值）")
    axes[1].set_ylabel("烘房水分浓度 / (kg/kg)")
    axes[1].set_xlabel("时间 / h")
    axes[1].legend(frameon=False, loc="lower right")
    t2 = a2["时间"].to_numpy(float) / HOUR
    hist4 = _load(ctx.dep("q4") / "history.json")
    radius = Radius.from_series(a2["时间"].to_numpy(float), a2["半径"].to_numpy(float))
    t_fine = np.linspace(0.0, float(a2["时间"].max()), 2000)
    axes[2].plot(t_fine / HOUR, np.asarray(radius.value(t_fine)) * 100.0, color=CAT[2], lw=1.2, label="PCHIP 插值 R(t)")
    axes[2].plot(
        t2[::4], a2["半径"].to_numpy(float)[::4], "o", ms=3, color=INK, mfc="white", label="附件 2 样本（每 2 h）"
    )
    axes[2].set_ylabel("半径 / cm")
    axes[2].set_xlabel("时间 / h")
    axes[2].legend(frameon=False)
    for ax, tag in zip(axes, "abc"):
        ax.set_title(f"({tag})", loc="left", fontsize=9)
    _finish(ctx, fig, "fig_chamber", index, "附件 1 烘房温湿度（a,b）与附件 2 半径及其 PCHIP 插值（c）；边界条件数据")

    # 2. property laws ------------------------------------------------------------------------------
    fig, axes = plotting.new_figure(6.3, 2.9, ncols=2)
    c = np.linspace(0.02, 2.6, 400)
    axes[0].plot(c, APPENDIX2.diffusivity(c, np.full_like(c, 28.0)), color=CAT[0], label="附录 2")
    for j, (props, name) in enumerate(((APPENDIX3, "附录 3"), (APPENDIX4, "附录 4"))):
        for temp, ls in ((28.0, "--"), (50.0, "-")):
            axes[0].plot(
                c, props.diffusivity(c, np.full_like(c, temp)), ls, color=CAT[j + 1], label=f"{name}, {temp:g} °C"
            )
    axes[0].set_yscale("log")
    axes[0].set_xlabel("水分浓度 $C$ / (kg/kg)")
    axes[0].set_ylabel("扩散系数 $D$ / (m$^2$/s)")
    axes[0].axvline(C_TARGET, color=MUTED, ls=":", lw=1)
    axes[0].legend(frameon=False, fontsize=7)
    for j, (props, name) in enumerate(((APPENDIX3, "附录 3"), (APPENDIX4, "附录 4"))):
        axes[1].plot(c, props.k(c) / (props.rho(c) * props.cp(c)) * 1e7, color=CAT[j + 1], label=name)
    axes[1].axhline(APPENDIX2.k0 / (APPENDIX2.rho0 * APPENDIX2.cp0) * 1e7, color=CAT[0], label="附录 2")
    axes[1].set_xlabel("水分浓度 $C$ / (kg/kg)")
    axes[1].set_ylabel(r"热扩散率 $\alpha$ / ($10^{-7}$ m$^2$/s)")
    axes[1].legend(frameon=False, fontsize=7)
    for ax, tag in zip(axes, "ab"):
        ax.set_title(f"({tag})", loc="left", fontsize=9)
    _finish(
        ctx,
        fig,
        "fig_properties",
        index,
        "三组附录物性：水分扩散系数随 C、T 的变化（a，对数坐标）与热扩散率（b）；说明干燥后期表面干层的刚性",
    )

    # 3-4. problem 1 ---------------------------------------------------------------------------------
    snaps = np.load(ctx.dep("q1") / "snapshots.npz")
    fig, axes = plotting.new_figure(6.3, 2.9, ncols=2)
    _profiles(axes[0], snaps, "temp", "温度 / °C", "Oranges", lambda t: f"{t:.0f} s")
    _profiles(axes[1], snaps, "moist", "水分浓度 / (kg/kg)", "Blues", lambda t: f"{t:.0f} s")
    for ax, tag in zip(axes, "ab"):
        ax.set_title(f"({tag})", loc="left", fontsize=9)
    _finish(ctx, fig, "fig_q1_profiles", index, "问题 1：表 1/2 各时刻的温度（a）与水分浓度（b）径向剖面")
    hist = _load(ctx.dep("q1") / "history.json")
    fig, axes = plotting.new_figure(6.3, 2.9, ncols=2)
    _history(
        axes[0], hist, [("centre_temp", "中心 r=0"), ("surface_temp", "表面 r=2 cm"), ("air_temp", "烘房")], "温度 / °C"
    )
    _history(
        axes[1],
        hist,
        [("centre_moist", "中心 r=0"), ("surface_moist", "表面 r=2 cm"), ("mean_moist", "体积平均")],
        "水分浓度 / (kg/kg)",
    )
    for ax, tag in zip(axes, "ab"):
        ax.set_xlabel("时间 / s")
        ax.set_title(f"({tag})", loc="left", fontsize=9)
    _finish(
        ctx,
        fig,
        "fig_q1_history",
        index,
        "问题 1：中心、表面与烘房温度（a）及中心、表面、平均水分浓度（b）随时间的变化",
    )

    # 5-6. problem 2 ---------------------------------------------------------------------------------
    snaps = np.load(ctx.dep("q2") / "snapshots.npz")
    fig, axes = plotting.new_figure(6.3, 2.9, ncols=2)
    _profiles(axes[0], snaps, "temp", "温度 / °C", "Oranges", lambda t: f"{t / HOUR:.1f} h")
    _profiles(axes[1], snaps, "moist", "水分浓度 / (kg/kg)", "Blues", lambda t: f"{t / HOUR:.1f} h")
    for ax, tag in zip(axes, "ab"):
        ax.set_title(f"({tag})", loc="left", fontsize=9)
    _finish(ctx, fig, "fig_q2_profiles", index, "问题 2：表 3/4 各时刻的温度（a）与水分浓度（b）径向剖面")
    hist = _load(ctx.dep("q2") / "history.json")
    fig, axes = plotting.new_figure(6.3, 2.9, ncols=2)
    _history(
        axes[0],
        hist,
        [("centre_temp", "中心 r=0"), ("surface_temp", "表面 r=2 cm"), ("air_temp", "烘房")],
        "温度 / °C",
        HOUR,
    )
    _history(
        axes[1],
        hist,
        [("centre_moist", "中心 r=0"), ("surface_moist", "表面 r=2 cm"), ("mean_moist", "体积平均")],
        "水分浓度 / (kg/kg)",
        HOUR,
    )
    for ax, tag in zip(axes, "ab"):
        ax.set_xlabel("时间 / h")
        ax.set_title(f"({tag})", loc="left", fontsize=9)
    _finish(
        ctx, fig, "fig_q2_history", index, "问题 2：3 h 内温度（a）与水分浓度（b）的时程；表面在预热阶段即接近烘房温度"
    )

    # 7. problem 3 -----------------------------------------------------------------------------------
    snaps = np.load(ctx.dep("q3") / "snapshots.npz")
    hist = _load(ctx.dep("q3") / "history.json")
    fig, axes = plotting.new_figure(6.3, 2.9, ncols=2)
    _profiles(axes[0], snaps, "moist", "水分浓度 / (kg/kg)", "Blues", lambda t: f"{t / HOUR:.0f} h")
    axes[0].axhline(C_TARGET, color=MUTED, ls=":", lw=1)
    axes[0].annotate("0.15", xy=(0.02, C_TARGET), xytext=(0.02, C_TARGET + 0.08), color=MUTED, fontsize=7)
    _history(
        axes[1],
        hist,
        [("centre_moist", "中心 r=0"), ("surface_moist", "表面"), ("mean_moist", "体积平均")],
        "水分浓度 / (kg/kg)",
        HOUR,
    )
    axes[1].axhline(C_TARGET, color=MUTED, ls=":", lw=1)
    axes[1].axvline(q3sum["t_dry_h"], color=MUTED, ls="--", lw=1)
    axes[1].annotate(
        f"烘干结束 {q3sum['t_dry_h']:.2f} h",
        xy=(q3sum["t_dry_h"], 1.2),
        xytext=(q3sum["t_dry_h"] * 0.55, 1.6),
        color=INK,
        fontsize=7,
        arrowprops={"arrowstyle": "-", "color": MUTED, "lw": 0.8},
    )
    axes[1].set_xlabel("时间 / h")
    for ax, tag in zip(axes, "ab"):
        ax.set_title(f"({tag})", loc="left", fontsize=9)
    _finish(
        ctx,
        fig,
        "fig_q3",
        index,
        "问题 3：每 6 h 与烘干结束时刻的水分剖面（a）及中心、表面、平均水分浓度时程（b），虚线为 0.15 kg/kg 判据",
    )

    # 8. problem 4 -----------------------------------------------------------------------------------
    snaps = np.load(ctx.dep("q4") / "snapshots.npz")
    q4sum = _load(ctx.dep("q4") / "summary.json")
    fig, axes = plotting.new_figure(6.3, 2.9, ncols=2)
    _profiles(axes[0], snaps, "moist", "水分浓度 / (kg/kg)", "Blues", lambda t: f"{t / HOUR:.0f} h")
    axes[0].axhline(C_TARGET, color=MUTED, ls=":", lw=1)
    surf_r = snaps["radius"] * 100.0
    axes[0].plot(surf_r, snaps["moist"][:, -1], "o", ms=3, color=INK, mfc="white", label="药材表面", zorder=5)
    axes[0].legend(frameon=False, ncol=2, fontsize=7, title="时间", title_fontsize=8)
    t4 = np.asarray(hist4["t"]) / HOUR
    axes[1].plot(
        t4, hist4["centre_moist"], color=CAT[0], label=f"收缩，拉格朗日（主）: {q4sum['main']['t_dry_h']:.1f} h"
    )
    alt_hist = pd.read_parquet(ctx.dep("q4") / "profiles_moist_alt.parquet")
    axes[1].plot(
        alt_hist["t"].to_numpy(float) / HOUR,
        alt_hist["0.0"].to_numpy(float),
        "--",
        color=CAT[1],
        label=f"收缩，欧拉（备选）: {q4sum['alternative']['t_dry_h']:.1f} h",
    )
    fixed_t = q4sum["fixed_radius"]["t_dry_h"]
    axes[1].axhline(C_TARGET, color=MUTED, ls=":", lw=1)
    axes[1].axvline(fixed_t, color=CAT[2], ls="-.", lw=1, label=f"半径固定 2 cm: {fixed_t:.1f} h")
    axes[1].set_xlabel("时间 / h")
    axes[1].set_ylabel("中心水分浓度 / (kg/kg)")
    axes[1].legend(frameon=False, fontsize=7)
    for ax, tag in zip(axes, "ab"):
        ax.set_title(f"({tag})", loc="left", fontsize=9)
    _finish(
        ctx,
        fig,
        "fig_q4",
        index,
        "问题 4：物理坐标下每 6 h 的水分剖面与表面位置（a）；中心水分时程：主模型、欧拉备选、固定半径（b）",
    )

    # 9. convergence and analytic errors ------------------------------------------------------------
    conv = _load(ctx.dep("convergence") / "convergence.json")
    ver = _load(ctx.dep("verification") / "verification_report.json")
    fig, axes = plotting.new_figure(6.3, 2.9, ncols=2)
    keys = [
        ("q1_temp", "问题 1 温度"),
        ("q1_moist", "问题 1 水分"),
        ("q2_temp", "问题 2 温度"),
        ("q2_moist", "问题 2 水分"),
        ("q3_moist", "问题 3 水分"),
        ("q4_moist", "问题 4 水分"),
    ]
    markers = "osD^v<"
    for j, (key, label) in enumerate(keys):
        rows = conv["grid"][key]
        axes[0].loglog(
            [r["dr_cm"] for r in rows],
            [r["max_abs_error"] for r in rows],
            marker=markers[j],
            ms=4,
            color=CAT[j],
            label=label,
            lw=1.2,
        )
    dr = np.array([0.1 / 2**k for k in range(conv["levels"] - 1)])
    e0 = conv["grid"]["q1_temp"][0]["max_abs_error"]
    axes[0].loglog(dr, e0 * (dr / dr[0]) ** 2, color=MUTED, ls=":", lw=1, label="二阶参考")
    axes[0].set_xlabel(r"$\Delta r$ / cm")
    axes[0].set_ylabel("相对最细网格的最大误差")
    axes[0].legend(frameon=False, fontsize=6.5)
    for j, (field, label) in enumerate((("temp", "温度（贝塞尔级数）"), ("moist", "水分（贝塞尔级数）"))):
        rows = ver["analytic"][field]["levels"]
        axes[1].loglog(
            [r["dr_cm"] for r in rows],
            [r["max_abs_error"] for r in rows],
            marker=markers[j],
            ms=4,
            color=CAT[j],
            label=label,
            lw=1.2,
        )
    drv = np.array([r["dr_cm"] for r in ver["analytic"]["temp"]["levels"]])
    e0 = ver["analytic"]["temp"]["levels"][0]["max_abs_error"]
    axes[1].loglog(drv, e0 * (drv / drv[0]) ** 2, color=MUTED, ls=":", lw=1, label="二阶参考")
    axes[1].set_xlabel(r"$\Delta r$ / cm")
    axes[1].set_ylabel("与解析解的最大误差")
    axes[1].legend(frameon=False, fontsize=6.5)
    for ax, tag in zip(axes, "ab"):
        ax.set_title(f"({tag})", loc="left", fontsize=9)
    _finish(
        ctx,
        fig,
        "fig_convergence",
        index,
        "网格收敛：各问题相对最细网格的误差（a）与常物性算例相对贝塞尔级数解的误差（b），虚线为二阶斜率",
    )

    # 10. conservation residual histories ----------------------------------------------------------
    fig, ax = plotting.new_figure(6.3, 2.6)
    for j, prob in enumerate(("q1", "q2", "q3", "q4")):
        h = _load(ctx.dep(prob) / "history.json")
        t = np.asarray(h["t"])
        w = np.asarray(h["mean_moist"]) / 2.0
        rad = np.asarray(h["radius_cm"]) / 100.0
        cum = np.asarray(h["cum_moist_loss"])  # int hm (C_s - C_air)/R dt accumulated inside the integrator
        if np.allclose(rad, rad[0]):  # d/dt (R^2 W) = -R hm (C_s - C_air)
            lhs = rad[0] ** 2 * (w - w[0])
            cum = -(rad[0] ** 2) * cum
        else:  # Lagrangian: dW/dt = -hm (C_s - C_air)/R
            lhs = w - w[0]
            cum = -cum
        res = np.abs(lhs - cum) / max(abs(lhs[-1]), 1e-300)
        ax.semilogy(t / HOUR, np.maximum(res, 1e-17), color=CAT[j], label=f"问题 {prob[1]}", lw=1.2)
    ax.set_xlabel("时间 / h")
    ax.set_ylabel("水分守恒相对残差")
    ax.legend(frameon=False)
    _finish(
        ctx,
        fig,
        "fig_conservation",
        index,
        "离散水分总量变化与积分器内累积的表面通量之差（相对最终失水量）：四个问题均在机器精度量级",
    )

    # 11. sensitivity ------------------------------------------------------------------------------
    sens = _load(ctx.dep("sensitivity") / "sensitivity.json")
    names = {
        "h": "$h$",
        "hm": "$h_m$",
        "D0": "$D_0$ (问题 3)",
        "Cair": "$C_{\\mathrm{air}}$",
        "Tplateau": "$T_{\\mathrm{air}}$ 平台",
        "shrink": "收缩幅度 $s$ (问题 4)",
        "D0q4": "$D_0$ (问题 4)",
    }
    rows = [r for r in sens["table"] if r.get("elasticity") is not None]
    fig, ax = plotting.new_figure(6.3, 2.6)
    labels = [names.get(r["parameter"], r["parameter"]) for r in rows]
    vals = [r["elasticity"] for r in rows]
    colors = [CAT[1] if v > 0 else CAT[0] for v in vals]
    ax.barh(labels, vals, color=colors, height=0.55)
    for y, v in enumerate(vals):
        ax.text(
            v + (0.02 if v >= 0 else -0.02),
            y,
            f"{v:+.2f}",
            va="center",
            ha="left" if v >= 0 else "right",
            fontsize=7,
            color=INK,
        )
    ax.axvline(0, color=INK, lw=0.8)
    ax.set_xlabel("烘干时间的弹性 $E$")
    lim = max(abs(v) for v in vals) * 1.35 + 0.05
    ax.set_xlim(-lim, lim)
    ax.invert_yaxis()
    _finish(ctx, fig, "fig_sensitivity", index, "烘干时间对各参数的弹性（±10% 中心差分）：正值表示参数增大使烘干变慢")

    # 12. end effects (2-D axisymmetric) -----------------------------------------------------------
    ax2 = np.load(ctx.dep("verification") / "axisym_snapshots.npz")
    ee = ver["end_effect_2d"]
    fig, axes = plotting.new_figure(6.3, 2.7, ncols=2)
    axes[0].plot(
        ee["axial_profile_at_1d_dry"]["z_cm"],
        ee["axial_profile_at_1d_dry"]["centre_line_moist"],
        color=CAT[0],
        label="二维轴对称，中心线 r=0",
    )
    axes[0].axhline(C_TARGET, color=MUTED, ls=":", lw=1)
    axes[0].set_xlabel("到中截面的轴向距离 $z$ / cm")
    axes[0].set_ylabel("水分浓度 / (kg/kg)")
    axes[0].legend(frameon=False, fontsize=7)
    t = ax2["t"] / HOUR
    axes[1].plot(t, ax2["moist_1d"][:, 0], color=CAT[0], label="一维模型，中心")
    axes[1].plot(t, ax2["moist_mid"][:, 0], "--", color=CAT[1], label="二维模型，中截面中心")
    axes[1].plot(t, ax2["moist_end"][:, 0], "-.", color=CAT[2], label="二维模型，端面中心")
    axes[1].axhline(C_TARGET, color=MUTED, ls=":", lw=1)
    axes[1].set_xlabel("时间 / h")
    axes[1].set_ylabel("水分浓度 / (kg/kg)")
    axes[1].legend(frameon=False, fontsize=7)
    for ax, tag in zip(axes, "ab"):
        ax.set_title(f"({tag})", loc="left", fontsize=9)
    _finish(
        ctx,
        fig,
        "fig_endeffect",
        index,
        "端部效应：一维烘干结束时刻二维模型中心线沿轴向的水分（a）；一维与二维（中截面、端面）中心水分时程（b）",
    )

    # 13. model-evaluation extension (MDR-0008) ---------------------------------------------------
    ext = _load(ctx.dep("extension") / "extension.json")
    dg = np.load(ctx.dep("extension") / "diagnostic.npz")
    va = np.load(ctx.dep("extension") / "variants.npz")
    hist3 = _load(ctx.dep("q3") / "history.json")
    fig, axes = plotting.new_figure(6.3, 2.9, ncols=2)
    t = dg["t"] / HOUR
    axes[0].semilogy(t, np.maximum(dg["dt_lat"], 1e-3), color=CAT[1], lw=1.4, label="潜热等效温降 $L_v j_w/h$")
    axes[0].semilogy(
        t, np.maximum(dg["model_dt"], 1e-3), color=CAT[0], lw=1.4, label="题给模型 $T_{\\mathrm{air}}-T_s$"
    )
    axes[0].axhline(1.0, color=MUTED, ls=":", lw=1)
    axes[0].set_ylim(1e-3, 2e2)
    axes[0].set_xlabel("时间 / h")
    axes[0].set_ylabel("温差 / K")
    axes[0].legend(frameon=False, fontsize=7)
    axes[1].plot(
        np.asarray(hist3["t"]) / HOUR,
        hist3["centre_moist"],
        color=CAT[0],
        lw=1.4,
        label=f"题给模型: {ext['base_t_dry_h']:.1f} h",
    )
    styles = {"cap": ("--", CAT[1]), "ceq": ("-.", CAT[2])}
    shown = [v for v in ext["variants"] if (v["case"] == "cap" and v["value"] == 1.0) or v["case"] == "ceq"]
    ceq_colors = [CAT[2], CAT[3], CAT[4], CAT[5]]
    k = 0
    for v in shown:
        ls, color = styles[v["case"]]
        if v["case"] == "ceq":
            color = ceq_colors[k % len(ceq_colors)]
            k += 1
            label = f"$C_{{\\mathrm{{eq}}}}={v['value']:g}$: {v['t_dry_h']:.1f} h"
        else:
            label = f"能量限制（湿球极限）: {v['t_dry_h']:.1f} h"
        axes[1].plot(
            va[f"{v['slug']}__t"] / HOUR, va[f"{v['slug']}__centre_moist"], ls, color=color, lw=1.3, label=label
        )
    axes[1].axhline(C_TARGET, color=MUTED, ls=":", lw=1)
    axes[1].set_xlabel("时间 / h")
    axes[1].set_ylabel("中心水分浓度 / (kg/kg)")
    axes[1].legend(frameon=False, fontsize=7)
    for ax, tag in zip(axes, "ab"):
        ax.set_title(f"({tag})", loc="left", fontsize=9)
    _finish(
        ctx,
        fig,
        "fig_extension",
        index,
        "模型评价扩展：题给模型蒸发率的潜热等效温降与模型实际温差（a）；能量限制与平衡含水率驱动力下的中心水分时程（b）",
    )

    ctx.write_json("figure_index.json", index)
    return {"figures": len(index)}
