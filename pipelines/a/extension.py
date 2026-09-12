"""Model-evaluation extension (MDR-0008): the physics the problem's stylised model leaves out, quantified.

* a-posteriori latent-heat diagnostic of the problem-3 solution: the latent-heat demand of the model's own
  evaporation rate, expressed as the equivalent surface cooling L_v j_w / h, against the convective supply;
* psychrometric reading of the chamber state (relative humidity, wet-bulb temperature);
* energy-limited variant: the evaporation flux capped by the wet-bulb limit h (T_air - T_wb)/L_v (the
  constant-rate period of classical drying theory) and the pure energy lower bound on the drying time;
* equilibrium-moisture driving force: C_air replaced by a constant C_eq in the mass boundary condition.

The main results of problems 1-4 are untouched; everything here feeds the paper's model-evaluation section.
"""

from __future__ import annotations

from dataclasses import replace
from typing import Any

import numpy as np

from forge.context import StageContext
from forge.runner import stage
from pipelines.a.common import build_spec, production_n, spec_to_recipe, time_options
from pipelines.a.physics import (
    C_INIT,
    C_TARGET,
    R0,
    FluxCap,
    dry_solid_density,
    latent_heat,
    relative_humidity,
    wet_bulb_temperature,
)
from pipelines.a.science import HOUR, _parallel
from pipelines.a.solver import ProblemSpec

ORDINALS = ["One", "Two", "Three", "Four", "Five", "Six"]


def slug(label: str) -> str:
    """File-system and npz friendly variant key: 'cap:1.0' -> 'cap_1p0'."""
    return label.replace(":", "_").replace(".", "p")


def _first_time_below(t: np.ndarray, values: np.ndarray, threshold: float) -> float | None:
    idx = np.where(values < threshold)[0]
    return None if len(idx) == 0 else float(t[idx[0]])


def latent_diagnostic(hist: dict[str, Any], spec: ProblemSpec, rho_s: float) -> dict[str, Any]:
    """Latent-heat demand implied by the model's own evaporation rate versus its convective heat supply.

    j_w = rho_s hm (C_s - C_air) is the evaporation mass flux behind the moisture boundary condition; the energy
    equation of the problem ignores its latent heat L_v j_w. The equivalent surface cooling L_v j_w / h is the
    temperature depression a coupled model would need to carry that load by convection alone.
    """
    t = np.asarray(hist["t"], dtype=float)
    ts = np.asarray(hist["surface_temp"], dtype=float)
    ta = np.asarray(hist["air_temp"], dtype=float)
    cs = np.asarray(hist["surface_moist"], dtype=float)
    ca = np.asarray(hist["air_hum"], dtype=float)
    mean_c = np.asarray(hist["mean_moist"], dtype=float)
    flux = np.asarray(spec.surface_moisture_loss(t, cs, ca), dtype=float)  # (kg/kg) m/s
    j_w = rho_s * flux  # kg/(m^2 s)
    q_lat = latent_heat(ts) * j_w  # W/m^2
    q_conv = spec.props.h * (ta - ts)  # W/m^2 delivered by the model's thermal boundary condition
    dt_lat = q_lat / spec.props.h  # K
    circ = 2.0 * np.pi * R0
    e_lat = np.concatenate([[0.0], np.cumsum(0.5 * (q_lat[1:] + q_lat[:-1]) * np.diff(t))]) * circ  # J/m
    e_conv = np.concatenate([[0.0], np.cumsum(0.5 * (q_conv[1:] + q_conv[:-1]) * np.diff(t))]) * circ
    water_removed = rho_s * np.pi * R0**2 * (mean_c[0] - mean_c)  # kg per m length
    at: dict[str, dict[str, float]] = {}
    for label, hours in (("three_h", 3.0), ("twelve_h", 12.0), ("twenty_four_h", 24.0)):
        i = int(np.argmin(np.abs(t - hours * HOUR)))
        at[label] = {
            "t_h": float(t[i] / HOUR),
            "dt_lat_K": float(dt_lat[i]),
            "q_lat_W_m2": float(q_lat[i]),
            "q_conv_W_m2": float(q_conv[i]),
            "e_lat_MJ_m": float(e_lat[i] / 1e6),
            "e_conv_kJ_m": float(e_conv[i] / 1e3),
        }
    below_1k = _first_time_below(t, dt_lat, 1.0)
    below_01k = _first_time_below(t, dt_lat, 0.1)
    return {
        "arrays": {
            "t": t,
            "dt_lat": dt_lat,
            "model_dt": ta - ts,
            "q_lat": q_lat,
            "q_conv": q_conv,
            "e_lat": e_lat,
            "e_conv": e_conv,
            "water_removed": water_removed,
        },
        "summary": {
            "rho_s0": rho_s,
            "j_w_initial": float(j_w[0]),
            "dt_lat_peak_K": float(dt_lat.max()),
            "dt_lat_peak_t_h": float(t[int(np.argmax(dt_lat))] / HOUR),
            "dt_lat_below_1K_h": None if below_1k is None else below_1k / HOUR,
            "dt_lat_below_0p1K_h": None if below_01k is None else below_01k / HOUR,
            "at": at,
            "e_lat_total_MJ_per_m": float(e_lat[-1] / 1e6),
            "e_conv_total_kJ_per_m": float(e_conv[-1] / 1e3),
            "latent_to_convective_ratio": float(e_lat[-1] / max(e_conv[-1], 1e-300)),
            "water_removed_kg_per_m": float(water_removed[-1]),
        },
    }


def chamber_psychrometrics(spec: ProblemSpec, cap: FluxCap) -> dict[str, Any]:
    ch = spec.chamber
    t_plateau = ch.temp_plateau * ch.temp_scale
    w_plateau = ch.hum_plateau * ch.hum_scale
    wb_plateau = wet_bulb_temperature(t_plateau, w_plateau)
    t0, w0 = float(ch.temp[0]), float(ch.hum_scale * ch.hum[0])
    wb0 = wet_bulb_temperature(t0, w0)
    return {
        "assumption": "chamber moisture concentration read as the humidity ratio of the air (kg water per kg dry air)"
        " at 101325 Pa; Magnus saturation pressure; ASHRAE adiabatic-saturation wet bulb",
        "initial": {
            "temp": t0,
            "w": w0,
            "rh": float(relative_humidity(t0, w0)),
            "wet_bulb": wb0,
            "wet_bulb_depression": t0 - wb0,
            "flux_cap": float(cap.cap[0]),
        },
        "plateau": {
            "temp": t_plateau,
            "w": w_plateau,
            "rh": float(relative_humidity(t_plateau, w_plateau)),
            "wet_bulb": wb_plateau,
            "wet_bulb_depression": t_plateau - wb_plateau,
            "flux_cap": cap.plateau,
            "latent_heat": float(latent_heat(wb_plateau)),
        },
    }


@stage(
    "extension",
    deps=("ingest", "q3"),
    description="Model-evaluation extension: latent-heat consistency, wet-bulb energy limit, C_eq driving force",
)
def extension(ctx: StageContext) -> dict[str, Any]:
    n = production_n(ctx)
    opts = time_options(ctx)
    horizon = float(ctx.cfg("a.drying.horizon_s", 432000))
    cap_factors = [float(f) for f in ctx.cfg("a.extension.cap_factors", [0.5, 1.0, 2.0])]
    c_eq_values = [float(c) for c in ctx.cfg("a.extension.c_eq", [0.075, 0.10, 0.125])]
    base = build_spec(ctx, "q3")
    chamber, props = base.chamber, base.props
    rho_s0 = float(dry_solid_density(props, C_INIT))
    q3 = ctx.dep_data("q3", "summary.json")
    t3 = float(q3["t_dry_s"])
    hist = ctx.dep_data("q3", "history.json")

    # 1. a-posteriori latent-heat diagnostic of the problem-3 solution
    diag = latent_diagnostic(hist, base, rho_s0)
    np.savez_compressed(ctx.out("diagnostic.npz"), **diag["arrays"])
    ctx.log.info("extension.diagnostic", **{k: v for k, v in diag["summary"].items() if not isinstance(v, dict)})

    # 2. psychrometrics of the chamber and the wet-bulb energy limit on the evaporation flux
    cap = FluxCap.wet_bulb(chamber, props.h, rho_s0)
    psychro = chamber_psychrometrics(base, cap)
    mean_end = float(q3["mean_moist_at_dry"])
    t_bound = rho_s0 * R0 * (C_INIT - mean_end) / (2.0 * cap.plateau)  # lateral surface only, plateau limit
    ctx.log.info("extension.psychrometrics", **psychro["plateau"], energy_bound_h=t_bound / HOUR)

    # 3. variants: energy-capped evaporation (x factors) and equilibrium-moisture driving force
    jobs: list[dict[str, Any]] = []
    for f in cap_factors:
        spec = replace(base, flux_cap=replace(cap, scale=f))
        jobs.append(
            {
                "label": f"cap:{f}",
                "recipe": spec_to_recipe(spec),
                "n": n,
                "opts": opts,
                "mode": "dry",
                "horizon": horizon,
                "checks": True,
            }
        )
    local_spec = replace(base, flux_cap=replace(cap, scale=1.0, local_rho_s=True))
    jobs.append(
        {
            "label": "caplocal:1.0",
            "recipe": spec_to_recipe(local_spec),
            "n": n,
            "opts": opts,
            "mode": "dry",
            "horizon": horizon,
            "checks": True,
        }
    )
    for c in c_eq_values:
        spec = replace(base, chamber=replace(chamber, hum=np.full_like(chamber.hum, c), hum_plateau=c, hum_scale=1.0))
        jobs.append(
            {
                "label": f"ceq:{c}",
                "recipe": spec_to_recipe(spec),
                "n": n,
                "opts": opts,
                "mode": "dry",
                "horizon": horizon,
                "checks": True,
            }
        )
    results = {r["label"]: r for r in _parallel(ctx, jobs)}

    variants: list[dict[str, Any]] = []
    arrays: dict[str, np.ndarray] = {}
    checks: dict[str, Any] = {}
    for label, res in results.items():
        kind, raw = label.split(":")
        value = float(raw)
        t_dry = float(res["t_dry"])
        row: dict[str, Any] = {
            "case": kind,
            "value": value,
            "label": label,
            "slug": slug(label),
            "t_dry_h": t_dry / HOUR,
            "t_dry_grid_h": float(res["t_dry_grid"]) / HOUR,
            "diff_pct": 100.0 * (t_dry / t3 - 1.0),
            "checks_passed": bool(res["checks"]["passed"]),
            "moisture_residual": float(res["checks"]["moisture_balance"]["max_relative_residual"]),
            "nfev": int(res["stats"]["nfev"]),
        }
        tt = np.asarray(res["mean_moist_t"], dtype=float)
        if kind in ("cap", "caplocal"):
            scaled = replace(cap, scale=value)
            cs = np.asarray(res["surface_moist"], dtype=float)
            ca = np.asarray(res["air_hum"], dtype=float)
            rho_s = np.asarray(dry_solid_density(props, cs), dtype=float) if kind == "caplocal" else rho_s0
            active = props.hm * (cs - ca) > np.asarray(scaled.value(tt), dtype=float) / rho_s * (1.0 + 1e-9)
            row["cap_release_h"] = float(tt[np.where(active)[0][-1]] / HOUR) if active.any() else 0.0
            row["cap_plateau_kg_m2_s"] = cap.plateau * value
            row["cap_plateau_kg_m2_h"] = cap.plateau * value * HOUR
        checks[label] = res["checks"]
        variants.append(row)
        for key in ("mean_moist", "centre_moist", "surface_moist"):
            arrays[f"{slug(label)}__{key}"] = np.asarray(res[key], dtype=float)
        arrays[f"{slug(label)}__t"] = tt
    np.savez_compressed(ctx.out("variants.npz"), **arrays)
    checks["passed"] = bool(all(v["checks_passed"] for v in variants))
    ctx.write_json("verification_extension.json", checks)
    report = {
        "base_t_dry_h": t3 / HOUR,
        "rho_s0": rho_s0,
        "psychrometrics": psychro,
        "diagnostic": diag["summary"],
        "energy_bound_h": t_bound / HOUR,
        "cap_factors": cap_factors,
        "c_eq": c_eq_values,
        "variants": variants,
        "n": n,
    }
    ctx.write_json("extension.json", report)

    # 4. paper numbers
    d = diag["summary"]
    ctx.number("ExtDrySolidDensity", rho_s0, ".0f")
    ctx.number("ExtRhInitPct", 100.0 * psychro["initial"]["rh"], ".0f")
    ctx.number("ExtRhPlateauPct", 100.0 * psychro["plateau"]["rh"], ".0f")
    ctx.number("ExtWetBulbInit", psychro["initial"]["wet_bulb"], ".1f")
    ctx.number("ExtWetBulbPlateau", psychro["plateau"]["wet_bulb"], ".1f")
    ctx.number("ExtWetBulbDepression", psychro["plateau"]["wet_bulb_depression"], ".1f")
    ctx.number("ExtLatentHeatMJ", psychro["plateau"]["latent_heat"] / 1e6, ".3f")
    ctx.number("ExtFluxCapPlateau", cap.plateau, ".2e")
    ctx.number("ExtFluxCapPlateauKgPerHour", cap.plateau * HOUR, ".3f")
    ctx.number("ExtFluxCapInit", float(cap.cap[0]), ".2e")
    ctx.number("ExtFluxCapInitKgPerHour", float(cap.cap[0]) * HOUR, ".3f")
    ctx.number("ExtModelFluxInitKgPerHour", d["j_w_initial"] * HOUR, ".2f")
    ctx.number("ExtFluxRatioPlateau", d["j_w_initial"] / cap.plateau, ".1f")
    ctx.number("ExtDrySolidDensityTarget", float(dry_solid_density(props, np.array(C_TARGET))), ".0f")
    ctx.number("ExtModelFluxInit", d["j_w_initial"], ".2e")
    ctx.number("ExtFluxRatioInit", d["j_w_initial"] / float(cap.cap[0]), ".0f")
    ctx.number("ExtLatentDeltaTPeak", d["dt_lat_peak_K"], ".0f")
    ctx.number("ExtLatentDeltaTThreeH", d["at"]["three_h"]["dt_lat_K"], ".1f")
    ctx.number("ExtLatentDeltaTTwelveH", d["at"]["twelve_h"]["dt_lat_K"], ".1f")
    ctx.number("ExtLatentDeltaTTwentyFourH", d["at"]["twenty_four_h"]["dt_lat_K"], ".2f")
    if d["dt_lat_below_1K_h"] is not None:
        ctx.number("ExtLatentBelowOneKHours", d["dt_lat_below_1K_h"], ".1f")
    if d["dt_lat_below_0p1K_h"] is not None:
        ctx.number("ExtLatentBelowTenthKHours", d["dt_lat_below_0p1K_h"], ".1f")
    ctx.number("ExtLatentEnergyMJ", d["e_lat_total_MJ_per_m"], ".2f")
    ctx.number("ExtConvEnergyKJ", d["e_conv_total_kJ_per_m"], ".1f")
    ctx.number("ExtLatentToConvRatio", d["latent_to_convective_ratio"], ".0f")
    ctx.number("ExtWaterRemovedKg", d["water_removed_kg_per_m"], ".3f")
    ctx.number("ExtEnergyBoundHours", t_bound / HOUR, ".1f")
    ctx.number("ExtEnergyBoundPct", 100.0 * (t_bound / t3 - 1.0), ".1f")
    for row in [v for v in variants if v["case"] == "caplocal"]:
        ctx.number("ExtCapLocalHours", row["t_dry_h"], ".2f")
        ctx.number("ExtCapLocalPct", row["diff_pct"], ".1f")
        ctx.number("ExtCapLocalReleaseHours", row["cap_release_h"], ".1f")
    cap_rows = [v for v in variants if v["case"] == "cap"]
    ceq_rows = [v for v in variants if v["case"] == "ceq"]
    for i, row in enumerate(cap_rows):
        tag = ORDINALS[i]
        ctx.number(f"ExtCap{tag}Factor", row["value"], ".1f")
        ctx.number(f"ExtCap{tag}Hours", row["t_dry_h"], ".2f")
        ctx.number(f"ExtCap{tag}Pct", row["diff_pct"], ".1f")
        ctx.number(f"ExtCap{tag}ReleaseHours", row["cap_release_h"], ".1f")
        if row["value"] == 1.0:
            ctx.number("ExtCappedHours", row["t_dry_h"], ".2f")
            ctx.number("ExtCappedPct", row["diff_pct"], ".1f")
            ctx.number("ExtCappedReleaseHours", row["cap_release_h"], ".1f")
    for i, row in enumerate(ceq_rows):
        tag = ORDINALS[i]
        ctx.number(f"ExtCeq{tag}Value", row["value"], ".3f")
        ctx.number(f"ExtCeq{tag}Hours", row["t_dry_h"], ".2f")
        ctx.number(f"ExtCeq{tag}Pct", row["diff_pct"], ".1f")
    return {
        "jobs": len(jobs),
        "checks_passed": checks["passed"],
        "energy_bound_h": t_bound / HOUR,
        "variants": {v["label"]: v["t_dry_h"] for v in variants},
    }
