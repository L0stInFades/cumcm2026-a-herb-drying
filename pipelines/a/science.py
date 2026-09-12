"""Science stages of Problem A: q1–q4 (solutions + independent checks), convergence, verification, sensitivity."""

from __future__ import annotations

import math
import re
from concurrent.futures import ProcessPoolExecutor
from dataclasses import replace
from typing import Any

import numpy as np
import pandas as pd

from forge.context import StageContext
from forge.runner import stage
from pipelines.a import verify
from pipelines.a.common import (
    COLS,
    PROPS,
    R_OUT_CM,
    R_TABLE_CM,
    build_spec,
    load_chamber,
    load_radius,
    mean_moisture,
    physical_frame,
    production_n,
    profiles_frame,
    snapshot_arrays,
    spec_from_recipe,
    spec_to_recipe,
    table_cells,
    time_options,
)
from pipelines.a.physics import (
    APPENDIX2,
    C_INIT,
    C_TARGET,
    R0,
    T_INIT,
    Chamber,
    Radius,
    biot_numbers,
)
from pipelines.a.solver import ProblemSpec, Solution, sample_physical, solve, solve_until_dry, subgrid_indices

Q1_TABLE_T = np.array([100, 300, 600, 900, 1200, 1500, 1800], dtype=float)
Q2_TABLE_T = np.arange(1, 7) * 1800.0
HOUR = 3600.0


# ----------------------------------------------------------------------------------------------
# helpers
# ----------------------------------------------------------------------------------------------


def _store_solution(ctx: StageContext, sol: Solution, spec: ProblemSpec, snapshot_times: np.ndarray) -> None:
    moving = not spec.radius.is_constant
    for field in ("temp", "moist"):
        frame = physical_frame(sol, field) if moving else profiles_frame(sol, field)
        frame.to_parquet(ctx.out(f"profiles_{field}.parquet"), index=False)
    np.savez_compressed(ctx.out("snapshots.npz"), **snapshot_arrays(sol, snapshot_times))
    history = {
        "t": sol.t.tolist(),
        "air_temp": sol.air_temp.tolist(),
        "air_hum": sol.air_hum.tolist(),
        "radius_cm": (sol.radius * 100.0).tolist(),
        "centre_temp": sol.temp[:, 0].tolist(),
        "surface_temp": sol.temp[:, -1].tolist(),
        "centre_moist": sol.moist[:, 0].tolist(),
        "surface_moist": sol.moist[:, -1].tolist(),
        "mean_moist": mean_moisture(sol).tolist(),
        "surface_moist_flux": sol.surface_moisture_flux(spec.props).tolist(),
        "cum_moist_loss": None if sol.cum_moist_loss is None else sol.cum_moist_loss.tolist(),
    }
    ctx.write_json("history.json", history)


def _checks(spec: ProblemSpec, sol: Solution) -> dict[str, Any]:
    if spec.mass_form == "conservative":
        assert sol.cum_moist_loss is not None
        balance = verify.conservative_moisture_balance(
            sol.t, sol.xi, sol.moist, sol.radius, spec.props, sol.cum_moist_loss
        )
    else:
        balance = verify.moisture_balance(
            sol.t,
            sol.xi,
            sol.moist,
            sol.radius,
            sol.air_hum,
            spec.props.hm,
            cum_loss=sol.cum_moist_loss,
            flux=spec.surface_moisture_loss(sol.t, sol.moist[:, -1], sol.air_hum),
        )
    extremum = verify.extremum_checks(sol.temp, sol.moist, sol.air_temp, sol.air_hum, spec.t_init, spec.c_init)
    report: dict[str, Any] = {"moisture_balance": balance, "extremum": extremum}
    if spec.props.rho1 == 0 and spec.props.cp1 == 0 and spec.props.k1 == 0 and spec.radius.is_constant:
        report["energy_balance"] = verify.energy_balance(
            sol.t, sol.xi, sol.temp, float(sol.radius[0]), sol.air_temp, spec.props, cum_loss=sol.cum_heat_loss
        )
    if sol.t_dry is not None and sol.dry_profile is not None:
        report["drying"] = verify.drying_checks(sol.t, sol.moist, sol.t_dry, sol.dry_profile)
    passed = [
        balance["max_relative_residual"] < 1e-4,
        extremum["temp_within_bounds"],
        extremum["moist_within_bounds"],
        extremum["moist_monotone_in_r"],
        extremum["moist_argmax_is_centre"],
    ]
    if "energy_balance" in report:
        passed.append(report["energy_balance"]["max_relative_residual"] < 1e-4)
    if "drying" in report:
        d = report["drying"]
        passed += [d["rows_before_all_wet"], d["rows_after_all_dry"], d["profile_max_error_at_t_dry"] < 1e-6]
    report["passed"] = bool(all(passed))
    return report


def _run_recipe(job: dict[str, Any]) -> dict[str, Any]:
    """Worker: solve one recipe and return compact arrays (used by the parallel studies)."""
    spec = spec_from_recipe(job["recipe"])
    n = int(job["n"])
    opts = dict(job["opts"])
    out: dict[str, Any] = {"label": job["label"], "n": n, "opts": opts}
    if job["mode"] == "dry":
        sol = solve_until_dry(spec, n, float(job.get("dt_out", 60.0)), horizon=float(job["horizon"]), **opts)
        idx = subgrid_indices(n)
        snap_t = np.asarray(job.get("snapshot_times", []), dtype=float)
        lookup = {float(t): i for i, t in enumerate(sol.t)}
        rows = [lookup[float(t)] for t in snap_t if float(t) in lookup]
        out.update(
            {
                "t_dry": sol.t_dry,
                "t_dry_grid": sol.stats.get("t_dry_grid"),
                "snapshot_t": sol.t[rows].tolist(),
                "snapshot_moist": sol.moist[rows][:, idx].tolist(),
                "dry_profile_sub": None if sol.dry_profile is None else sol.dry_profile[idx].tolist(),
                "mean_moist_t": sol.t.tolist(),
                "mean_moist": mean_moisture(sol).tolist(),
                "centre_moist": sol.moist[:, 0].tolist(),
                "surface_moist": sol.moist[:, -1].tolist(),
                "air_hum": sol.air_hum.tolist(),
                "stats": sol.stats,
            }
        )
        if job.get("checks"):
            out["checks"] = _checks(spec, sol)
    else:
        t_out = np.asarray(job["t_out"], dtype=float)
        sol = solve(spec, n, t_out, **opts)
        idx = subgrid_indices(n)
        out.update({"t": sol.t.tolist(), "temp": sol.temp[:, idx].tolist(), "moist": sol.moist[:, idx].tolist()})
        out["stats"] = sol.stats
    return out


def _parallel(ctx: StageContext, jobs: list[dict[str, Any]], workers: int | None = None) -> list[dict[str, Any]]:
    limit = int(ctx.cfg("a.sensitivity.workers", 8)) if workers is None else int(workers)
    workers = max(1, min(limit, len(jobs)))
    ctx.log.info("parallel.start", jobs=len(jobs), workers=workers)
    results: list[dict[str, Any]] = []
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for i, res in enumerate(pool.map(_run_recipe, jobs)):
            results.append(res)
            ctx.log.info("parallel.progress", done=i + 1, total=len(jobs), label=res["label"])
    return results


def _tag(text: str) -> str:
    """Letters-only CamelCase key fragment for ctx.number (digits are spelled out)."""
    out = text
    for old, new in (("q1", "Qone"), ("q2", "Qtwo"), ("q3", "Qthree"), ("q4", "Qfour"), ("D0", "Dzero")):
        out = out.replace(old, new)
    parts = re.split(r"[^A-Za-z]+", out)
    return "".join(p[:1].upper() + p[1:] for p in parts if p)


def _table_times_hours(t_dry: float, step_h: float = 6.0) -> np.ndarray:
    last = math.floor(t_dry / HOUR / step_h)
    return np.arange(1, last + 1) * step_h * HOUR


# ----------------------------------------------------------------------------------------------
# q1: preheating stage, constant properties (appendix 2)
# ----------------------------------------------------------------------------------------------


@stage("q1", deps=("ingest",), description="Problem 1: preheating stage (0-1800 s), constant properties + D(C)")
def q1(ctx: StageContext) -> dict[str, Any]:
    spec = build_spec(ctx, "q1")
    n = production_n(ctx)
    opts = time_options(ctx)
    t_out = np.arange(0.0, 1800.0 + 0.5, 1.0)
    sol = solve(spec, n, t_out, **opts)
    ctx.log.info("q1.solved", **sol.stats)
    _store_solution(ctx, sol, spec, Q1_TABLE_T)
    temp = profiles_frame(sol, "temp")
    moist = profiles_frame(sol, "moist")
    tables = {"temp": table_cells(temp, Q1_TABLE_T), "moist": table_cells(moist, Q1_TABLE_T)}
    ctx.write_json("tables.json", tables)
    checks = _checks(spec, sol)
    ctx.write_json("verification_q1.json", checks)
    bi = biot_numbers(APPENDIX2, C_INIT, T_INIT)
    ctx.write_json("summary.json", {"spec": spec.describe(), "n": n, "stats": sol.stats, "biot": bi})
    mean_c = mean_moisture(sol)
    ctx.number("QoneCentreTempEnd", float(sol.temp[-1, 0]))
    ctx.number("QoneSurfaceTempEnd", float(sol.temp[-1, -1]))
    ctx.number("QoneCentreMoistEnd", float(sol.moist[-1, 0]))
    ctx.number("QoneSurfaceMoistEnd", float(sol.moist[-1, -1]))
    ctx.number("QoneMeanMoistEnd", float(mean_c[-1]))
    ctx.number("QoneMoistLossPct", float(100.0 * (1.0 - mean_c[-1] / mean_c[0])), ".2f")
    ctx.number("QoneBiotHeat", bi["Bi_heat"], ".3f")
    ctx.number("QoneBiotMass", bi["Bi_mass"], ".3f")
    ctx.number("QoneAlpha", bi["alpha"], ".3e")
    ctx.number("QoneDinit", bi["D"], ".3e")
    ctx.number("QoneAirTempEnd", float(sol.air_temp[-1]), ".3f")
    ctx.number("QoneSurfaceTempLag", float(sol.air_temp[-1] - sol.temp[-1, -1]), ".3f")
    ctx.number("GridNodes", n + 1)
    ctx.number("GridDrCm", 2.0 / n, ".4f")
    ctx.number("TimeRtol", opts["rtol"], ".0e")
    ctx.number("TimeAtol", opts["atol"], ".0e")
    return {"n": n, "checks_passed": checks["passed"], **{k: v for k, v in sol.stats.items() if k != "message"}}


# ----------------------------------------------------------------------------------------------
# q2: whole process, variable properties (appendix 3), first 3 h at 1 s
# ----------------------------------------------------------------------------------------------


@stage("q2", deps=("ingest",), description="Problem 2: whole-process model with appendix-3 properties, 0-3 h at 1 s")
def q2(ctx: StageContext) -> dict[str, Any]:
    spec = build_spec(ctx, "q2")
    n = production_n(ctx)
    opts = time_options(ctx)
    t_out = np.arange(0.0, 3.0 * HOUR + 0.5, 1.0)
    sol = solve(spec, n, t_out, **opts)
    ctx.log.info("q2.solved", **sol.stats)
    _store_solution(ctx, sol, spec, Q2_TABLE_T)
    temp = profiles_frame(sol, "temp")
    moist = profiles_frame(sol, "moist")
    tables = {"temp": table_cells(temp, Q2_TABLE_T), "moist": table_cells(moist, Q2_TABLE_T)}
    ctx.write_json("tables.json", tables)
    checks = _checks(spec, sol)
    ctx.write_json("verification_q2.json", checks)
    bi0 = biot_numbers(spec.props, C_INIT, T_INIT)
    bi_end = biot_numbers(spec.props, float(sol.moist[-1, -1]), float(sol.temp[-1, -1]))
    ctx.write_json(
        "summary.json", {"spec": spec.describe(), "n": n, "stats": sol.stats, "biot_init": bi0, "biot_end": bi_end}
    )
    mean_c = mean_moisture(sol)
    ctx.number("QtwoCentreTempThreeH", float(sol.temp[-1, 0]))
    ctx.number("QtwoSurfaceTempThreeH", float(sol.temp[-1, -1]))
    ctx.number("QtwoCentreMoistThreeH", float(sol.moist[-1, 0]))
    ctx.number("QtwoSurfaceMoistThreeH", float(sol.moist[-1, -1]))
    ctx.number("QtwoMeanMoistThreeH", float(mean_c[-1]))
    ctx.number("QtwoMoistLossPct", float(100.0 * (1.0 - mean_c[-1] / mean_c[0])), ".2f")
    ctx.number("QtwoDinit", bi0["D"], ".3e")
    ctx.number("QtwoDsurfaceThreeH", bi_end["D"], ".3e")
    ctx.number("QtwoBiotMassThreeH", bi_end["Bi_mass"], ".2f")
    ctx.number("QtwoAlphaInit", bi0["alpha"], ".3e")
    return {"n": n, "checks_passed": checks["passed"], **{k: v for k, v in sol.stats.items() if k != "message"}}


# ----------------------------------------------------------------------------------------------
# q3: drying time (appendix 3, fixed radius)
# ----------------------------------------------------------------------------------------------


def _drying_outputs(
    ctx: StageContext, sol: Solution, spec: ProblemSpec, key: str, moving: bool
) -> tuple[dict[str, Any], list[list[Any]]]:
    assert sol.t_dry is not None and sol.dry_profile is not None
    t_dry = sol.t_dry
    rows_t = _table_times_hours(t_dry)
    frame = physical_frame(sol, "moist") if moving else profiles_frame(sol, "moist")
    table: list[list[Any]] = []
    if moving:
        lookup = {float(t): i for i, t in enumerate(sol.t)}
        for t, row in zip(rows_t, table_cells(frame, rows_t)):
            surface = float(frame.at[lookup[float(t)], "surface"])
            table.append([t / HOUR, *[None if np.isnan(v) else float(v) for v in row[1:]], surface])
        radius_dry = float(spec.radius.value(t_dry))
        end_vals = sample_physical(sol.dry_profile, sol.xi, radius_dry, R_TABLE_CM / 100.0)
        end_row = [t_dry / HOUR, *[None if np.isnan(v) else float(v) for v in end_vals], float(sol.dry_profile[-1])]
    else:
        table = [[t / HOUR, *row[1:]] for t, row in zip(rows_t, table_cells(frame, rows_t))]
        idx = subgrid_indices(len(sol.xi) - 1)
        prof = sol.dry_profile[idx]
        end_row = [t_dry / HOUR, *[float(prof[round(r / 0.1)]) for r in R_TABLE_CM]]
    table.append(end_row)
    mean_c = mean_moisture(sol)
    below = np.where(mean_c < C_TARGET)[0]
    half = np.where(mean_c <= 0.5 * mean_c[0])[0]
    summary = {
        "t_dry_s": t_dry,
        "t_dry_h": t_dry / HOUR,
        "t_dry_grid_s": sol.stats.get("t_dry_grid"),
        "t_dry_grid_h": float(sol.stats.get("t_dry_grid", 0.0)) / HOUR,
        "t_mean_below_h": None if len(below) == 0 else float(sol.t[below[0]] / HOUR),
        "t_half_loss_h": None if len(half) == 0 else float(sol.t[half[0]] / HOUR),
        "surface_moist_at_dry": float(sol.dry_profile[-1]),
        "mean_moist_at_dry": float(2.0 * (sol.dry_profile @ verify.control_volumes(sol.xi))),
        "rows": int(len(sol.t) - 1),
        "table_times_h": (rows_t / HOUR).tolist(),
    }
    ctx.write_json(f"drying_{key}.json", summary)
    return summary, table


@stage("q3", deps=("ingest",), description="Problem 3: drying time until every point is below 0.15 kg/kg")
def q3(ctx: StageContext) -> dict[str, Any]:
    spec = build_spec(ctx, "q3")
    n = production_n(ctx)
    opts = time_options(ctx)
    sol = solve_until_dry(
        spec, n, float(ctx.cfg("a.drying.dt_out_s", 60.0)), horizon=float(ctx.cfg("a.drying.horizon_s", 432000)), **opts
    )
    assert sol.t_dry is not None
    ctx.log.info("q3.solved", t_dry_h=sol.t_dry / HOUR, **sol.stats)
    snaps = np.concatenate([_table_times_hours(sol.t_dry), [sol.t[-1]]])
    _store_solution(ctx, sol, spec, snaps)
    np.savez_compressed(ctx.out("dry_profile.npz"), xi=sol.xi, moist=sol.dry_profile, temp=sol.dry_temp, t=sol.t_dry)
    summary, table = _drying_outputs(ctx, sol, spec, "q3", moving=False)
    ctx.write_json("tables.json", {"moist": table})
    checks = _checks(spec, sol)
    ctx.write_json("verification_q3.json", checks)
    ctx.write_json("summary.json", {"spec": spec.describe(), "n": n, "stats": sol.stats, **summary})
    ctx.number("QthreeDryHours", summary["t_dry_h"], ".2f")
    ctx.number("QthreeDryHoursFour", summary["t_dry_h"], ".4f")
    ctx.number("QthreeDryHoursGrid", summary["t_dry_grid_h"], ".4f")
    ctx.number("QthreeDrySeconds", summary["t_dry_s"], ".1f")
    ctx.number("QthreeDryDays", summary["t_dry_h"] / 24.0, ".2f")
    ctx.number("QthreeSurfaceMoistAtDry", summary["surface_moist_at_dry"])
    ctx.number("QthreeMeanMoistAtDry", summary["mean_moist_at_dry"])
    if summary["t_mean_below_h"] is not None:
        ctx.number("QthreeMeanCriterionHours", summary["t_mean_below_h"], ".2f")
    if summary["t_half_loss_h"] is not None:
        ctx.number("QthreeHalfLossHours", summary["t_half_loss_h"], ".2f")
    ctx.number("QthreeResultRows", summary["rows"])
    ctx.number("QthreeAirTempPlateau", spec.chamber.temp_plateau, ".3f")
    ctx.number("QthreeAirHumPlateau", spec.chamber.hum_plateau, ".5f")
    return {"n": n, "t_dry_h": summary["t_dry_h"], "checks_passed": checks["passed"], "nfev": sol.stats["nfev"]}


# ----------------------------------------------------------------------------------------------
# q4: shrinking radius (attachment 2), appendix-4 properties, Lagrangian main + Eulerian alternative
# ----------------------------------------------------------------------------------------------


@stage("q4", deps=("ingest",), description="Problem 4: shrinking cylinder (attachment 2) with appendix-4 properties")
def q4(ctx: StageContext) -> dict[str, Any]:
    n = production_n(ctx)
    opts = time_options(ctx)
    dt_out = float(ctx.cfg("a.drying.dt_out_s", 60.0))
    horizon = float(ctx.cfg("a.drying.horizon_s", 432000))
    main_form = str(ctx.cfg("a.scheme.formulation", "lagrangian"))
    alt_form = "eulerian" if main_form == "lagrangian" else "lagrangian"
    spec = build_spec(ctx, "q4", formulation=main_form)
    sol = solve_until_dry(spec, n, dt_out, horizon=horizon, **opts)
    assert sol.t_dry is not None and sol.dry_profile is not None
    ctx.log.info("q4.solved", formulation=main_form, t_dry_h=sol.t_dry / HOUR, **sol.stats)
    snaps = np.concatenate([_table_times_hours(sol.t_dry), [sol.t[-1]]])
    _store_solution(ctx, sol, spec, snaps)
    np.savez_compressed(
        ctx.out("dry_profile.npz"),
        xi=sol.xi,
        moist=sol.dry_profile,
        temp=sol.dry_temp,
        t=sol.t_dry,
        radius=float(spec.radius.value(sol.t_dry)),
    )
    summary, table = _drying_outputs(ctx, sol, spec, "q4", moving=True)
    ctx.write_json("tables.json", {"moist": table})
    checks = _checks(spec, sol)

    # alternative formulation (MDR-0006) and the fixed-radius reference with the same properties
    alt_spec = replace(spec, formulation=alt_form)
    alt = solve_until_dry(alt_spec, n, dt_out, horizon=horizon, **opts)
    assert alt.t_dry is not None
    ctx.log.info("q4.alternative", formulation=alt_form, t_dry_h=alt.t_dry / HOUR, **alt.stats)
    alt_summary, alt_table = _drying_outputs(ctx, alt, alt_spec, "q4_alt", moving=True)
    physical_frame(alt, "moist").to_parquet(ctx.out("profiles_moist_alt.parquet"), index=False)
    ctx.write_json("tables_alt.json", {"moist": alt_table, "formulation": alt_form})
    alt_checks = _checks(alt_spec, alt)
    fixed_spec = replace(spec, radius=Radius.constant(R0))
    fixed = solve_until_dry(fixed_spec, n, dt_out, horizon=horizon, **opts)
    assert fixed.t_dry is not None
    fixed_summary, _ = _drying_outputs(ctx, fixed, fixed_spec, "q4_fixed", moving=False)

    # degenerate check: both moving-boundary code paths with a (formally non-constant) constant radius
    pseudo = Radius(np.array([0.0, 1.0, 2.0]), np.array([R0, R0, R0]), force_interpolant=True)
    t_short = np.arange(0.0, 6.0 * HOUR + 0.5, 600.0)
    ref = solve(fixed_spec, n, t_short, **opts)
    degenerate: dict[str, Any] = {}
    for form in ("lagrangian", "eulerian"):
        test = solve(replace(spec, radius=pseudo, formulation=form), n, t_short, **opts)
        degenerate[form] = {
            "max_abs_diff_moist": float(np.max(np.abs(test.moist - ref.moist))),
            "max_abs_diff_temp": float(np.max(np.abs(test.temp - ref.temp))),
        }
    degenerate["passed"] = bool(
        all(v["max_abs_diff_moist"] < 1e-6 and v["max_abs_diff_temp"] < 1e-4 for v in degenerate.values())
    )
    checks["degenerate_constant_radius"] = degenerate
    checks["passed"] = bool(checks["passed"] and degenerate["passed"])
    checks["alternative_formulation"] = alt_checks
    ctx.write_json("verification_q4.json", checks)

    radius_dry = float(spec.radius.value(sol.t_dry)) * 100.0
    ctx.write_json(
        "summary.json",
        {
            "spec": spec.describe(),
            "n": n,
            "stats": sol.stats,
            "main": summary,
            "alternative": {"formulation": alt_form, **alt_summary},
            "fixed_radius": fixed_summary,
            "radius_at_dry_cm": radius_dry,
        },
    )
    ctx.number("QfourDryHours", summary["t_dry_h"], ".2f")
    ctx.number("QfourDryHoursFour", summary["t_dry_h"], ".4f")
    ctx.number("QfourDryHoursGrid", summary["t_dry_grid_h"], ".4f")
    ctx.number("QfourDrySeconds", summary["t_dry_s"], ".1f")
    ctx.number("QfourDryDays", summary["t_dry_h"] / 24.0, ".2f")
    ctx.number("QfourDryHoursAlt", alt_summary["t_dry_h"], ".2f")
    ctx.number("QfourDryHoursFixed", fixed_summary["t_dry_h"], ".2f")
    ctx.number("QfourAltDiffPct", 100.0 * (alt_summary["t_dry_h"] / summary["t_dry_h"] - 1.0), ".2f")
    ctx.number("QfourFixedDiffPct", 100.0 * (fixed_summary["t_dry_h"] / summary["t_dry_h"] - 1.0), ".2f")
    ctx.number("QfourRadiusAtDryCm", radius_dry, ".3f")
    ctx.number("QfourSurfaceMoistAtDry", summary["surface_moist_at_dry"])
    ctx.number("QfourMeanMoistAtDry", summary["mean_moist_at_dry"])
    ctx.number("QfourResultRows", summary["rows"])
    worst = max(v["max_abs_diff_moist"] for v in degenerate.values() if isinstance(v, dict))
    ctx.number("QfourDegenerateMaxDiff", worst, ".1e")
    return {
        "n": n,
        "t_dry_h": summary["t_dry_h"],
        "t_dry_alt_h": alt_summary["t_dry_h"],
        "t_dry_fixed_h": fixed_summary["t_dry_h"],
        "checks_passed": checks["passed"],
    }


# ----------------------------------------------------------------------------------------------
# convergence: grid refinement, time tolerance, integrator cross-check
# ----------------------------------------------------------------------------------------------


def _output_grid_study(ctx: StageContext, refine: int, opts: dict[str, Any]) -> dict[str, Any]:
    """Refinement study on the *full* result-file grid (every 1 s x 0.1 cm), not only the table instants.

    The table instants start at 100 s (problem 1) and 0.5 h (problem 2), but result1/result2.xlsx start at
    t = 1 s, where the surface boundary layer sqrt(D t) is thinner than a few production cells. This study
    therefore compares the *delivered* arrays against a Richardson limit built from two finer grids on every
    one of the 1800 x 21 (resp. 10800 x 21) output cells, and counts the cells whose fourth decimal differs.
    """
    levels = [refine + 1, refine + 2]
    grids = {"q1": np.arange(0.0, 1800.0 + 0.5, 1.0), "q2": np.arange(0.0, 3.0 * HOUR + 0.5, 1.0)}
    jobs = [
        {
            "label": f"{prob}:out:k{k}",
            "recipe": spec_to_recipe(build_spec(ctx, prob)),
            "n": 20 * 2**k,
            "opts": opts,
            "mode": "fixed",
            "t_out": t_out.tolist(),
        }
        for prob, t_out in grids.items()
        for k in levels
    ]
    results = {r["label"]: r for r in _parallel(ctx, jobs, workers=2)}
    out: dict[str, Any] = {}
    for prob, t_out in grids.items():
        for field in ("temp", "moist"):
            delivered = pd.read_parquet(ctx.dep(prob) / f"profiles_{field}.parquet")[COLS].to_numpy(float)
            fine = {k: np.asarray(results[f"{prob}:out:k{k}"][field], dtype=float) for k in levels}
            d1 = float(np.max(np.abs(delivered - fine[levels[0]])))
            d2 = float(np.max(np.abs(fine[levels[0]] - fine[levels[1]])))
            order = float(np.log2(d1 / d2)) if d2 > 0.0 else float("nan")
            used = order if 1.0 <= order <= 3.0 else 2.0
            limit = fine[levels[1]] + (fine[levels[1]] - fine[levels[0]]) / (2.0**used - 1.0)
            err = np.abs(delivered - limit)[1:]  # t = 0 is the exact initial condition
            times = t_out[1:]
            i, j = np.unravel_index(int(np.argmax(err)), err.shape)
            differing = np.round(delivered[1:], 4) != np.round(limit[1:], 4)
            bad_rows = np.where(differing.any(axis=1))[0]
            last_t = float(times[bad_rows[-1]]) if len(bad_rows) else 0.0
            after = err[times > last_t]
            # the digit-flip count is dominated by values that happen to sit on a rounding boundary; the
            # decisive statistic is how many cells carry a discretisation error above half an ulp of the
            # fourth decimal (5e-5), i.e. cells whose fourth decimal is genuinely not resolved
            coarse = err > 5e-5
            coarse_rows = np.where(coarse.any(axis=1))[0]
            last_coarse_t = float(times[coarse_rows[-1]]) if len(coarse_rows) else 0.0
            beyond = err[times > last_coarse_t]
            out[f"{prob}_{field}"] = {
                "cells_above_half_ulp": int(coarse.sum()),
                "last_above_half_ulp_t_s": last_coarse_t,
                "error_after_half_ulp": float(beyond.max()) if beyond.size else 0.0,
                "grid_levels": levels,
                "observed_order": order,
                "order_used": used,
                "production_error": float(err.max()),
                "max_at_t_s": float(times[i]),
                "max_at_r_cm": float(R_OUT_CM[j]),
                "cells_total": int(err.size),
                "cells_differing": int(differing.sum()),
                "last_differing_t_s": last_t,
                "error_after_last_differing": float(after.max()) if after.size else 0.0,
            }
            ctx.log.info("convergence.output_grid", quantity=f"{prob}_{field}", **out[f"{prob}_{field}"])
    return out


@stage(
    "convergence",
    deps=("ingest", "q1", "q2", "q3", "q4"),
    description="Grid/time convergence of all four problems (observed orders, drying-time convergence)",
)
def convergence(ctx: StageContext) -> dict[str, Any]:
    levels = int(ctx.cfg("a.grid.convergence_levels", 5))
    refine = int(ctx.cfg("a.grid.refine", 3))
    opts = time_options(ctx)
    horizon = float(ctx.cfg("a.drying.horizon_s", 432000))
    t3 = ctx.dep_data("q3", "summary.json")["t_dry_s"]
    t4 = ctx.dep_data("q4", "summary.json")["main"]["t_dry_s"]
    jobs: list[dict[str, Any]] = []
    for k in range(levels):
        n = 20 * 2**k
        for prob, t_out in (("q1", Q1_TABLE_T), ("q2", Q2_TABLE_T)):
            jobs.append(
                {
                    "label": f"{prob}:k{k}",
                    "recipe": spec_to_recipe(build_spec(ctx, prob)),
                    "n": n,
                    "opts": opts,
                    "mode": "fixed",
                    "t_out": np.concatenate([[0.0], t_out]).tolist(),
                }
            )
        for prob, t_dry in (("q3", t3), ("q4", t4)):
            jobs.append(
                {
                    "label": f"{prob}:k{k}",
                    "recipe": spec_to_recipe(build_spec(ctx, prob)),
                    "n": n,
                    "opts": opts,
                    "mode": "dry",
                    "horizon": horizon,
                    "snapshot_times": _table_times_hours(t_dry).tolist(),
                }
            )
    # time-tolerance, integrator and maximum-step cross-checks on the production grid (problem 3, the longest run)
    n_prod = 20 * 2**refine
    time_cases = (
        ("tol:strict", {"rtol": 1e-10, "atol": 1e-12}),
        ("tol:radau", {"method": "Radau"}),
        ("tol:loose", {"rtol": 1e-7, "atol": 1e-9}),
        ("dt:3600", {"max_step": 3600.0}),
        ("dt:900", {"max_step": 900.0}),
        ("dt:300", {"max_step": 300.0}),
    )
    for label, extra in time_cases:
        jobs.append(
            {
                "label": label,
                "recipe": spec_to_recipe(build_spec(ctx, "q3")),
                "n": n_prod,
                "opts": {**opts, **extra},
                "mode": "dry",
                "horizon": horizon,
                "snapshot_times": _table_times_hours(t3).tolist(),
            }
        )
    results = {r["label"]: r for r in _parallel(ctx, jobs)}
    report: dict[str, Any] = {"levels": levels, "production_refine": refine, "grid": {}, "time": {}}
    for prob in ("q1", "q2"):
        for field in ("temp", "moist"):
            lv = [
                {"k": k, "dr_cm": 0.1 / 2**k, "v": np.asarray(results[f"{prob}:k{k}"][field])[1:, ::5]}
                for k in range(levels)
            ]
            report["grid"][f"{prob}_{field}"] = verify.convergence_table(lv, "v")
    for prob in ("q3", "q4"):
        common = set.intersection(*(set(results[f"{prob}:k{k}"]["snapshot_t"]) for k in range(levels)))
        lv = []
        for k in range(levels):
            res = results[f"{prob}:k{k}"]
            keep = [i for i, t in enumerate(res["snapshot_t"]) if t in common]
            lv.append({"k": k, "dr_cm": 0.1 / 2**k, "v": np.asarray(res["snapshot_moist"])[keep][:, ::5]})
        report["grid"][f"{prob}_moist"] = verify.convergence_table(lv, "v")
        report["grid"][f"{prob}_moist_times_h"] = sorted(t / HOUR for t in common)
        dry = [results[f"{prob}:k{k}"]["t_dry"] for k in range(levels)]
        rows = [{"k": k, "dr_cm": 0.1 / 2**k, "t_dry_h": d / HOUR} for k, d in enumerate(dry)]
        orders_t = verify.triplet_orders([d / HOUR for d in dry])
        for i in range(levels - 1):
            rows[i]["diff_to_finest_min"] = (dry[i] - dry[-1]) / 60.0
        for i, row in enumerate(rows):
            row["observed_order_triplet"] = orders_t[i]
        report["grid"][f"{prob}_t_dry"] = rows
    base = results[f"q3:k{refine}"]
    for label, _extra in time_cases:
        r = results[label]
        common = set(r["snapshot_t"]) & set(base["snapshot_t"])
        a = np.asarray([row for t, row in zip(r["snapshot_t"], r["snapshot_moist"]) if t in common])
        b = np.asarray([row for t, row in zip(base["snapshot_t"], base["snapshot_moist"]) if t in common])
        report["time"][label] = {
            "t_dry_diff_s": float(r["t_dry"] - base["t_dry"]),
            "max_abs_diff_moist": float(np.max(np.abs(a - b))),
            "nfev": r["stats"]["nfev"],
            "options": r["opts"],
        }
    report["production"] = {
        "q1_temp_err": report["grid"]["q1_temp"][refine]["max_abs_error"] if refine < levels - 1 else None,
        "q1_moist_err": report["grid"]["q1_moist"][refine]["max_abs_error"] if refine < levels - 1 else None,
        "q2_temp_err": report["grid"]["q2_temp"][refine]["max_abs_error"] if refine < levels - 1 else None,
        "q2_moist_err": report["grid"]["q2_moist"][refine]["max_abs_error"] if refine < levels - 1 else None,
        "q3_moist_err": report["grid"]["q3_moist"][refine]["max_abs_error"] if refine < levels - 1 else None,
        "q4_moist_err": report["grid"]["q4_moist"][refine]["max_abs_error"] if refine < levels - 1 else None,
    }
    report["output_grid"] = _output_grid_study(ctx, refine, opts)
    ctx.write_json("convergence.json", report)
    ctx.write_json("raw_t_dry.json", {k: v.get("t_dry") for k, v in results.items() if "t_dry" in v})
    for key, row in report["output_grid"].items():
        tag = _tag(key)
        ctx.number(f"OutGridErr{tag}", row["production_error"], ".1e")
        ctx.number(f"OutGridOrder{tag}", row["observed_order"], ".2f")
        ctx.number(f"OutGridCells{tag}", row["cells_differing"])
        ctx.number(f"OutGridCellPct{tag}", 100.0 * row["cells_differing"] / row["cells_total"], ".2f")
        ctx.number(f"OutGridLastDiffSec{tag}", row["last_differing_t_s"], ".0f")
        ctx.number(f"OutGridErrAfter{tag}", row["error_after_last_differing"], ".1e")
        ctx.number(f"OutGridCoarse{tag}", row["cells_above_half_ulp"])
        ctx.number(f"OutGridCoarsePct{tag}", 100.0 * row["cells_above_half_ulp"] / row["cells_total"], ".2f")
        ctx.number(f"OutGridCoarseSec{tag}", row["last_above_half_ulp_t_s"], ".0f")
        ctx.number(f"OutGridErrCoarse{tag}", row["error_after_half_ulp"], ".1e")
    orders = {}
    for key, rows in report["grid"].items():
        if key.endswith(("_t_dry", "_times_h")):
            continue
        vals = [r["observed_order"] for r in rows if r.get("observed_order") is not None]
        orders[key] = vals[-1] if vals else None
        tag = _tag(key)
        if vals:
            ctx.number(f"ConvOrder{tag}", vals[-1], ".2f")
        triplet = [r["observed_order_triplet"] for r in rows if r.get("observed_order_triplet") is not None]
        if triplet:
            ctx.number(f"ConvTripletOrder{tag}", triplet[-1], ".2f")
        if refine < levels - 1:
            ctx.number(f"ConvErr{tag}", rows[refine]["max_abs_error"], ".1e")
    for prob in ("q3", "q4"):
        rows = report["grid"][f"{prob}_t_dry"]
        tag = "Qthree" if prob == "q3" else "Qfour"
        ctx.number(f"ConvDryDiff{tag}Min", abs(rows[refine]["diff_to_finest_min"]), ".2f")
        triplet = [r["observed_order_triplet"] for r in rows if r.get("observed_order_triplet") is not None]
        if triplet:
            ctx.number(f"ConvDryOrder{tag}", triplet[-1], ".2f")
    ctx.number("TimeStrictDiffMoist", report["time"]["tol:strict"]["max_abs_diff_moist"], ".1e")
    ctx.number("TimeRadauDiffMoist", report["time"]["tol:radau"]["max_abs_diff_moist"], ".1e")
    ctx.number("TimeStrictDryDiffSec", abs(report["time"]["tol:strict"]["t_dry_diff_s"]), ".2f")
    ctx.number("TimeRadauDryDiffSec", abs(report["time"]["tol:radau"]["t_dry_diff_s"]), ".2f")
    ctx.number("TimeStepCoarseExtraNfev", int(report["time"]["dt:3600"]["nfev"]) - int(base["stats"]["nfev"]))
    ctx.number("TimeStepCoarseDiffMoist", report["time"]["dt:3600"]["max_abs_diff_moist"], ".1e")
    ctx.number("TimeStepCoarseDryDiffSec", abs(report["time"]["dt:3600"]["t_dry_diff_s"]), ".2f")
    ctx.number("TimeStepFineDiffMoist", report["time"]["dt:300"]["max_abs_diff_moist"], ".1e")
    ctx.number("TimeStepFineDryDiffSec", abs(report["time"]["dt:300"]["t_dry_diff_s"]), ".2f")
    return {"orders": orders, "jobs": len(jobs)}


# ----------------------------------------------------------------------------------------------
# verification: analytic solution, independent 2-D model (cross-check and end effects)
# ----------------------------------------------------------------------------------------------


def _analytic_case(ctx: StageContext, field: str, levels: int, opts: dict[str, Any]) -> dict[str, Any]:
    """Constant-property Robin cylinder with constant ambient; solver vs Bessel series on several grids."""
    t_air, c_air = 50.0, 0.05
    d_const = float(APPENDIX2.diffusivity(np.array(C_INIT), np.array(T_INIT)))
    props = replace(APPENDIX2, d_const=d_const)
    chamber = Chamber.constant(t_air, c_air)
    spec = ProblemSpec(props, chamber, Radius.constant(R0))
    if field == "temp":
        t_out = np.concatenate([[0.0], Q1_TABLE_T])
        diffusivity, transfer, coeff, u0, u_inf = props.k0 / (props.rho0 * props.cp0), props.h, props.k0, T_INIT, t_air
    else:
        t_out = np.array([0.0, 300.0, 600.0, 1800.0, 3600.0, 7200.0, 14400.0, 28800.0, 43200.0, 86400.0])
        diffusivity, transfer, coeff, u0, u_inf = d_const, props.hm, d_const, C_INIT, c_air
    rows = []
    for k in range(levels):
        n = 20 * 2**k
        sol = solve(spec, n, t_out, **opts)
        num = getattr(sol, field)[1:, :: n // 20]
        exact = verify.analytic_field(np.linspace(0, 1, 21), t_out[1:], diffusivity, transfer, coeff, u0, u_inf)
        rows.append({"k": k, "dr_cm": 0.1 / 2**k, "max_abs_error": float(np.max(np.abs(num - exact)))})
    for i in range(len(rows) - 1):
        rows[i]["observed_order"] = float(np.log2(rows[i]["max_abs_error"] / rows[i + 1]["max_abs_error"]))
    return {
        "field": field,
        "biot": transfer * R0 / coeff,
        "diffusivity": diffusivity,
        "times_s": t_out[1:].tolist(),
        "levels": rows,
    }


def _duhamel_case(ctx: StageContext, levels: int, opts: dict[str, Any]) -> dict[str, Any]:
    """Problem 1's *actual* temperature field against the exact Duhamel series (time-dependent ambient).

    The ambient of problem 1 is not constant (28 -> 41.5 degC over 1800 s), so the separation-of-variables
    series is not its solution; Duhamel superposition over the piecewise-linear chamber history is. This gives
    a true-solution benchmark for the delivered temperature sheet of result1.xlsx, not only for a constructed
    constant-ambient case.
    """
    chamber = load_chamber(ctx)
    props = APPENDIX2
    alpha = props.k0 / (props.rho0 * props.cp0)
    spec = ProblemSpec(props, chamber, Radius.constant(R0))
    xi_table = R_TABLE_CM / (100.0 * R0)
    t_table = np.concatenate([[0.0], Q1_TABLE_T])
    exact_table = verify.duhamel_field(
        xi_table, t_table, alpha, props.h, props.k0, T_INIT, chamber.t, chamber.temp, chamber.temp_plateau
    )
    rows = []
    for k in range(levels):
        n = 20 * 2**k
        sol = solve(spec, n, t_table, **opts)
        num = sol.temp[1:, :: n // 4]
        rows.append({"k": k, "dr_cm": 0.1 / 2**k, "max_abs_error": float(np.max(np.abs(num - exact_table[1:])))})
    for i in range(len(rows) - 1):
        rows[i]["observed_order"] = float(np.log2(rows[i]["max_abs_error"] / rows[i + 1]["max_abs_error"]))
    delivered = pd.read_parquet(ctx.dep("q1") / "profiles_temp.parquet")
    t_full = delivered["t"].to_numpy(float)
    exact_full = verify.duhamel_field(
        R_OUT_CM / (100.0 * R0),
        t_full,
        alpha,
        props.h,
        props.k0,
        T_INIT,
        chamber.t,
        chamber.temp,
        chamber.temp_plateau,
    )
    err_full = np.abs(delivered[COLS].to_numpy(float) - exact_full)[1:]
    times = t_full[1:]
    i, j = np.unravel_index(int(np.argmax(err_full)), err_full.shape)
    return {
        "field": "temp",
        "biot": props.h * R0 / props.k0,
        "diffusivity": alpha,
        "ambient": "piecewise linear, attachment 1 (28.0 -> 41.5 degC)",
        "times_s": Q1_TABLE_T.tolist(),
        "levels": rows,
        "full_grid": {
            "cells": int(err_full.size),
            "max_abs_error": float(err_full.max()),
            "max_at_t_s": float(times[i]),
            "max_at_r_cm": float(R_OUT_CM[j]),
            "cells_differing_fourth_decimal": int(
                (np.round(delivered[COLS].to_numpy(float)[1:], 4) != np.round(exact_full[1:], 4)).sum()
            ),
        },
    }


@stage(
    "verification",
    deps=("ingest", "q1", "q2", "q3"),
    description="Analytic Bessel-series comparison and independent 2-D axisymmetric cross-check / end effects",
)
def verification(ctx: StageContext) -> dict[str, Any]:
    opts = time_options(ctx)
    refine = int(ctx.cfg("a.grid.refine", 3))
    levels = min(int(ctx.cfg("a.grid.convergence_levels", 5)), refine + 2)
    analytic = {f: _analytic_case(ctx, f, levels, opts) for f in ("temp", "moist")}
    analytic["duhamel"] = _duhamel_case(ctx, levels, opts)
    ctx.log.info("verification.analytic", **{f: a["levels"][refine]["max_abs_error"] for f, a in analytic.items()})
    ctx.log.info("verification.duhamel", **analytic["duhamel"]["full_grid"])
    ctx.write_json("analytic.json", analytic)

    # problems 2 and 3 are the same initial-boundary-value problem (appendix 3, fixed radius) integrated over
    # different horizons; result2.xlsx must therefore be the first 3 h of the whole-process solution (MDR-0005)
    q2m = pd.read_parquet(ctx.dep("q2") / "profiles_moist.parquet")
    q2t = pd.read_parquet(ctx.dep("q2") / "profiles_temp.parquet")
    q3m = pd.read_parquet(ctx.dep("q3") / "profiles_moist.parquet")
    q3t = pd.read_parquet(ctx.dep("q3") / "profiles_temp.parquet")
    shared = np.intersect1d(q2m["t"].to_numpy(float), q3m["t"].to_numpy(float))
    sel2, sel3 = np.isin(q2m["t"].to_numpy(float), shared), np.isin(q3m["t"].to_numpy(float), shared)
    horizon_check = {
        "shared_instants": len(shared),
        "max_abs_diff_moist": float(np.max(np.abs(q2m[COLS].to_numpy(float)[sel2] - q3m[COLS].to_numpy(float)[sel3]))),
        "max_abs_diff_temp": float(np.max(np.abs(q2t[COLS].to_numpy(float)[sel2] - q3t[COLS].to_numpy(float)[sel3]))),
    }
    horizon_check["passed"] = bool(
        horizon_check["max_abs_diff_moist"] < 1e-6 and horizon_check["max_abs_diff_temp"] < 1e-4
    )
    ctx.log.info("verification.horizon_consistency", **horizon_check)

    # independent 2-D axisymmetric model, coarse grid; insulated ends must reproduce the 1-D scheme exactly
    chamber = load_chamber(ctx)
    props = PROPS["q3"]
    nz = int(ctx.cfg("a.verification.axisym_nz", 25))
    k2 = int(ctx.cfg("a.verification.axisym_refine", 0))
    nr = 20 * 2**k2
    t3 = float(ctx.dep_data("q3", "summary.json")["t_dry_s"])
    t_out = np.arange(0.0, 12.0 * HOUR + 0.5, 600.0)
    spec = ProblemSpec(props, chamber, Radius.constant(R0))
    one_d = solve(spec, nr, t_out, **opts)
    ins = verify.solve_axisymmetric(props, chamber, t_out, nr=nr, nz=4, end_bc="insulated", **opts)
    cross = {
        "max_abs_diff_moist": float(np.max(np.abs(ins["moist"][:, :, 0] - one_d.moist))),
        "max_abs_diff_temp": float(np.max(np.abs(ins["temp"][:, :, 0] - one_d.temp))),
        "z_uniformity_moist": float(np.max(np.abs(ins["moist"] - ins["moist"][:, :, :1]))),
    }
    cross["passed"] = bool(cross["max_abs_diff_moist"] < 1e-6 and cross["max_abs_diff_temp"] < 1e-4)
    ctx.log.info("verification.cross_check", **cross)

    # end effects: Robin conditions on the end faces, whole drying period, coarse grid
    t_end = np.arange(0.0, math.ceil(1.05 * t3 / HOUR) * HOUR + 0.5, 1800.0)
    one_d_long = solve(spec, nr, t_end, detect_dry=True, **opts)
    two_d = verify.solve_axisymmetric(props, chamber, t_end, nr=nr, nz=nz, end_bc="robin", detect_dry=True, **opts)
    mid = two_d["moist"][:, :, 0]
    end_face = two_d["moist"][:, :, -1]
    diff_mid = np.abs(mid - one_d_long.moist)
    end_effect = {
        "nz": nz,
        "dz_cm": 12.5 / nz,
        "nr": nr,
        "t_dry_1d_h": None if one_d_long.t_dry is None else one_d_long.t_dry / HOUR,
        "t_dry_2d_h": two_d.get("t_dry", np.nan) / HOUR,
        "max_abs_diff_midplane_moist": float(diff_mid.max()),
        "max_abs_diff_midplane_temp": float(np.max(np.abs(two_d["temp"][:, :, 0] - one_d_long.temp))),
        "centre_moist_diff_at_1d_dry": float(
            abs(
                np.interp(one_d_long.t_dry or t3, t_end, mid[:, 0])
                - np.interp(one_d_long.t_dry or t3, t_end, one_d_long.moist[:, 0])
            )
        ),
        "end_face_centre_moist_at_1d_dry": float(np.interp(one_d_long.t_dry or t3, t_end, end_face[:, 0])),
        "axial_profile_at_1d_dry": {
            "z_cm": (two_d["z"] * 100).tolist(),
            "centre_line_moist": [
                float(np.interp(one_d_long.t_dry or t3, t_end, two_d["moist"][:, 0, j])) for j in range(nz + 1)
            ],
        },
    }
    if one_d_long.t_dry is not None and "t_dry" in two_d:
        end_effect["t_dry_rel_diff_pct"] = 100.0 * (two_d["t_dry"] / one_d_long.t_dry - 1.0)
    ctx.log.info("verification.end_effect", **{k: v for k, v in end_effect.items() if not isinstance(v, dict)})
    np.savez_compressed(
        ctx.out("axisym_snapshots.npz"),
        t=t_end,
        r=two_d["r"],
        z=two_d["z"],
        moist_mid=mid,
        moist_end=end_face,
        moist_1d=one_d_long.moist,
        centre_line=two_d["moist"][:, 0, :],
    )
    report = {
        "analytic": analytic,
        "cross_check_2d": cross,
        "end_effect_2d": end_effect,
        "horizon_consistency": horizon_check,
    }
    report["passed"] = bool(
        cross["passed"]
        and horizon_check["passed"]
        and all(a["levels"][refine]["max_abs_error"] < 1e-4 for a in analytic.values())
    )
    ctx.write_json("verification_report.json", report)
    ctx.number("AnalyticErrTemp", analytic["temp"]["levels"][refine]["max_abs_error"], ".1e")
    ctx.number("AnalyticErrMoist", analytic["moist"]["levels"][refine]["max_abs_error"], ".1e")
    ctx.number("AnalyticOrderTemp", analytic["temp"]["levels"][refine - 1]["observed_order"], ".2f")
    ctx.number("AnalyticOrderMoist", analytic["moist"]["levels"][refine - 1]["observed_order"], ".2f")
    ctx.number("AnalyticBiotHeat", analytic["temp"]["biot"], ".3f")
    ctx.number("AnalyticBiotMass", analytic["moist"]["biot"], ".3f")
    duh = analytic["duhamel"]
    ctx.number("DuhamelErrTable", duh["levels"][refine]["max_abs_error"], ".1e")
    ctx.number("DuhamelOrder", duh["levels"][refine - 1]["observed_order"], ".2f")
    ctx.number("DuhamelErrFull", duh["full_grid"]["max_abs_error"], ".1e")
    ctx.number("DuhamelCellsDiffer", duh["full_grid"]["cells_differing_fourth_decimal"])
    ctx.number("DuhamelCells", duh["full_grid"]["cells"])
    ctx.number("HorizonDiffMoist", horizon_check["max_abs_diff_moist"], ".1e")
    ctx.number("HorizonDiffTemp", horizon_check["max_abs_diff_temp"], ".1e")
    ctx.number("HorizonSharedInstants", horizon_check["shared_instants"])
    ctx.number("CrossCheckDiffMoist", cross["max_abs_diff_moist"], ".1e")
    ctx.number("CrossCheckDiffTemp", cross["max_abs_diff_temp"], ".1e")
    ctx.number("EndEffectMidplaneDiff", end_effect["max_abs_diff_midplane_moist"], ".1e")
    ctx.number("EndEffectCentreDiffAtDry", end_effect["centre_moist_diff_at_1d_dry"], ".1e")
    if "t_dry_rel_diff_pct" in end_effect:
        ctx.number("EndEffectDryDiffPct", end_effect["t_dry_rel_diff_pct"], ".3f")
        ctx.number("EndEffectDryHoursTwoD", end_effect["t_dry_2d_h"], ".2f")
        ctx.number("EndEffectDryHoursOneD", end_effect["t_dry_1d_h"], ".2f")
    ctx.number("EndEffectDzCm", end_effect["dz_cm"], ".2f")
    return {
        "passed": report["passed"],
        "cross_check": cross["passed"],
        "duhamel_full_grid_error": duh["full_grid"]["max_abs_error"],
        "horizon_consistency": horizon_check["passed"],
        "end_effect_pct": end_effect.get("t_dry_rel_diff_pct"),
    }


# ----------------------------------------------------------------------------------------------
# sensitivity: parameter elasticities of the drying time and alternative interpretations
# ----------------------------------------------------------------------------------------------


@stage(
    "sensitivity",
    deps=("ingest", "q3", "q4"),
    description="Drying-time elasticities (h, h_m, D0, C_air, T plateau, shrinkage) and alternative interpretations",
)
def sensitivity(ctx: StageContext) -> dict[str, Any]:
    factors = [float(f) for f in ctx.cfg("a.sensitivity.factors", [0.8, 0.9, 1.1, 1.2])]
    n = production_n(ctx)
    opts = time_options(ctx)
    horizon = float(ctx.cfg("a.drying.horizon_s", 432000))
    base3 = build_spec(ctx, "q3")
    base4 = build_spec(ctx, "q4")
    t3 = float(ctx.dep_data("q3", "summary.json")["t_dry_s"])
    t4 = float(ctx.dep_data("q4", "summary.json")["main"]["t_dry_s"])

    def job(label: str, spec: ProblemSpec) -> dict[str, Any]:
        return {"label": label, "recipe": spec_to_recipe(spec), "n": n, "opts": opts, "mode": "dry", "horizon": horizon}

    jobs: list[dict[str, Any]] = []
    params = {
        "h": lambda f: replace(base3, props=base3.props.scaled(h=f)),
        "hm": lambda f: replace(base3, props=base3.props.scaled(hm=f)),
        "D0": lambda f: replace(base3, props=base3.props.scaled(d0=f)),
        "Cair": lambda f: replace(base3, chamber=replace(base3.chamber, hum_scale=f)),
        "Tplateau": lambda f: replace(base3, chamber=replace(base3.chamber, temp_scale=f)),
        "shrink": lambda f: replace(base4, radius=load_radius(ctx, f)),
        "D0q4": lambda f: replace(base4, props=base4.props.scaled(d0=f)),
    }
    for name, make in params.items():
        for f in factors:
            jobs.append(job(f"{name}:{f}", make(f)))
    alternatives = {
        "plateau_last_sample": replace(base3, chamber=load_chamber(ctx, "last_sample")),
        "plateau_nominal": replace(base3, chamber=load_chamber(ctx, "nominal")),
        "face_arithmetic": replace(base3, face_scheme="arithmetic"),
        "face_harmonic": replace(base3, face_scheme="harmonic"),
        "mass_conservative": replace(base3, mass_form="conservative"),
    }
    for name, spec in alternatives.items():
        jobs.append(job(f"alt:{name}", spec))
    results = {r["label"]: r for r in _parallel(ctx, jobs)}

    table: list[dict[str, Any]] = []
    elasticities: dict[str, float] = {}
    for name in params:
        base_t = t4 if name in ("shrink", "D0q4") else t3
        row: dict[str, Any] = {"parameter": name, "base_h": base_t / HOUR}
        for f in factors:
            row[f"x{f}"] = results[f"{name}:{f}"]["t_dry"] / HOUR
        if 0.9 in factors and 1.1 in factors:
            e = ((results[f"{name}:1.1"]["t_dry"] - results[f"{name}:0.9"]["t_dry"]) / base_t) / 0.2
            row["elasticity"] = e
            elasticities[name] = e
        table.append(row)
    alt_rows = [
        {
            "case": name,
            "t_dry_h": results[f"alt:{name}"]["t_dry"] / HOUR,
            "diff_pct": 100.0 * (results[f"alt:{name}"]["t_dry"] / t3 - 1.0),
        }
        for name in alternatives
    ]
    q4s = ctx.dep_data("q4", "summary.json")
    q3s = ctx.dep_data("q3", "summary.json")
    alt_rows.append(
        {
            "case": "q4_" + q4s["alternative"]["formulation"],
            "t_dry_h": q4s["alternative"]["t_dry_h"],
            "diff_pct": 100.0 * (q4s["alternative"]["t_dry_h"] / q4s["main"]["t_dry_h"] - 1.0),
        }
    )
    alt_rows.append(
        {
            "case": "q4_fixed_radius",
            "t_dry_h": q4s["fixed_radius"]["t_dry_h"],
            "diff_pct": 100.0 * (q4s["fixed_radius"]["t_dry_h"] / q4s["main"]["t_dry_h"] - 1.0),
        }
    )
    if q3s.get("t_mean_below_h") is not None:
        alt_rows.append(
            {
                "case": "q3_mean_criterion",
                "t_dry_h": q3s["t_mean_below_h"],
                "diff_pct": 100.0 * (q3s["t_mean_below_h"] / q3s["t_dry_h"] - 1.0),
            }
        )
    ctx.write_json(
        "sensitivity.json", {"factors": factors, "table": table, "elasticities": elasticities, "alternatives": alt_rows}
    )
    for name, e in elasticities.items():
        ctx.number("Elas" + _tag(name), e, ".3f")
    for row in alt_rows:
        key = "Alt" + _tag(row["case"])
        ctx.number(key + "Hours", row["t_dry_h"], ".2f")
        ctx.number(key + "Pct", row["diff_pct"], ".2f")
    return {"jobs": len(jobs), "elasticities": elasticities}
