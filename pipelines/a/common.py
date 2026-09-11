"""Shared helpers for the Problem-A science stages: loading inputs, building specs, storing solutions."""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd

from forge.context import StageContext
from pipelines.a.physics import (
    APPENDIX2,
    APPENDIX3,
    APPENDIX4,
    R0,
    Chamber,
    Properties,
    Radius,
)
from pipelines.a.solver import ProblemSpec, Solution, sample_physical, subgrid_indices

R_OUT_CM = np.round(np.arange(0, 21) * 0.1, 1)  # 0.0 … 2.0 cm
R_TABLE_CM = np.array([0.0, 0.5, 1.0, 1.5, 2.0])
COLS = [f"{r:.1f}" for r in R_OUT_CM]
PROPS = {"q1": APPENDIX2, "q2": APPENDIX3, "q3": APPENDIX3, "q4": APPENDIX4}


def load_chamber(ctx: StageContext, plateau: str | None = None) -> Chamber:
    df = pd.read_parquet(ctx.dep("ingest") / "data" / "附件1__Sheet1.parquet")
    rule = plateau or str(ctx.cfg("a.chamber.plateau", "mean_last_hour"))
    return Chamber.from_series(
        df["时间"].to_numpy(float), df["温度"].to_numpy(float), df["水分浓度"].to_numpy(float), rule
    )


def load_radius(ctx: StageContext, shrink_scale: float = 1.0) -> Radius:
    df = pd.read_parquet(ctx.dep("ingest") / "data" / "附件2__Sheet1.parquet")
    return Radius.from_series(df["时间"].to_numpy(float), df["半径"].to_numpy(float), shrink_scale)


def production_n(ctx: StageContext) -> int:
    return 20 * 2 ** int(ctx.cfg("a.grid.refine", 3))


def time_options(ctx: StageContext) -> dict[str, Any]:
    return {
        "method": str(ctx.cfg("a.time.method", "BDF")),
        "rtol": float(ctx.cfg("a.time.rtol", 1e-7)),
        "atol": float(ctx.cfg("a.time.atol", 1e-9)),
    }


def build_spec(
    ctx: StageContext,
    problem: str,
    *,
    props: Properties | None = None,
    chamber: Chamber | None = None,
    radius: Radius | None = None,
    formulation: str | None = None,
    face_scheme: str | None = None,
) -> ProblemSpec:
    """Default specification of problem q1…q4 with optional overrides (sensitivity / alternatives)."""
    if radius is None:
        radius = load_radius(ctx) if problem == "q4" else Radius.constant(R0)
    return ProblemSpec(
        props=props or PROPS[problem],
        chamber=chamber or load_chamber(ctx),
        radius=radius,
        formulation=formulation or str(ctx.cfg("a.scheme.formulation", "lagrangian")),
        face_scheme=face_scheme or str(ctx.cfg("a.scheme.face", "midpoint")),
    )


def profiles_frame(sol: Solution, field: str) -> pd.DataFrame:
    """Field values on the 0.1 cm output grid (fixed radius: coincident nodes) as a DataFrame with a 't' column."""
    values = getattr(sol, field)
    idx = subgrid_indices(len(sol.xi) - 1)
    frame = pd.DataFrame(values[:, idx], columns=COLS)
    frame.insert(0, "t", sol.t)
    return frame


def physical_frame(sol: Solution, field: str) -> pd.DataFrame:
    """Field values at fixed physical radii 0…2 cm (NaN beyond R(t)) plus the surface value, shrinking domain."""
    values = getattr(sol, field)
    rows = np.full((len(sol.t), len(R_OUT_CM)), np.nan)
    for i in range(len(sol.t)):
        rows[i] = sample_physical(values[i], sol.xi, float(sol.radius[i]), R_OUT_CM / 100.0)
    frame = pd.DataFrame(rows, columns=COLS)
    frame.insert(0, "t", sol.t)
    frame["surface"] = values[:, -1]
    frame["radius_cm"] = sol.radius * 100.0
    return frame


def table_cells(frame: pd.DataFrame, times: np.ndarray, radii_cm: np.ndarray = R_TABLE_CM) -> list[list[float]]:
    """Rows [t, v(r_1), …] of a profiles frame at the requested instants (must exist exactly in 't')."""
    out: list[list[float]] = []
    lookup = {float(t): i for i, t in enumerate(frame["t"].to_numpy(float))}
    for t in times:
        i = lookup[float(t)]
        out.append([float(t)] + [float(frame.at[i, f"{r:.1f}"]) for r in radii_cm])
    return out


def snapshot_arrays(sol: Solution, times: np.ndarray) -> dict[str, np.ndarray]:
    lookup = {float(t): i for i, t in enumerate(sol.t)}
    idx = [lookup[float(t)] for t in times if float(t) in lookup]
    return {
        "t": sol.t[idx],
        "xi": sol.xi,
        "temp": sol.temp[idx],
        "moist": sol.moist[idx],
        "radius": sol.radius[idx],
    }


def mean_moisture(sol: Solution) -> np.ndarray:
    """Volume-averaged (dry-matter-weighted for the Lagrangian case) moisture history: 2 * sum vol_i C_i."""
    from pipelines.a.verify import control_volumes

    return 2.0 * (sol.moist @ control_volumes(sol.xi))


def spec_to_recipe(spec: ProblemSpec) -> dict[str, Any]:
    """Picklable description of a spec (workers rebuild the interpolants themselves)."""
    from dataclasses import asdict

    ch = spec.chamber
    rd = spec.radius
    return {
        "props": asdict(spec.props),
        "chamber": {
            "t": ch.t.tolist(),
            "temp": ch.temp.tolist(),
            "hum": ch.hum.tolist(),
            "temp_plateau": ch.temp_plateau,
            "hum_plateau": ch.hum_plateau,
            "temp_scale": ch.temp_scale,
            "hum_scale": ch.hum_scale,
        },
        "radius": {"t": rd.t.tolist(), "r": rd.r.tolist(), "shrink_scale": rd.shrink_scale},
        "formulation": spec.formulation,
        "face_scheme": spec.face_scheme,
        "t_init": spec.t_init,
        "c_init": spec.c_init,
    }


def spec_from_recipe(recipe: dict[str, Any]) -> ProblemSpec:
    ch = recipe["chamber"]
    rd = recipe["radius"]
    return ProblemSpec(
        props=Properties(**recipe["props"]),
        chamber=Chamber(
            np.asarray(ch["t"], float),
            np.asarray(ch["temp"], float),
            np.asarray(ch["hum"], float),
            float(ch["temp_plateau"]),
            float(ch["hum_plateau"]),
            float(ch["temp_scale"]),
            float(ch["hum_scale"]),
        ),
        radius=Radius(np.asarray(rd["t"], float), np.asarray(rd["r"], float), float(rd["shrink_scale"])),
        formulation=str(recipe["formulation"]),
        face_scheme=str(recipe["face_scheme"]),
        t_init=float(recipe["t_init"]),
        c_init=float(recipe["c_init"]),
    )
