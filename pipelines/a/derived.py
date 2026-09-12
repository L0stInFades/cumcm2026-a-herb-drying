"""Derived quantities quoted in the paper: dimensionless groups, time scales, extrapolated errors, solver statistics.

Every number the manuscript quotes must come from a stage's ``ctx.number`` registration (docs/ENGINEERING_STANDARD.md,
MDR-0009). This stage computes the a-priori estimates used in the analysis and verification sections from the stored
outputs of the science stages only (no new simulations), so that the text of the paper carries no hand-typed arithmetic.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from forge.context import StageContext
from forge.runner import stage
from pipelines.a.physics import APPENDIX3, C_INIT, C_TARGET, LENGTH, R0, thermal_diffusivity

HOUR = 3600.0
Q1_END_S = 1800.0


def fourier_number(diffusivity: float, t: float, length: float) -> float:
    """Fo = D t / L^2."""
    return float(diffusivity * t / length**2)


def richardson_error(values: list[float], order: float = 2.0) -> tuple[float, float]:
    """Richardson extrapolation of a sequence computed on grids halved at every level (finest last).

    Returns ``(limit, error)``: the limit estimated from the two finest levels and the absolute error of the
    second-finest level (the production grid) with respect to that limit.
    """
    if len(values) < 2:
        raise ValueError("need at least two levels")
    fine, coarse = float(values[-1]), float(values[-2])
    limit = fine + (fine - coarse) / (2.0**order - 1.0)
    return limit, abs(coarse - limit)


def first_time_within(t: np.ndarray, a: np.ndarray, b: np.ndarray, tol: float) -> float | None:
    """First instant after which |a - b| stays below ``tol`` for the rest of the record (None if it never does)."""
    ok = np.abs(np.asarray(a, dtype=float) - np.asarray(b, dtype=float)) < tol
    if not ok.any():
        return None
    bad = np.where(~ok)[0]
    start = 0 if len(bad) == 0 else int(bad[-1]) + 1
    return None if start >= len(t) else float(np.asarray(t, dtype=float)[start])


@stage(
    "derived",
    deps=("q1", "q2", "q3", "q4", "convergence"),
    description="Dimensionless groups, time scales, extrapolated grid errors and solver statistics quoted in the paper",
)
def derived(ctx: StageContext) -> dict[str, Any]:
    q1 = ctx.dep_data("q1", "summary.json")
    q2h = ctx.dep_data("q2", "history.json")
    q3 = ctx.dep_data("q3", "summary.json")
    q4 = ctx.dep_data("q4", "summary.json")
    conv = ctx.dep_data("convergence", "convergence.json")

    # problem 1 (appendix 2, initial state): time scales, Fourier numbers, diffusion length, stiffness ratio
    alpha1 = float(q1["biot"]["alpha"])
    d1 = float(q1["biot"]["D"])
    t_heat = R0**2 / alpha1
    t_mass = R0**2 / d1
    fo_heat = fourier_number(alpha1, Q1_END_S, R0)
    fo_mass = fourier_number(d1, Q1_END_S, R0)
    diffusion_length = float(np.sqrt(d1 * Q1_END_S))

    # problem 2: thermal equilibration of the centre with the chamber
    t2 = np.asarray(q2h["t"], dtype=float)
    gap_1k = first_time_within(t2, q2h["centre_temp"], q2h["air_temp"], 1.0)
    gap_half = first_time_within(t2, q2h["centre_temp"], q2h["air_temp"], 0.5)

    # problem 3 (appendix 3): axial Fourier numbers for the 1-D reduction, diffusivity range (stiffness)
    t3 = float(q3["t_dry_s"])
    t_plateau = float(q3["spec"]["chamber"]["temp_plateau"])
    alpha3 = thermal_diffusivity(APPENDIX3, C_INIT)
    d_max = float(APPENDIX3.diffusivity(np.array(C_INIT), np.array(t_plateau)))
    d_target = float(APPENDIX3.diffusivity(np.array(C_TARGET), np.array(t_plateau)))
    d_surface_dry = float(APPENDIX3.diffusivity(np.array(q3["surface_moist_at_dry"]), np.array(t_plateau)))
    half_length = LENGTH / 2.0
    fo_axial_heat = fourier_number(alpha3, t3, half_length)
    fo_axial_mass = fourier_number(d_max, t3, half_length)

    # drying-time discretisation error of the production grid (Richardson, order 2, two finest levels)
    rows3 = conv["grid"]["q3_t_dry"]
    rows4 = conv["grid"]["q4_t_dry"]
    limit3, err3 = richardson_error([float(r["t_dry_h"]) for r in rows3])
    limit4, err4 = richardson_error([float(r["t_dry_h"]) for r in rows4])

    # problem 4: shrinkage geometry
    r_end_cm = float(q4["spec"]["radius"]["r_end_cm"])
    r0_cm = float(q4["spec"]["radius"]["r_start_cm"])
    area_ratio = (r_end_cm / r0_cm) ** 2

    n = int(q1["n"])
    report = {
        "aspect_ratio": LENGTH / R0,
        "system_size": 2 * (n + 1) + 2,
        "q1": {
            "alpha": alpha1,
            "D0": d1,
            "t_heat_s": t_heat,
            "t_mass_h": t_mass / HOUR,
            "stiffness_ratio": alpha1 / d1,
            "Fo_heat_1800s": fo_heat,
            "Fo_mass_1800s": fo_mass,
            "diffusion_length_cm": diffusion_length * 100.0,
        },
        "q2": {"centre_within_1K_h": None if gap_1k is None else gap_1k / HOUR, "centre_within_0p5K_h": gap_half},
        "q3": {
            "alpha3_C0": alpha3,
            "D_max": d_max,
            "D_at_target": d_target,
            "D_surface_at_dry": d_surface_dry,
            "D_ratio_max_to_surface": d_max / d_surface_dry,
            "Fo_axial_heat": fo_axial_heat,
            "Fo_axial_mass_upper": fo_axial_mass,
            "richardson_limit_h": limit3,
            "production_error_min": err3 * 60.0,
            "stats": q3["stats"],
        },
        "q4": {
            "richardson_limit_h": limit4,
            "production_error_min": err4 * 60.0,
            "r_end_cm": r_end_cm,
            "area_ratio": area_ratio,
            "inverse_radius_ratio_squared": 1.0 / area_ratio,
            "stats": q4["stats"],
        },
    }
    ctx.write_json("derived.json", report)

    ctx.number("AspectRatio", LENGTH / R0, ".1f")
    ctx.number("SystemSize", 2 * (n + 1) + 2)
    ctx.number("QoneHeatTimeScaleS", t_heat, ".0f")
    ctx.number("QoneMassTimeScaleH", t_mass / HOUR, ".1f")
    ctx.number("QoneStiffnessRatio", alpha1 / d1, ".0f")
    ctx.number("QoneFourierHeat", fo_heat, ".2f")
    ctx.number("QoneFourierMass", fo_mass, ".3f")
    ctx.number("QoneDiffusionLengthCm", diffusion_length * 100.0, ".2f")
    if gap_1k is not None:
        ctx.number("QtwoCentreWithinOneKelvinHours", gap_1k / HOUR, ".2f")
    if gap_half is not None:
        ctx.number("QtwoCentreWithinHalfKelvinHours", gap_half / HOUR, ".2f")
    ctx.number("FoAxialHeat", fo_axial_heat, ".2f")
    ctx.number("FoAxialMassUpper", fo_axial_mass, ".2f")
    ctx.number("DmaxAppendixThree", d_max, ".2e")
    ctx.number("DtargetAppendixThree", d_target, ".1e")
    ctx.number("DsurfaceAtDryQthree", d_surface_dry, ".1e")
    ctx.number("DratioQthree", d_max / d_surface_dry, ".0f")
    ctx.number("QthreeDryHoursExtrapolated", limit3, ".4f")
    ctx.number("QthreeGridErrorMin", err3 * 60.0, ".2f")
    ctx.number("QfourDryHoursExtrapolated", limit4, ".4f")
    ctx.number("QfourGridErrorMin", err4 * 60.0, ".2f")
    ctx.number("QthreeNfev", int(q3["stats"]["nfev"]))
    ctx.number("QthreeNlu", int(q3["stats"]["nlu"]))
    ctx.number("QthreeSegments", int(q3["stats"]["segments"]))
    ctx.number("QfourNfev", int(q4["stats"]["nfev"]))
    ctx.number("QfourRadiusFinalCm", r_end_cm, ".3f")
    ctx.number("QfourAreaRatioPct", 100.0 * area_ratio, ".1f")
    ctx.number("QfourRadiusRatioSquared", 1.0 / area_ratio, ".2f")
    ctx.number("ChamberSamples", int(q3["spec"]["chamber"]["samples"]))
    ctx.number("ChamberEndHours", float(q3["spec"]["chamber"]["t_end"]) / HOUR, ".0f")
    ctx.number("RadiusSamples", int(q4["spec"]["radius"]["samples"]))
    ctx.number("RadiusEndHours", float(q4["spec"]["radius"]["t_end"]) / HOUR, ".0f")
    return {"numbers": 33, "fo_axial_heat": fo_axial_heat, "fo_axial_mass_upper": fo_axial_mass}
