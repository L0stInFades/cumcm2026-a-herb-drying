"""Tests of the model-evaluation extension (MDR-0008): psychrometrics, the energy flux cap, the C_eq driving force."""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest

from pipelines.a import verify
from pipelines.a.common import spec_from_recipe, spec_to_recipe
from pipelines.a.extension import latent_diagnostic, slug
from pipelines.a.physics import (
    APPENDIX2,
    APPENDIX3,
    C_INIT,
    Chamber,
    FluxCap,
    Radius,
    dry_solid_density,
    humidity_ratio,
    latent_heat,
    relative_humidity,
    saturation_pressure,
    wet_bulb_temperature,
)
from pipelines.a.solver import ProblemSpec, solve, solve_until_dry


def test_psychrometrics_reproduce_reference_values() -> None:
    assert saturation_pressure(50.0) == pytest.approx(12352.0, rel=0.01)  # steam tables: 12.35 kPa
    assert saturation_pressure(20.0) == pytest.approx(2339.0, rel=0.01)
    assert relative_humidity(50.0, 0.05) == pytest.approx(0.61, abs=0.02)
    tw = wet_bulb_temperature(50.0, 0.05)
    assert 40.5 < tw < 43.0
    w_sat = float(humidity_ratio(saturation_pressure(30.0)))
    assert wet_bulb_temperature(30.0, w_sat) == pytest.approx(30.0)  # saturated air: wet bulb = dry bulb
    assert wet_bulb_temperature(50.0, 0.02) < wet_bulb_temperature(50.0, 0.05) < 50.0
    assert latent_heat(50.0) == pytest.approx(2.383e6, rel=0.005)  # steam tables: 2382 kJ/kg
    assert dry_solid_density(APPENDIX3, C_INIT) == pytest.approx((650 + 128 * 2.55) / 3.55)


def test_flux_cap_limits_early_loss_and_keeps_the_balance_exact() -> None:
    chamber = Chamber.constant(50.0, 0.05)
    rho_s = float(dry_solid_density(APPENDIX3, C_INIT))
    cap = FluxCap.wet_bulb(chamber, APPENDIX3.h, rho_s)
    assert 5e-5 < cap.plateau < 1.5e-4  # ~ 25 W/(m^2 K) x 8 K / 2.4 MJ/kg
    assert cap.value(1e9) == pytest.approx(cap.plateau)
    base = ProblemSpec(APPENDIX3, chamber, Radius.constant())
    capped = replace(base, flux_cap=cap)
    t_out = np.arange(0.0, 3601.0, 60.0)
    sol_b = solve(base, 20, t_out)
    sol_c = solve(capped, 20, t_out)
    vol = verify.control_volumes(sol_b.xi)
    assert sol_c.moist[-1] @ vol > sol_b.moist[-1] @ vol  # less water leaves when the flux is capped
    flux = capped.surface_moisture_loss(sol_c.t, sol_c.moist[:, -1], sol_c.air_hum)
    assert np.all(flux <= cap.plateau / rho_s * (1 + 1e-12))
    mb = verify.moisture_balance(
        sol_c.t,
        sol_c.xi,
        sol_c.moist,
        sol_c.radius,
        sol_c.air_hum,
        APPENDIX3.hm,
        cum_loss=sol_c.cum_moist_loss,
        flux=flux,
    )
    assert mb["exact"]["max_relative_residual"] < 1e-8
    # below the release moisture the loss law is the problem's own
    c_release = 0.05 + cap.plateau / (rho_s * APPENDIX3.hm)
    low = np.array([0.5 * c_release])
    assert capped.surface_moisture_loss(0.0, low, 0.05) == pytest.approx(APPENDIX3.hm * (low - 0.05))
    # the cap survives the worker recipe round trip
    again = spec_from_recipe(spec_to_recipe(capped))
    assert again.flux_cap is not None and again.flux_cap.value(30.0) == pytest.approx(cap.value(30.0))
    assert slug("cap:1.0") == "cap_1p0"


def test_equilibrium_moisture_driving_force_bounds_the_surface_and_is_monotone_for_constant_d() -> None:
    # constant D: the problem is linear, max C < 0.15 <=> theta < (0.15 - C_eq)/(2.55 - C_eq), so a larger C_eq
    # must take longer; with the exponential D(C) of the appendices a drier surface forms a low-D crust and the
    # ordering is not guaranteed (case hardening), which is exactly what the extension stage quantifies.
    props = replace(APPENDIX2, d_const=2e-8)
    base = ProblemSpec(props, Chamber.constant(50.0, 0.05), Radius.constant())
    slow = replace(base, chamber=Chamber.constant(50.0, 0.10))
    sol_base = solve_until_dry(base, 20, 60.0, horizon=3600.0)
    sol_slow = solve_until_dry(slow, 20, 60.0, horizon=3600.0)
    assert sol_base.t_dry is not None and sol_slow.t_dry is not None
    assert sol_slow.t_dry > sol_base.t_dry
    assert sol_slow.moist[:, -1].min() >= 0.10 - 1e-9  # the surface never drops below the driving-force floor
    assert sol_slow.dry_profile is not None and sol_slow.dry_profile.max() == pytest.approx(0.15, abs=1e-6)


def test_latent_diagnostic_scales_with_the_evaporation_flux() -> None:
    spec = ProblemSpec(APPENDIX3, Chamber.constant(50.0, 0.05), Radius.constant())
    sol = solve(spec, 20, np.arange(0.0, 3601.0, 60.0))
    rho_s = float(dry_solid_density(APPENDIX3, C_INIT))
    hist = {
        "t": sol.t,
        "surface_temp": sol.temp[:, -1],
        "air_temp": sol.air_temp,
        "surface_moist": sol.moist[:, -1],
        "air_hum": sol.air_hum,
        "mean_moist": 2.0 * (sol.moist @ verify.control_volumes(sol.xi)),
    }
    diag = latent_diagnostic(hist, spec, rho_s)
    j0 = rho_s * APPENDIX3.hm * (C_INIT - 0.05)
    assert diag["summary"]["j_w_initial"] == pytest.approx(j0)
    assert diag["summary"]["dt_lat_peak_K"] == pytest.approx(float(latent_heat(28.0)) * j0 / APPENDIX3.h)
    assert diag["arrays"]["e_lat"][-1] > 0 and np.all(np.diff(diag["arrays"]["e_lat"]) >= 0)
    assert diag["summary"]["water_removed_kg_per_m"] > 0
