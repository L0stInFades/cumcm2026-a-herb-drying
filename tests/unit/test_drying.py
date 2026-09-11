"""Tests of the Problem-A physics, solver and independent verifiers (small grids, short horizons)."""

from __future__ import annotations

import math
from dataclasses import replace

import numpy as np
import pytest
from hypothesis import given
from hypothesis import strategies as st
from scipy.special import j0, j1

from pipelines.a import verify
from pipelines.a.physics import (
    APPENDIX2,
    APPENDIX3,
    APPENDIX4,
    C_INIT,
    R0,
    T_INIT,
    Chamber,
    Radius,
)
from pipelines.a.solver import (
    ProblemSpec,
    RadialSystem,
    make_grid,
    sample_physical,
    solve,
    solve_until_dry,
    subgrid_indices,
)

# ---- physics --------------------------------------------------------------------------------------


def test_property_laws_match_the_appendices() -> None:
    c = np.array(2.55)
    assert APPENDIX2.diffusivity(c, np.array(28.0)) == pytest.approx(7e-9 * math.exp(-0.89 / 2.55))
    assert APPENDIX3.rho(c) == pytest.approx(650 + 128 * 2.55)
    assert APPENDIX3.cp(c) == pytest.approx(1450 + 2736 * 2.55 / 3.55)
    assert APPENDIX3.k(c) == pytest.approx(0.21 + 0.38 * 2.55 / 3.55)
    d3 = 2.4e-3 * math.exp(-0.45 / 2.55) * math.exp(-3850 / (28 + 273.15))
    assert APPENDIX3.diffusivity(c, np.array(28.0)) == pytest.approx(d3, rel=1e-12)
    d4 = 4.2e-4 * math.exp(-0.30 / 0.2) * math.exp(-3850 / (50 + 273.15))
    assert APPENDIX4.diffusivity(np.array(0.2), np.array(50.0)) == pytest.approx(d4, rel=1e-12)
    assert APPENDIX4.h == 25.0 and APPENDIX4.hm == 8e-7
    assert APPENDIX2.diffusivity(np.array(0.0), np.array(28.0)) == 0.0  # floor: no division by zero


def test_chamber_interpolates_linearly_and_holds_the_plateau() -> None:
    t = np.arange(0, 14401, 60, dtype=float)
    temp = 28 + 22 * np.minimum(t / 7200, 1.0)
    hum = np.full_like(t, 0.05)
    ch = Chamber.from_series(t, temp, hum, "mean_last_hour")
    assert ch.air_temperature(30.0) == pytest.approx(28 + 22 * 30 / 7200)
    assert ch.air_temperature(1e6) == pytest.approx(50.0)
    assert Chamber.from_series(t, temp, hum, "nominal").hum_plateau == 0.05
    assert Chamber.from_series(t, temp, hum, "last_sample").temp_plateau == pytest.approx(temp[-1])


def test_radius_is_monotone_clamped_and_scaled() -> None:
    t = np.arange(0, 7201, 1800, dtype=float)
    r_cm = np.array([2.0, 1.8, 1.7, 1.7, 1.65])
    rad = Radius.from_series(t, r_cm)
    fine = np.linspace(0, 9000, 500)
    values = rad.value(fine)
    assert np.all(np.diff(values) <= 1e-15)
    assert rad.value(0.0) == pytest.approx(0.02) and rad.value(1e5) == pytest.approx(0.0165)
    assert rad.rate(1e5) == 0.0 and rad.rate(900.0) < 0
    half = Radius.from_series(t, r_cm, shrink_scale=0.5)
    assert half.value(7200.0) == pytest.approx(0.02 - 0.5 * (0.02 - 0.0165))
    assert Radius.constant().is_constant and not rad.is_constant


# ---- solver ---------------------------------------------------------------------------------------


def test_grid_measures_sum_to_half() -> None:
    g = make_grid(20)
    assert g.vol.sum() == pytest.approx(0.5)
    assert g.vol[0] == pytest.approx(g.dxi**2 / 8)
    np.testing.assert_allclose(verify.control_volumes(g.xi), g.vol)


def test_equilibrium_is_a_steady_state() -> None:
    spec = ProblemSpec(APPENDIX2, Chamber.constant(T_INIT, C_INIT), Radius.constant())
    system = RadialSystem(spec, make_grid(10))
    y = np.concatenate([np.full(11, T_INIT), np.full(11, C_INIT)])
    np.testing.assert_allclose(system(0.0, y), 0.0, atol=1e-14)


def test_solver_matches_the_bessel_series_for_constant_properties() -> None:
    d_const = 5e-9
    props = replace(APPENDIX2, d_const=d_const)
    spec = ProblemSpec(props, Chamber.constant(50.0, 0.05), Radius.constant())
    t_out = np.array([0.0, 300.0, 900.0, 1800.0])
    sol = solve(spec, 40, t_out, rtol=1e-8, atol=1e-10)
    alpha = props.k0 / (props.rho0 * props.cp0)
    exact_t = verify.analytic_field(sol.xi, t_out, alpha, props.h, props.k0, T_INIT, 50.0)
    exact_c = verify.analytic_field(sol.xi, t_out, d_const, props.hm, d_const, C_INIT, 0.05)
    assert np.max(np.abs(sol.temp - exact_t)) < 2e-3
    assert np.max(np.abs(sol.moist - exact_c)) < 5e-3  # thin early boundary layer on the coarse test grid


def test_discrete_moisture_and_energy_balances_hold() -> None:
    spec = ProblemSpec(APPENDIX2, Chamber.constant(50.0, 0.05), Radius.constant())
    sol = solve(spec, 20, np.arange(0.0, 1801.0, 2.0))
    mb = verify.moisture_balance(
        sol.t, sol.xi, sol.moist, sol.radius, sol.air_hum, spec.props.hm, cum_loss=sol.cum_moist_loss
    )
    eb = verify.energy_balance(sol.t, sol.xi, sol.temp, R0, sol.air_temp, spec.props, cum_loss=sol.cum_heat_loss)
    assert mb["sampled"]["max_relative_residual"] < 1e-4 and eb["sampled"]["max_relative_residual"] < 1e-4
    assert mb["exact"]["max_relative_residual"] < 1e-8 and eb["exact"]["max_relative_residual"] < 1e-8
    ex = verify.extremum_checks(sol.temp, sol.moist, sol.air_temp, sol.air_hum, T_INIT, C_INIT)
    assert ex["temp_within_bounds"] and ex["moist_within_bounds"] and ex["moist_monotone_in_r"]


def test_moving_boundary_reduces_to_the_fixed_domain_for_constant_radius() -> None:
    chamber = Chamber.constant(50.0, 0.05)
    t_out = np.arange(0.0, 3601.0, 600.0)
    ref = solve(ProblemSpec(APPENDIX4, chamber, Radius.constant()), 20, t_out, rtol=1e-10, atol=1e-12)
    pseudo = Radius(np.array([0.0, 1.0, 2.0]), np.array([R0, R0, R0]), force_interpolant=True)
    for form in ("lagrangian", "eulerian"):
        test = solve(ProblemSpec(APPENDIX4, chamber, pseudo, formulation=form), 20, t_out, rtol=1e-10, atol=1e-12)
        assert np.max(np.abs(test.moist - ref.moist)) < 1e-8
        assert np.max(np.abs(test.temp - ref.temp)) < 1e-6


def test_lagrangian_shrinkage_conserves_dry_basis_water() -> None:
    t = np.array([0.0, 1800.0, 3600.0, 7200.0])
    rad = Radius.from_series(t, np.array([2.0, 1.8, 1.7, 1.6]))
    spec = ProblemSpec(APPENDIX4, Chamber.constant(50.0, 0.05), rad, formulation="lagrangian")
    sol = solve(spec, 20, np.arange(0.0, 7201.0, 5.0))
    mb = verify.moisture_balance(
        sol.t, sol.xi, sol.moist, sol.radius, sol.air_hum, spec.props.hm, cum_loss=sol.cum_moist_loss
    )
    assert mb["sampled"]["max_relative_residual"] < 2e-4
    assert mb["exact"]["max_relative_residual"] < 1e-8
    # the Eulerian alternative loses water through the receding boundary: its balance residual is not small
    alt = solve(replace(spec, formulation="eulerian"), 20, np.arange(0.0, 7201.0, 5.0))
    mb_alt = verify.moisture_balance(
        alt.t, alt.xi, alt.moist, alt.radius, alt.air_hum, spec.props.hm, cum_loss=alt.cum_moist_loss
    )
    assert mb_alt["exact"]["max_relative_residual"] > 1e-3


def test_solve_until_dry_detects_the_crossing_consistently() -> None:
    props = replace(APPENDIX2, d0=7e-7)  # fast drying for the test
    spec = ProblemSpec(props, Chamber.constant(50.0, 0.05), Radius.constant())
    sol = solve_until_dry(spec, 20, 60.0, horizon=3600.0)
    assert sol.t_dry is not None and sol.dry_profile is not None
    checks = verify.drying_checks(sol.t, sol.moist, sol.t_dry, sol.dry_profile)
    assert checks["rows_before_all_wet"] and checks["rows_after_all_dry"]
    assert checks["profile_max_error_at_t_dry"] < 1e-6
    assert sol.t[-1] == pytest.approx(60.0 * math.ceil(sol.t_dry / 60.0))


def test_subgrid_and_physical_sampling() -> None:
    idx = subgrid_indices(160)
    assert len(idx) == 21 and idx[1] == 8
    xi = np.linspace(0, 1, 41)
    profile = 1.0 - xi**2
    out = sample_physical(profile, xi, 0.015, np.array([0.0, 0.0075, 0.015, 0.016]))
    assert out[0] == pytest.approx(1.0) and out[1] == pytest.approx(0.75, abs=1e-9)
    assert out[2] == pytest.approx(0.0, abs=1e-9) and np.isnan(out[3])


def test_piecewise_integration_is_exact_for_piecewise_linear_forcing() -> None:
    """y' = -y + f(t) with f piecewise linear (kinks every 60 s): restarting at the kinks keeps BDF accurate."""
    from pipelines.a.solver import integrate_piecewise

    knots = np.arange(0.0, 601.0, 60.0)
    values = np.array([0, 1, 0, 2, 0, 3, 0, 4, 0, 5, 0], dtype=float)

    def fun(t: float, y: np.ndarray) -> np.ndarray:
        return -y + np.interp(t, knots, values)

    t_out = np.arange(0.0, 601.0, 30.0)
    states, stats, event_t, _ = integrate_piecewise(
        fun, np.array([0.0]), t_out, knots[1:], method="BDF", rtol=1e-10, atol=1e-12, jac_sparsity=None
    )
    # exact solution: y(t) = int_0^t e^{-(t-s)} f(s) ds, evaluated with a fine trapezoid rule on each segment
    fine = np.linspace(0.0, 600.0, 600001)
    f_fine = np.interp(fine, knots, values)
    exact = [np.trapezoid(np.exp(-(t - fine[fine <= t])) * f_fine[fine <= t], fine[fine <= t]) for t in t_out]
    assert stats["segments"] == 10 and event_t is None
    np.testing.assert_allclose(states[:, 0], exact, atol=1e-7)


# ---- verifiers ------------------------------------------------------------------------------------


@given(st.floats(min_value=0.05, max_value=50.0))
def test_robin_eigenvalues_satisfy_the_characteristic_equation(bi: float) -> None:
    lam = verify.robin_eigenvalues(bi, 6)
    assert np.all(np.diff(lam) > 0)
    np.testing.assert_allclose(lam * j1(lam) - bi * j0(lam), 0.0, atol=1e-9)


def test_analytic_series_has_the_right_limits() -> None:
    xi = np.linspace(0, 1, 11)
    theta = verify.analytic_robin_cylinder(xi, np.array([1e-4, 5.0]), 1.0)
    assert np.all(theta[0, :-2] > 0.99)  # essentially unchanged interior at tiny Fourier number
    assert np.all(np.abs(theta[1]) < 1e-3)  # fully equilibrated


def test_convergence_table_reports_second_order() -> None:
    levels = [{"k": k, "dr_cm": 0.1 / 2**k, "v": np.array([1.0 + (0.1 / 2**k) ** 2])} for k in range(4)]
    rows = verify.convergence_table(levels, "v")
    assert rows[0]["observed_order"] == pytest.approx(2.0, abs=0.2)


def test_axisymmetric_model_with_insulated_ends_reproduces_the_radial_scheme() -> None:
    chamber = Chamber.constant(50.0, 0.05)
    t_out = np.array([0.0, 300.0, 900.0])
    one_d = solve(ProblemSpec(APPENDIX3, chamber, Radius.constant()), 20, t_out)
    two_d = verify.solve_axisymmetric(APPENDIX3, chamber, t_out, nr=20, nz=3, end_bc="insulated")
    assert np.max(np.abs(two_d["moist"][:, :, 0] - one_d.moist)) < 1e-7
    assert np.max(np.abs(two_d["temp"][:, :, 0] - one_d.temp)) < 1e-5
    assert np.max(np.abs(two_d["moist"] - two_d["moist"][:, :, :1])) < 1e-9
