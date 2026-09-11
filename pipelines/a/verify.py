"""Independent checks of the drying solutions (MDR-0007). Nothing here calls the radial solver's internals.

* analytic Bessel-series solution of the constant-property Robin cylinder;
* discrete balances (moisture, energy) recomputed from the solution arrays only;
* extremum-principle, monotonicity and drying-criterion checks;
* grid-convergence bookkeeping (errors against the finest grid, observed orders);
* an independently written two-dimensional axisymmetric finite-volume model (physical coordinates)
  used both as a code cross-check (insulated ends) and to quantify end effects (Robin ends).
"""

from __future__ import annotations

from itertools import pairwise
from typing import Any

import numpy as np
from scipy import sparse
from scipy.optimize import brentq
from scipy.special import j0, j1, jn_zeros

from pipelines.a.physics import C_TARGET, LENGTH, R0, Chamber, Properties

# ----------------------------------------------------------------------------------------------
# 1. Analytic solution: infinite cylinder, uniform initial value, Robin boundary, constant ambient
# ----------------------------------------------------------------------------------------------


def robin_eigenvalues(bi: float, count: int) -> np.ndarray:
    """Roots of lambda J1(lambda) = Bi J0(lambda); the n-th root lies between consecutive zeros of J0."""
    zeros = np.concatenate([[0.0], jn_zeros(0, count)])
    roots = []
    for lo, hi in pairwise(zeros):
        f = lambda x: x * j1(x) - bi * j0(x)  # noqa: E731
        roots.append(brentq(f, lo + 1e-12, hi - 1e-12, xtol=1e-14, rtol=1e-14, maxiter=500))
    return np.asarray(roots)


def analytic_robin_cylinder(xi: np.ndarray, fourier: np.ndarray, bi: float, n_terms: int = 400) -> np.ndarray:
    """Dimensionless excess theta(xi, Fo) = (u - u_inf)/(u_0 - u_inf); returns shape (len(fourier), len(xi))."""
    lam = robin_eigenvalues(bi, n_terms)
    coeff = 2.0 * bi / ((lam**2 + bi**2) * j0(lam))
    basis = j0(np.outer(lam, xi))  # (n_terms, nxi)
    decay = np.exp(-np.outer(np.asarray(fourier, dtype=float), lam**2))  # (nt, n_terms)
    return (decay * coeff) @ basis


def analytic_field(
    xi: np.ndarray, t: np.ndarray, diffusivity: float, transfer: float, coeff: float, u0: float, u_inf: float
) -> np.ndarray:
    """Physical field u(xi, t) for a constant-property cylinder of radius R0 with Robin coefficient ``transfer``."""
    bi = transfer * R0 / coeff
    fo = diffusivity * np.asarray(t, dtype=float) / R0**2
    theta = analytic_robin_cylinder(xi, fo, bi)
    theta[np.asarray(t) == 0.0] = 1.0  # the series converges slowly at Fo = 0; the initial value is exact
    return u_inf + (u0 - u_inf) * theta


# ----------------------------------------------------------------------------------------------
# 2. Balances and qualitative properties computed from output arrays
# ----------------------------------------------------------------------------------------------


def control_volumes(xi: np.ndarray) -> np.ndarray:
    faces = 0.5 * (xi[:-1] + xi[1:])
    lower = np.concatenate([[0.0], faces])
    upper = np.concatenate([faces, [1.0]])
    return 0.5 * (upper**2 - lower**2)


def _balance(lhs: np.ndarray, cumulative: np.ndarray) -> dict[str, float]:
    residual = lhs - cumulative
    scale = max(abs(lhs[-1]), 1e-300)
    return {
        "max_abs_residual": float(np.max(np.abs(residual))),
        "final_relative_residual": float(abs(residual[-1]) / scale),
        "max_relative_residual": float(np.max(np.abs(residual)) / scale),
    }


def moisture_balance(
    t: np.ndarray,
    xi: np.ndarray,
    moist: np.ndarray,
    radius: np.ndarray,
    air_hum: np.ndarray,
    hm: float,
    cum_loss: np.ndarray | None = None,
) -> dict[str, Any]:
    """Discrete water inventory W = sum_i vol_i C_i versus the time-integrated surface flux.

    Fixed radius: d/dt (R^2 W) = -R hm (C_s - C_air); shrinking (Lagrangian): d/dt W = -hm (C_s - C_air)/R(t).
    ``sampled``: flux from the output samples (trapezoid rule, limited by the output spacing);
    ``exact``: flux accumulated inside the integrator (``cum_loss`` = int hm (C_s - C_air)/R dt), which turns the
    balance into a linear invariant of the ODE system and therefore tests the scheme itself.
    """
    vol = control_volumes(xi)
    inventory = moist @ vol
    flux = hm * (moist[:, -1] - air_hum)
    if np.allclose(radius, radius[0]):
        lhs = radius[0] ** 2 * (inventory - inventory[0])
        integrand = -radius[0] * flux
        exact_scale = radius[0] ** 2
    else:
        lhs = inventory - inventory[0]
        integrand = -flux / radius
        exact_scale = 1.0
    cumulative = np.concatenate([[0.0], np.cumsum(0.5 * (integrand[1:] + integrand[:-1]) * np.diff(t))])
    report: dict[str, Any] = {
        "inventory_initial": float(inventory[0]),
        "inventory_final": float(inventory[-1]),
        "loss_fraction": float(1.0 - inventory[-1] / inventory[0]),
        "sampled": _balance(lhs, cumulative),
    }
    if cum_loss is not None:
        report["exact"] = _balance(lhs, -exact_scale * cum_loss)
    best = report.get("exact", report["sampled"])
    report["max_relative_residual"] = best["max_relative_residual"]
    report["final_relative_residual"] = best["final_relative_residual"]
    return report


def energy_balance(
    t: np.ndarray,
    xi: np.ndarray,
    temp: np.ndarray,
    radius: float,
    air_temp: np.ndarray,
    props: Properties,
    cum_loss: np.ndarray | None = None,
) -> dict[str, Any]:
    """Constant-property energy balance: rho cp R^2 d/dt sum vol_i T_i = -R h (T_s - T_air)."""
    vol = control_volumes(xi)
    rho_cp = float(props.rho(np.array(0.0)) * props.cp(np.array(0.0)))
    inventory = rho_cp * radius**2 * (temp @ vol)
    integrand = -radius * props.h * (temp[:, -1] - air_temp)
    cumulative = np.concatenate([[0.0], np.cumsum(0.5 * (integrand[1:] + integrand[:-1]) * np.diff(t))])
    report: dict[str, Any] = {"sampled": _balance(inventory - inventory[0], cumulative)}
    if cum_loss is not None:  # rho cp R^2 dW/dt = -R h (T_s - T_air) and Q_h = int h (T_s - T_air)/R dt
        report["exact"] = _balance(inventory - inventory[0], -(radius**2) * cum_loss)
    best = report.get("exact", report["sampled"])
    report["max_relative_residual"] = best["max_relative_residual"]
    report["final_relative_residual"] = best["final_relative_residual"]
    return report


def extremum_checks(
    temp: np.ndarray, moist: np.ndarray, air_temp: np.ndarray, air_hum: np.ndarray, t0: float, c0: float
) -> dict[str, Any]:
    t_lo, t_hi = min(t0, float(air_temp.min())), max(t0, float(air_temp.max()))
    c_lo, c_hi = min(c0, float(air_hum.min())), max(c0, float(air_hum.max()))
    tol = 1e-9
    dc = np.diff(moist, axis=1)  # along r: should be <= 0 (drying from the surface)
    return {
        "temp_min": float(temp.min()),
        "temp_max": float(temp.max()),
        "temp_bounds": [t_lo, t_hi],
        "temp_within_bounds": bool(temp.min() >= t_lo - tol and temp.max() <= t_hi + tol),
        "moist_min": float(moist.min()),
        "moist_max": float(moist.max()),
        "moist_bounds": [c_lo, c_hi],
        "moist_within_bounds": bool(moist.min() >= c_lo - tol and moist.max() <= c_hi + tol),
        "moist_monotone_in_r": bool(np.all(dc <= 1e-8)),
        "moist_max_positive_slope": float(max(dc.max(), 0.0)),
        "moist_argmax_is_centre": bool(np.all(moist.max(axis=1) - moist[:, 0] <= 1e-9)),  # ties up to roundoff
    }


def drying_checks(
    t: np.ndarray, moist: np.ndarray, t_dry: float, dry_profile: np.ndarray, target: float = C_TARGET
) -> dict[str, Any]:
    """The exact drying instant must separate 'not yet dry' rows from 'dry' rows on the output grid."""
    max_c = moist.max(axis=1)
    before = t < t_dry
    return {
        "t_dry": float(t_dry),
        "t_dry_hours": float(t_dry / 3600.0),
        "profile_max_at_t_dry": float(dry_profile.max()),
        "profile_max_error_at_t_dry": float(abs(dry_profile.max() - target)),
        "rows_before_all_wet": bool(np.all(max_c[before] >= target)),
        "rows_after_all_dry": bool(np.all(max_c[~before] < target)),
        "last_row_time": float(t[-1]),
        "last_row_max": float(max_c[-1]),
    }


# ----------------------------------------------------------------------------------------------
# 3. Convergence bookkeeping
# ----------------------------------------------------------------------------------------------


def convergence_table(levels: list[dict[str, Any]], key: str) -> list[dict[str, Any]]:
    """levels: [{"k": k, "dr_cm": ..., key: array}] with arrays sampled on common points; finest is last."""
    reference = levels[-1][key]
    rows: list[dict[str, Any]] = []
    errors: list[float] = []
    for level in levels[:-1]:
        err = float(np.max(np.abs(level[key] - reference)))
        errors.append(err)
        rows.append({"k": level["k"], "dr_cm": level["dr_cm"], "max_abs_error": err})
    for i, row in enumerate(rows):
        if i + 1 < len(rows) and errors[i + 1] > 0:
            row["observed_order"] = float(np.log2(errors[i] / errors[i + 1]))
        else:
            row["observed_order"] = None
    return rows


def richardson_estimate(coarse_err: float, fine_err: float, order: float = 2.0) -> float:
    """Error estimate of the fine level from two successive errors against a common finer reference."""
    return float(abs(fine_err) / (2**order - 1))


# ----------------------------------------------------------------------------------------------
# 4. Independent two-dimensional axisymmetric finite-volume model (fixed radius, physical coordinates)
# ----------------------------------------------------------------------------------------------


class AxisymmetricSystem:
    """Vertex-centred FV on (r, z) in [0, R0] x [0, L/2]; z = 0 is the mid-plane (symmetry)."""

    def __init__(self, props: Properties, chamber: Chamber, nr: int, nz: int, end_bc: str = "robin") -> None:
        self.props = props
        self.chamber = chamber
        self.nr, self.nz = nr, nz
        self.end_bc = end_bc
        self.dr = R0 / nr
        self.dz = (LENGTH / 2.0) / nz
        r = np.linspace(0.0, R0, nr + 1)
        faces = 0.5 * (r[:-1] + r[1:])
        self.r_face = faces
        lower = np.concatenate([[0.0], faces])
        upper = np.concatenate([faces, [R0]])
        self.m_r = 0.5 * (upper**2 - lower**2)  # radial measure per node
        self.w_z = np.full(nz + 1, self.dz)
        self.w_z[0] = self.w_z[-1] = 0.5 * self.dz  # half cells at the mid-plane and the end face
        self.shape = (nr + 1, nz + 1)
        self.m = (nr + 1) * (nz + 1)

    def __call__(self, t: float, y: np.ndarray) -> np.ndarray:
        p = self.props
        temp = y[: self.m].reshape(self.shape)
        moist = y[self.m :].reshape(self.shape)
        t_air = float(self.chamber.air_temperature(t))
        c_air = float(self.chamber.air_humidity(t))
        # radial faces (between i and i+1), shape (nr, nz+1)
        cf = 0.5 * (moist[:-1, :] + moist[1:, :])
        tf = 0.5 * (temp[:-1, :] + temp[1:, :])
        fc_r = self.r_face[:, None] * p.diffusivity(cf, tf) * (moist[1:, :] - moist[:-1, :]) / self.dr
        ft_r = self.r_face[:, None] * p.k(cf) * (temp[1:, :] - temp[:-1, :]) / self.dr
        # axial faces (between j and j+1), shape (nr+1, nz)
        cz = 0.5 * (moist[:, :-1] + moist[:, 1:])
        tz = 0.5 * (temp[:, :-1] + temp[:, 1:])
        fc_z = p.diffusivity(cz, tz) * (moist[:, 1:] - moist[:, :-1]) / self.dz
        ft_z = p.k(cz) * (temp[:, 1:] - temp[:, :-1]) / self.dz
        div_c = np.zeros(self.shape)
        div_t = np.zeros(self.shape)
        div_c[:-1, :] += fc_r
        div_c[1:, :] -= fc_r
        div_t[:-1, :] += ft_r
        div_t[1:, :] -= ft_r
        div_c[-1, :] += -R0 * p.hm * (moist[-1, :] - c_air)
        div_t[-1, :] += -R0 * p.h * (temp[-1, :] - t_air)
        div_c *= self.w_z[None, :]
        div_t *= self.w_z[None, :]
        ax_c = np.zeros(self.shape)
        ax_t = np.zeros(self.shape)
        ax_c[:, :-1] += fc_z
        ax_c[:, 1:] -= fc_z
        ax_t[:, :-1] += ft_z
        ax_t[:, 1:] -= ft_z
        if self.end_bc == "robin":
            ax_c[:, -1] += -p.hm * (moist[:, -1] - c_air)
            ax_t[:, -1] += -p.h * (temp[:, -1] - t_air)
        div_c += ax_c * self.m_r[:, None]
        div_t += ax_t * self.m_r[:, None]
        measure = self.m_r[:, None] * self.w_z[None, :]
        dc = div_c / measure
        dt = div_t / (measure * p.rho(moist) * p.cp(moist))
        return np.concatenate([dt.ravel(), dc.ravel()])

    def jacobian_sparsity(self) -> sparse.csr_matrix:
        nr1, nz1 = self.shape
        idx = np.arange(self.m).reshape(self.shape)
        pairs = [(idx, idx)]
        pairs.append((idx[:-1, :], idx[1:, :]))
        pairs.append((idx[1:, :], idx[:-1, :]))
        pairs.append((idx[:, :-1], idx[:, 1:]))
        pairs.append((idx[:, 1:], idx[:, :-1]))
        rows = np.concatenate([a.ravel() for a, _ in pairs])
        cols = np.concatenate([b.ravel() for _, b in pairs])
        blocks = []
        for ro in (0, self.m):
            for co in (0, self.m):
                blocks.append((rows + ro, cols + co))
        r_all = np.concatenate([b[0] for b in blocks])
        c_all = np.concatenate([b[1] for b in blocks])
        return sparse.csr_matrix((np.ones(len(r_all)), (r_all, c_all)), shape=(2 * self.m, 2 * self.m))


def solve_axisymmetric(
    props: Properties,
    chamber: Chamber,
    t_out: np.ndarray,
    nr: int = 20,
    nz: int = 25,
    end_bc: str = "robin",
    t_init: float = 28.0,
    c_init: float = 2.55,
    rtol: float = 1e-7,
    atol: float = 1e-9,
    detect_dry: bool = False,
    method: str = "BDF",
) -> dict[str, Any]:
    system = AxisymmetricSystem(props, chamber, nr, nz, end_bc)
    y0 = np.concatenate([np.full(system.m, t_init), np.full(system.m, c_init)])

    def event(t: float, y: np.ndarray) -> float:
        return float(np.max(y[system.m :]) - C_TARGET)

    event.terminal = False  # type: ignore[attr-defined]
    event.direction = -1  # type: ignore[attr-defined]
    from pipelines.a.solver import integrate_piecewise  # generic ODE driver only (not the radial scheme)

    states, stats, t_dry, _ = integrate_piecewise(
        system,
        y0,
        np.asarray(t_out, dtype=float),
        chamber.knots(),
        method=method,
        rtol=rtol,
        atol=atol,
        jac_sparsity=system.jacobian_sparsity(),
        event=event if detect_dry else None,
    )
    temp = states[:, : system.m].reshape(len(t_out), *system.shape)
    moist = states[:, system.m :].reshape(len(t_out), *system.shape)
    out: dict[str, Any] = {
        "t": np.asarray(t_out, dtype=float),
        "r": np.linspace(0.0, R0, nr + 1),
        "z": np.linspace(0.0, LENGTH / 2.0, nz + 1),
        "temp": temp,
        "moist": moist,
        "nfev": stats["nfev"],
        "nlu": stats["nlu"],
    }
    if t_dry is not None:
        out["t_dry"] = t_dry
    return out
