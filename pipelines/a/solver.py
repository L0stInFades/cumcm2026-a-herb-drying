"""Radial finite-volume solver for the coupled heat/moisture equations on a (possibly shrinking) cylinder.

Discretisation (MDR-0004): vertex-centred finite volumes on the Landau coordinate xi = r/R(t) in [0, 1],
nodes xi_i = i/n, control volumes [xi_{i-1/2}, xi_{i+1/2}] with the radial measure (xi^2 differences)/2,
Robin conditions applied as face fluxes on the outer half-cell, symmetry at the axis. Time integration
by scipy's adaptive BDF (or Radau) with a block-tridiagonal Jacobian sparsity pattern; every requested
output instant is sampled from the integrator's dense interpolant.

Two formulations for a moving boundary (MDR-0006):
  * ``lagrangian`` - material coordinate, dry matter conserved, no convective term, 1/R(t)^2 scaling;
  * ``eulerian``   - laboratory-frame PDE on 0 < r < R(t); the transformation adds +xi (Rdot/R) d/dxi.
For a constant radius both coincide with the fixed-domain scheme.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any

import numpy as np
from scipy import sparse
from scipy.integrate import solve_ivp

from pipelines.a.physics import C_INIT, C_TARGET, T_INIT, Chamber, Properties, Radius

FACE_SCHEMES = ("midpoint", "arithmetic", "harmonic")
FORMULATIONS = ("lagrangian", "eulerian")


@dataclass(frozen=True)
class Grid:
    n: int
    xi: np.ndarray  # nodes, n+1
    dxi: float
    xi_face: np.ndarray  # interior faces xi_{i+1/2}, n values
    vol: np.ndarray  # control-volume measure (xi_{i+1/2}^2 - xi_{i-1/2}^2)/2, n+1 values


def make_grid(n: int) -> Grid:
    xi = np.linspace(0.0, 1.0, n + 1)
    faces = 0.5 * (xi[:-1] + xi[1:])
    lower = np.concatenate([[0.0], faces])
    upper = np.concatenate([faces, [1.0]])
    return Grid(n, xi, 1.0 / n, faces, 0.5 * (upper**2 - lower**2))


@dataclass(frozen=True)
class ProblemSpec:
    """Everything that defines one initial-boundary-value problem."""

    props: Properties
    chamber: Chamber
    radius: Radius
    formulation: str = "lagrangian"
    face_scheme: str = "midpoint"
    t_init: float = T_INIT
    c_init: float = C_INIT

    def describe(self) -> dict[str, Any]:
        return {
            "properties": self.props.describe(),
            "chamber": self.chamber.describe(),
            "radius": self.radius.describe(),
            "formulation": self.formulation,
            "face_scheme": self.face_scheme,
            "t_init": self.t_init,
            "c_init": self.c_init,
        }


@dataclass
class Solution:
    t: np.ndarray  # output instants, s
    xi: np.ndarray  # nodes
    temp: np.ndarray  # (len(t), n+1) degC
    moist: np.ndarray  # (len(t), n+1) kg/kg
    radius: np.ndarray  # R(t) at output instants, m
    air_temp: np.ndarray
    air_hum: np.ndarray
    stats: dict[str, Any] = field(default_factory=dict)
    cum_moist_loss: np.ndarray | None = None  # integral of hm (C_s - C_air)/R dt, integrated as an ODE state
    cum_heat_loss: np.ndarray | None = None  # integral of h (T_s - T_air)/R dt
    t_dry: float | None = None  # exact instant at which max_i C_i first drops below C_TARGET
    dry_profile: np.ndarray | None = None  # moisture profile at t_dry
    dry_temp: np.ndarray | None = None

    @property
    def r(self) -> np.ndarray:
        """Physical node positions (len(t), n+1) in m."""
        return np.outer(self.radius, self.xi)

    def surface_moisture_flux(self, props: Properties) -> np.ndarray:
        """h_m (C_s - C_air) at the output instants, (kg/kg) m/s."""
        return props.hm * (self.moist[:, -1] - self.air_hum)

    def surface_heat_flux(self, props: Properties) -> np.ndarray:
        return props.h * (self.temp[:, -1] - self.air_temp)


class RadialSystem:
    """Right-hand side f(t, y) of the semi-discrete system, y = [T_0..T_n, C_0..C_n, Q_m, Q_h].

    The two trailing states accumulate the boundary losses hm (C_s - C_air)/R and h (T_s - T_air)/R so that
    the discrete balances sum_i vol_i C_i + Q_m = const (and its thermal analogue for constant properties)
    are linear invariants of the ODE system, which BDF/Radau preserve to integrator accuracy.
    """

    def __init__(self, spec: ProblemSpec, grid: Grid) -> None:
        if spec.face_scheme not in FACE_SCHEMES:
            raise ValueError(f"face_scheme must be one of {FACE_SCHEMES}")
        if spec.formulation not in FORMULATIONS:
            raise ValueError(f"formulation must be one of {FORMULATIONS}")
        self.spec = spec
        self.grid = grid
        self.n1 = grid.n + 1
        self.convective = spec.formulation == "eulerian" and not spec.radius.is_constant

    def face_coefficients(self, temp: np.ndarray, moist: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        props = self.spec.props
        scheme = self.spec.face_scheme
        if scheme == "midpoint":
            cf = 0.5 * (moist[:-1] + moist[1:])
            tf = 0.5 * (temp[:-1] + temp[1:])
            return props.diffusivity(cf, tf), props.k(cf)
        dn = props.diffusivity(moist, temp)
        kn = props.k(moist)
        if scheme == "arithmetic":
            return 0.5 * (dn[:-1] + dn[1:]), 0.5 * (kn[:-1] + kn[1:])
        return 2.0 * dn[:-1] * dn[1:] / (dn[:-1] + dn[1:]), 2.0 * kn[:-1] * kn[1:] / (kn[:-1] + kn[1:])

    def __call__(self, t: float, y: np.ndarray) -> np.ndarray:
        g = self.grid
        n1 = self.n1
        temp = y[:n1]
        moist = y[n1 : 2 * n1]
        props = self.spec.props
        radius = float(self.spec.radius.value(t))
        r2 = radius * radius
        t_air = float(self.spec.chamber.air_temperature(t))
        c_air = float(self.spec.chamber.air_humidity(t))
        d_face, k_face = self.face_coefficients(temp, moist)
        flux_c = g.xi_face * d_face * (moist[1:] - moist[:-1]) / g.dxi  # xi D dC/dxi at interior faces
        flux_t = g.xi_face * k_face * (temp[1:] - temp[:-1]) / g.dxi
        div_c = np.empty(n1)
        div_t = np.empty(n1)
        div_c[:-1] = flux_c
        div_t[:-1] = flux_t
        div_c[-1] = -radius * props.hm * (moist[-1] - c_air)  # outer face: xi D dC/dxi = R (-hm (C_s - C_air))
        div_t[-1] = -radius * props.h * (temp[-1] - t_air)
        div_c[1:] -= flux_c
        div_t[1:] -= flux_t
        dc = div_c / (r2 * g.vol)
        dt = div_t / (r2 * g.vol * props.rho(moist) * props.cp(moist))
        if self.convective:
            rate = float(self.spec.radius.rate(t)) / radius
            grad_c = np.empty(n1)
            grad_t = np.empty(n1)
            grad_c[0] = 0.0
            grad_t[0] = 0.0
            grad_c[1:-1] = (moist[2:] - moist[:-2]) / (2.0 * g.dxi)
            grad_t[1:-1] = (temp[2:] - temp[:-2]) / (2.0 * g.dxi)
            grad_c[-1] = (3.0 * moist[-1] - 4.0 * moist[-2] + moist[-3]) / (2.0 * g.dxi)
            grad_t[-1] = (3.0 * temp[-1] - 4.0 * temp[-2] + temp[-3]) / (2.0 * g.dxi)
            dc = dc + g.xi * rate * grad_c
            dt = dt + g.xi * rate * grad_t
        losses = [props.hm * (moist[-1] - c_air) / radius, props.h * (temp[-1] - t_air) / radius]
        return np.concatenate([dt, dc, losses])

    def jacobian_sparsity(self) -> sparse.csr_matrix:
        n1 = self.n1
        rows: list[int] = []
        cols: list[int] = []
        for i in range(n1):
            neighbours = [j for j in (i - 1, i, i + 1) if 0 <= j < n1]
            if self.convective and i == n1 - 1:
                neighbours.append(n1 - 3)
            for j in neighbours:
                for row_off in (0, n1):
                    for col_off in (0, n1):
                        rows.append(i + row_off)
                        cols.append(j + col_off)
        rows += [2 * n1, 2 * n1 + 1]
        cols += [2 * n1 - 1, n1 - 1]  # Q_m depends on C_n, Q_h on T_n
        data = np.ones(len(rows))
        return sparse.csr_matrix((data, (rows, cols)), shape=(2 * n1 + 2, 2 * n1 + 2))


def _dry_event(n1: int, target: float) -> Any:
    def event(t: float, y: np.ndarray) -> float:
        return float(np.max(y[n1 : 2 * n1]) - target)

    event.terminal = False  # type: ignore[attr-defined]
    event.direction = -1  # type: ignore[attr-defined]
    return event


def solve(
    spec: ProblemSpec,
    n: int,
    t_out: np.ndarray,
    *,
    method: str = "BDF",
    rtol: float = 1e-7,
    atol: float = 1e-9,
    detect_dry: bool = False,
    target: float = C_TARGET,
    max_step: float = np.inf,
) -> Solution:
    """Integrate the semi-discrete system and sample the dense solution at ``t_out`` (must start at 0)."""
    grid = make_grid(n)
    system = RadialSystem(spec, grid)
    t_out = np.asarray(t_out, dtype=float)
    if t_out[0] != 0.0 or np.any(np.diff(t_out) <= 0):
        raise ValueError("t_out must start at 0 and be strictly increasing")
    y0 = np.concatenate([np.full(grid.n + 1, spec.t_init), np.full(grid.n + 1, spec.c_init), [0.0, 0.0]])
    events = [_dry_event(grid.n + 1, target)] if detect_dry else None
    result = solve_ivp(
        system,
        (0.0, float(t_out[-1])),
        y0,
        method=method,
        t_eval=t_out,
        rtol=rtol,
        atol=atol,
        jac_sparsity=system.jacobian_sparsity(),
        dense_output=detect_dry,
        events=events,
        max_step=max_step,
    )
    if not result.success:
        raise RuntimeError(f"solve_ivp failed: {result.message}")
    n1 = grid.n + 1
    temp = result.y[:n1].T.copy()
    moist = result.y[n1 : 2 * n1].T.copy()
    radius = np.asarray(spec.radius.value(t_out), dtype=float)
    sol = Solution(
        t=t_out,
        xi=grid.xi,
        temp=temp,
        moist=moist,
        radius=radius,
        air_temp=np.asarray(spec.chamber.air_temperature(t_out), dtype=float),
        air_hum=np.asarray(spec.chamber.air_humidity(t_out), dtype=float),
        stats={
            "method": method,
            "rtol": rtol,
            "atol": atol,
            "n_intervals": n,
            "nfev": int(result.nfev),
            "njev": int(result.njev),
            "nlu": int(result.nlu),
            "message": str(result.message),
        },
        cum_moist_loss=result.y[2 * n1].copy(),
        cum_heat_loss=result.y[2 * n1 + 1].copy(),
    )
    if detect_dry and result.t_events and len(result.t_events[0]):
        t_dry = float(result.t_events[0][0])
        state = result.sol(t_dry)
        sol.t_dry = t_dry
        sol.dry_temp = np.asarray(state[:n1], dtype=float)
        sol.dry_profile = np.asarray(state[n1 : 2 * n1], dtype=float)
    return sol


def solve_until_dry(
    spec: ProblemSpec,
    n: int,
    dt_out: float = 60.0,
    *,
    horizon: float = 5 * 86400.0,
    max_horizon: float = 30 * 86400.0,
    **kwargs: Any,
) -> Solution:
    """Integrate until max_i C_i < C_TARGET, keep output rows up to the first ``dt_out`` multiple after it."""
    while True:
        t_out = np.arange(0.0, horizon + 0.5 * dt_out, dt_out)
        sol = solve(spec, n, t_out, detect_dry=True, **kwargs)
        if sol.t_dry is not None:
            break
        horizon *= 2.0
        if horizon > max_horizon:
            raise RuntimeError("drying criterion not met within the maximum horizon")
    last = math.ceil(sol.t_dry / dt_out - 1e-9)
    keep = slice(0, last + 1)
    sol.t = sol.t[keep]
    sol.temp = sol.temp[keep]
    sol.moist = sol.moist[keep]
    sol.radius = sol.radius[keep]
    sol.air_temp = sol.air_temp[keep]
    sol.air_hum = sol.air_hum[keep]
    if sol.cum_moist_loss is not None and sol.cum_heat_loss is not None:
        sol.cum_moist_loss = sol.cum_moist_loss[keep]
        sol.cum_heat_loss = sol.cum_heat_loss[keep]
    sol.stats["t_dry_grid"] = float(last * dt_out)
    return sol


def subgrid_indices(n: int, n_out: int = 20) -> np.ndarray:
    """Indices of the nodes coinciding with the coarse output grid (n must be an integer multiple of n_out)."""
    if n % n_out:
        raise ValueError(f"n={n} is not a multiple of the output grid {n_out}")
    return np.arange(0, n + 1, n // n_out)


def sample_physical(profile: np.ndarray, xi: np.ndarray, radius: float, r_targets: np.ndarray) -> np.ndarray:
    """Cubic-spline sample of a nodal profile at physical radii; NaN where r > R(t)."""
    from scipy.interpolate import CubicSpline

    spline = CubicSpline(xi, profile, bc_type=((1, 0.0), "not-a-knot"))
    xi_t = np.asarray(r_targets, dtype=float) / radius
    out = np.full(xi_t.shape, np.nan)
    inside = xi_t <= 1.0 + 1e-12
    out[inside] = spline(np.minimum(xi_t[inside], 1.0))
    return out
