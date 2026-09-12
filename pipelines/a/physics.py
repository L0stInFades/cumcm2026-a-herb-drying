"""Physical laws and boundary data of the herb-drying problem (SI units; temperatures in degC).

Property laws follow the problem's appendices (MDR-0001/0002), the drying-chamber history and its
constant-temperature plateau follow MDR-0003, the shrinking radius follows MDR-0006.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Any

import numpy as np
from scipy.interpolate import PchipInterpolator

R0 = 0.02  # initial radius, m
LENGTH = 0.25  # cylinder length, m
T_INIT = 28.0  # degC
C_INIT = 2.55  # kg/kg (dry basis)
C_TARGET = 0.15  # kg/kg, drying finished when every point is below this value
KELVIN = 273.15
C_FLOOR = 1e-6  # only used inside D(C, T) to avoid division by zero
T_DATA_END = 14400.0  # last chamber sample, s
PLATEAU_WINDOW = 3600.0  # averaging window for the constant-temperature plateau, s


@dataclass(frozen=True)
class Properties:
    """rho = rho0 + rho1*C; cp = cp0 + cp1*C/(C+1); k = k0 + k1*C/(C+1); D = d0*exp(-da/C)*exp(-de/T_K)."""

    name: str
    rho0: float
    rho1: float
    cp0: float
    cp1: float
    k0: float
    k1: float
    d0: float
    da: float
    de: float
    h: float = 25.0  # W/(m^2 K), appendix 2 (carried over to problems 2-4, MDR-0002)
    hm: float = 8e-7  # m/s, appendix 2
    d_const: float | None = None  # constant-D override (analytic verification only)

    def rho(self, c: np.ndarray) -> np.ndarray:
        return self.rho0 + self.rho1 * c

    def cp(self, c: np.ndarray) -> np.ndarray:
        return self.cp0 + self.cp1 * c / (c + 1.0)

    def k(self, c: np.ndarray) -> np.ndarray:
        return self.k0 + self.k1 * c / (c + 1.0)

    def diffusivity(self, c: np.ndarray, t_degc: np.ndarray) -> np.ndarray:
        if self.d_const is not None:
            return np.full_like(np.asarray(c, dtype=float), self.d_const)
        cs = np.maximum(c, C_FLOOR)
        out = self.d0 * np.exp(-self.da / cs)
        if self.de:
            out = out * np.exp(-self.de / (np.asarray(t_degc, dtype=float) + KELVIN))
        return out

    def scaled(self, **factors: float) -> Properties:
        """Return a copy with the named fields multiplied by the given factors (sensitivity runs)."""
        return replace(self, **{key: getattr(self, key) * value for key, value in factors.items()})

    def describe(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "rho": f"{self.rho0} + {self.rho1} C",
            "cp": f"{self.cp0} + {self.cp1} C/(C+1)",
            "k": f"{self.k0} + {self.k1} C/(C+1)",
            "D": f"{self.d0} exp(-{self.da}/C) exp(-{self.de}/T_K)" if self.d_const is None else f"{self.d_const}",
            "h": self.h,
            "hm": self.hm,
        }


APPENDIX2 = Properties("appendix2", 820.0, 0.0, 2600.0, 0.0, 0.36, 0.0, 7e-9, 0.89, 0.0)
APPENDIX3 = Properties("appendix3", 650.0, 128.0, 1450.0, 2736.0, 0.21, 0.38, 2.4e-3, 0.45, 3850.0)
APPENDIX4 = Properties("appendix4", 760.0, 90.0, 1850.0, 2150.0, 0.12, 0.20, 4.2e-4, 0.30, 3850.0)
PROPERTIES = {"appendix2": APPENDIX2, "appendix3": APPENDIX3, "appendix4": APPENDIX4}


@dataclass(frozen=True)
class Chamber:
    """Drying-chamber temperature and humidity: linear interpolation of attachment 1, constant plateau after it."""

    t: np.ndarray
    temp: np.ndarray
    hum: np.ndarray
    temp_plateau: float
    hum_plateau: float
    temp_scale: float = 1.0  # multiplies the plateau only (sensitivity)
    hum_scale: float = 1.0  # multiplies the whole humidity history (sensitivity)

    @classmethod
    def from_series(cls, t: np.ndarray, temp: np.ndarray, hum: np.ndarray, plateau: str = "mean_last_hour") -> Chamber:
        t = np.asarray(t, dtype=float)
        temp = np.asarray(temp, dtype=float)
        hum = np.asarray(hum, dtype=float)
        if plateau == "mean_last_hour":
            window = t >= t[-1] - PLATEAU_WINDOW
            tp, hp = float(temp[window].mean()), float(hum[window].mean())
        elif plateau == "last_sample":
            tp, hp = float(temp[-1]), float(hum[-1])
        elif plateau == "nominal":
            tp, hp = 50.0, 0.05
        else:
            raise ValueError(f"unknown plateau rule {plateau!r}")
        return cls(t, temp, hum, tp, hp)

    @classmethod
    def constant(cls, temp: float, hum: float) -> Chamber:
        return cls(np.array([0.0, 1.0]), np.array([temp, temp]), np.array([hum, hum]), temp, hum)

    def air_temperature(self, t: Any) -> Any:
        return np.interp(t, self.t, self.temp, right=self.temp_plateau * self.temp_scale)

    def air_humidity(self, t: Any) -> Any:
        return self.hum_scale * np.interp(t, self.t, self.hum, right=self.hum_plateau)

    def knots(self) -> np.ndarray:
        """Instants where the forcing is not smooth (sample times incl. the switch to the plateau)."""
        return self.t[1:] if len(self.t) > 2 else np.array([])

    def describe(self) -> dict[str, Any]:
        return {
            "samples": len(self.t),
            "t_end": float(self.t[-1]),
            "temp_plateau": self.temp_plateau * self.temp_scale,
            "hum_plateau": self.hum_plateau * self.hum_scale,
            "temp_scale": self.temp_scale,
            "hum_scale": self.hum_scale,
        }


class Radius:
    """Radius history R(t): monotone PCHIP through attachment 2 (m), held at its last value afterwards."""

    def __init__(
        self, t: np.ndarray, r: np.ndarray, shrink_scale: float = 1.0, *, force_interpolant: bool = False
    ) -> None:
        self.t = np.asarray(t, dtype=float)
        self.r = np.asarray(r, dtype=float)
        self.shrink_scale = shrink_scale
        self.is_constant = bool(np.all(self.r == self.r[0])) and not force_interpolant
        self._t0, self._t1 = float(self.t[0]), float(self.t[-1])
        if not self.is_constant:
            self._pchip = PchipInterpolator(self.t, self.r, extrapolate=False)
            self._deriv = self._pchip.derivative()

    @classmethod
    def from_series(cls, t: np.ndarray, r_cm: np.ndarray, shrink_scale: float = 1.0) -> Radius:
        r = R0 - shrink_scale * (R0 - np.asarray(r_cm, dtype=float) / 100.0)
        return cls(np.asarray(t, dtype=float), r, shrink_scale)

    @classmethod
    def constant(cls, r: float = R0) -> Radius:
        return cls(np.array([0.0, 1.0]), np.array([r, r]), 0.0)

    def value(self, t: Any) -> Any:
        if self.is_constant:
            return np.full_like(np.asarray(t, dtype=float), self.r[0]) if np.ndim(t) else float(self.r[0])
        return self._pchip(np.clip(t, self._t0, self._t1))

    def knots(self) -> np.ndarray:
        """PCHIP knots (curvature jumps) including the switch to the constant tail."""
        return np.array([]) if self.is_constant else self.t[1:]

    def rate(self, t: Any) -> Any:
        if self.is_constant:
            return np.zeros_like(np.asarray(t, dtype=float)) if np.ndim(t) else 0.0
        tt = np.asarray(t, dtype=float)
        inside = (tt > self._t0) & (tt < self._t1)
        return np.where(inside, self._deriv(np.clip(tt, self._t0, self._t1)), 0.0)

    def describe(self) -> dict[str, Any]:
        return {
            "samples": len(self.t),
            "t_end": float(self.t[-1]),
            "r_start_cm": float(self.r[0] * 100),
            "r_end_cm": float(self.r[-1] * 100),
            "shrink_scale": self.shrink_scale,
        }


def thermal_diffusivity(props: Properties, c: float) -> float:
    cc = np.asarray(c, dtype=float)
    return float(props.k(cc) / (props.rho(cc) * props.cp(cc)))


def biot_numbers(props: Properties, c: float, t_degc: float, r: float = R0) -> dict[str, float]:
    cc = np.asarray(c, dtype=float)
    return {
        "Bi_heat": float(props.h * r / props.k(cc)),
        "Bi_mass": float(props.hm * r / props.diffusivity(cc, np.asarray(t_degc))),
        "alpha": thermal_diffusivity(props, c),
        "D": float(props.diffusivity(cc, np.asarray(t_degc))),
    }


# ----------------------------------------------------------------------------------------------
# Psychrometrics and the energy-limited evaporation cap (model-evaluation extension, MDR-0008)
# ----------------------------------------------------------------------------------------------

P_ATM = 101325.0  # Pa
L_REF = 2.501e6  # J/kg, latent heat of vaporisation of water at 0 degC
L_SLOPE = 2361.0  # J/(kg K), linear decrease of the latent heat with temperature (0-100 degC)


def latent_heat(t_degc: Any) -> Any:
    """Latent heat of vaporisation of water, J/kg (linear fit to the steam tables, 0-100 degC)."""
    return L_REF - L_SLOPE * np.asarray(t_degc, dtype=float)


def saturation_pressure(t_degc: Any) -> Any:
    """Saturation vapour pressure over liquid water, Pa (Magnus form, Alduchov and Eskridge 1996)."""
    t = np.asarray(t_degc, dtype=float)
    return 610.94 * np.exp(17.625 * t / (t + 243.04))


def humidity_ratio(p_vapour: Any, p_total: float = P_ATM) -> Any:
    """Humidity ratio (kg water per kg dry air) of moist air with vapour partial pressure ``p_vapour``."""
    pv = np.asarray(p_vapour, dtype=float)
    return 0.622 * pv / (p_total - pv)


def vapour_pressure(w: Any, p_total: float = P_ATM) -> Any:
    """Vapour partial pressure, Pa, of moist air with humidity ratio ``w``."""
    ww = np.asarray(w, dtype=float)
    return ww * p_total / (0.622 + ww)


def relative_humidity(t_degc: Any, w: Any) -> Any:
    """Relative humidity (0-1) of air at dry-bulb temperature ``t_degc`` with humidity ratio ``w``."""
    return vapour_pressure(w) / saturation_pressure(t_degc)


def wet_bulb_temperature(t_degc: float, w: float) -> float:
    """Thermodynamic wet-bulb (adiabatic-saturation) temperature, degC.

    ASHRAE Handbook - Fundamentals (2017), chapter 1, eq. (33), temperatures in degC and energies in kJ/kg:
    w = ((2501 - 2.326 t*) w_s(t*) - 1.006 (t - t*)) / (2501 + 1.86 t - 4.186 t*). Saturated air returns t.
    """
    from scipy.optimize import brentq

    t = float(t_degc)
    target = float(w)

    def residual(tw: float) -> float:
        ws = float(humidity_ratio(saturation_pressure(tw)))
        return ((2501.0 - 2.326 * tw) * ws - 1.006 * (t - tw)) / (2501.0 + 1.86 * t - 4.186 * tw) - target

    if residual(t) <= 0.0:
        return t
    return float(brentq(residual, -40.0, t, xtol=1e-10, maxiter=200))


def dry_solid_density(props: Properties, c: Any) -> Any:
    """Dry-matter density rho_s = rho(C)/(1 + C), kg/m^3 (bulk density per unit of dry-basis moisture)."""
    cc = np.asarray(c, dtype=float)
    return props.rho(cc) / (1.0 + cc)


@dataclass(frozen=True)
class FluxCap:
    """Upper bound on the surface evaporation mass flux, kg/(m^2 s), sampled on the chamber history.

    ``value(t)`` interpolates linearly between the samples and holds ``plateau`` afterwards; ``rho_s`` converts
    the bound into the (kg/kg) m/s units of the moisture boundary condition; ``scale`` multiplies the bound.
    """

    t: np.ndarray
    cap: np.ndarray
    plateau: float
    rho_s: float
    scale: float = 1.0

    @classmethod
    def wet_bulb(cls, chamber: Chamber, h: float, rho_s: float, scale: float = 1.0) -> FluxCap:
        """Energy limit h (T_air - T_wb)/L_v(T_wb) of a surface held at the wet-bulb temperature."""

        def bound(temp: float, hum: float) -> float:
            tw = wet_bulb_temperature(temp, hum)
            return float(h * (temp - tw) / latent_heat(tw))

        caps = np.array([bound(float(a), float(b)) for a, b in zip(chamber.temp, chamber.hum_scale * chamber.hum)])
        plateau = bound(chamber.temp_plateau * chamber.temp_scale, chamber.hum_plateau * chamber.hum_scale)
        return cls(np.asarray(chamber.t, dtype=float), caps, plateau, rho_s, scale)

    def value(self, t: Any) -> Any:
        return self.scale * np.interp(t, self.t, self.cap, right=self.plateau)

    def describe(self) -> dict[str, Any]:
        return {
            "cap_initial": float(self.cap[0]) * self.scale,
            "cap_plateau": self.plateau * self.scale,
            "rho_s": self.rho_s,
            "scale": self.scale,
        }
