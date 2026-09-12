"""Unit tests of the helpers behind the derived paper quantities (MDR-0009)."""

from __future__ import annotations

import numpy as np

from pipelines.a.derived import first_time_within, fourier_number, richardson_error


def test_fourier_number() -> None:
    assert fourier_number(2.0, 3.0, 2.0) == 1.5


def test_richardson_error_recovers_the_limit_of_a_second_order_sequence() -> None:
    limit = 10.0
    values = [limit + 4.0 / 4**k for k in range(4)]  # error /4 per level: grids halved, order 2
    est, err = richardson_error(values, order=2.0)
    assert abs(est - limit) < 1e-12
    assert abs(err - 4.0 / 4**2) < 1e-12


def test_first_time_within_requires_the_gap_to_stay_small() -> None:
    t = np.arange(6.0)
    a = np.array([0.0, 0.9, 0.2, 0.7, 0.1, 0.05])
    b = np.zeros(6)
    assert first_time_within(t, a, b, 0.5) == 4.0
    assert first_time_within(t, a + 10.0, b, 0.5) is None
    assert first_time_within(t, b, b, 0.5) == 0.0


def test_value_at_picks_the_nearest_recorded_instant() -> None:
    from pipelines.a.derived import value_at

    t = np.array([0.0, 60.0, 120.0, 180.0])
    v = np.array([1.0, 2.0, 3.0, 4.0])
    assert value_at(t, v, 120.0) == 3.0
    assert value_at(t, v, 119.0) == 3.0
    assert value_at(t, v, 1000.0) == 4.0
