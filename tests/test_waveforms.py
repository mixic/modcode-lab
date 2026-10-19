# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Tests for PSK and Zadoff–Chu waveform generation."""

import numpy as np
import pytest

from modcode.waveforms import (
    cyclic_cross_correlation,
    periodic_autocorrelation,
    psk_constellation,
    zadoff_chu,
)


@pytest.mark.parametrize("order", [2, 4, 8])
def test_psk_points_form_a_cyclic_group(order: int) -> None:
    points = psk_constellation(order)
    assert points.shape == (order,)
    assert np.allclose(np.abs(points), 1.0)
    for first in range(order):
        for second in range(order):
            product = points[first] * points[second]
            assert np.isclose(product, points[(first + second) % order])


def test_psk_phase_offset_is_a_coset() -> None:
    order = 4
    offset = np.pi / 4.0
    rotated = psk_constellation(order, phase_offset=offset)
    assert np.allclose(rotated, np.exp(1j * offset) * psk_constellation(order))

    step = np.exp(2j * np.pi / order)
    assert np.allclose(np.roll(rotated, -1), rotated * step)

    # A coset is not closed under multiplication.
    assert not np.any(np.isclose(rotated[0] * rotated[1], rotated))


@pytest.mark.parametrize(("length", "root"), [(31, 3), (13, 5), (9, 2)])
def test_zc_has_constant_amplitude_and_perfect_autocorrelation(
    length: int, root: int
) -> None:
    sequence = zadoff_chu(length, root)
    assert np.allclose(np.abs(sequence), 1.0)

    autocorrelation = periodic_autocorrelation(sequence)
    assert autocorrelation[0] == pytest.approx(1.0)
    assert np.max(np.abs(autocorrelation[1:])) < 1e-10


def test_zc_cross_correlation_is_one_over_sqrt_length() -> None:
    first = zadoff_chu(31, 1)
    second = zadoff_chu(31, 2)
    cross = cyclic_cross_correlation(first, second)
    assert np.max(np.abs(cross)) == pytest.approx(1.0 / np.sqrt(31.0), abs=1e-10)


def test_composite_length_cross_correlation_needs_coprime_difference() -> None:
    first = zadoff_chu(15, 1)
    coprime = cyclic_cross_correlation(first, zadoff_chu(15, 2))
    non_coprime = cyclic_cross_correlation(first, zadoff_chu(15, 4))

    assert np.max(np.abs(coprime)) == pytest.approx(1.0 / np.sqrt(15.0))
    assert np.max(np.abs(non_coprime)) > 1.0 / np.sqrt(15.0)


def test_periodic_autocorrelation_matches_shaping_helper() -> None:
    from modcode.shaping import periodic_autocorrelation as shaping_autocorrelation

    sequence = zadoff_chu(31, 3)
    assert np.allclose(
        periodic_autocorrelation(sequence),
        shaping_autocorrelation(sequence),
    )


def test_invalid_waveform_inputs() -> None:
    with pytest.raises(ValueError):
        psk_constellation(1)
    with pytest.raises(ValueError):
        psk_constellation(4, phase_offset=np.inf)

    with pytest.raises(ValueError):
        zadoff_chu(10, 3)
    with pytest.raises(ValueError):
        zadoff_chu(15, 3)
    with pytest.raises(ValueError):
        zadoff_chu(15, 15)

    with pytest.raises(ValueError):
        periodic_autocorrelation([])
    with pytest.raises(ValueError):
        cyclic_cross_correlation([1.0, 0.0], [1.0])
