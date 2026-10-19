# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Tests for APSK ring geometry and cyclic phase actions."""

import numpy as np
import pytest

from modcode.modulation.apsk import APSK


def test_16_apsk_has_unit_average_energy_and_round_trips() -> None:
    modulation = APSK((4, 12), (1.0, 2.0))
    constellation = modulation.constellation
    symbols = np.asarray(list(constellation.values()))
    bits = np.asarray(list(constellation), dtype=np.uint8).reshape(-1)
    assert len(constellation) == 16
    assert np.isclose(np.mean(np.abs(symbols) ** 2), 1.0)
    assert np.array_equal(modulation.demodulate(modulation.modulate(bits)), bits)


@pytest.mark.parametrize(("ring_index", "ring_order"), [(0, 4), (1, 12)])
def test_each_ring_is_an_orbit_under_cyclic_phase_steps(
    ring_index: int, ring_order: int
) -> None:
    modulation = APSK((4, 12), (1.0, 2.0))
    points = modulation.ring_points(ring_index)
    for shift in range(ring_order):
        expected = points * np.exp(2j * np.pi * shift / ring_order)
        assert np.allclose(modulation.rotate_ring(ring_index, shift), expected)
    assert np.allclose(modulation.rotate_ring(ring_index, ring_order), points)


def test_apsk_rejects_invalid_ring_definitions() -> None:
    with pytest.raises(ValueError):
        APSK((4, 10), (1.0, 2.0))
    with pytest.raises(ValueError):
        APSK((4, 12), (2.0, 1.0))
