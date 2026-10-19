# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Tests for square Gray-mapped QAM."""

import itertools

import numpy as np
import pytest

from modcode.modulation.qam import QAM


@pytest.mark.parametrize("order", [16, 64])
def test_horizontal_and_vertical_neighbors_are_gray_adjacent(order: int) -> None:
    modulation = QAM(order)
    points = modulation.constellation
    for first, second in itertools.combinations(points, 2):
        distance = abs(points[first] - points[second])
        hamming_distance = sum(a != b for a, b in zip(first, second, strict=True))
        minimum_distance = min(
            abs(point - other)
            for point in points.values()
            for other in points.values()
            if point != other
        )
        if np.isclose(
            distance,
            minimum_distance,
        ):
            assert hamming_distance == 1


@pytest.mark.parametrize("order", [16, 64])
def test_qam_round_trip_and_unit_average_power(order: int) -> None:
    modulation = QAM(order)
    bits = np.random.default_rng(order).integers(
        0, 2, size=modulation.bits_per_symbol * 100, dtype=np.uint8
    )
    symbols = modulation.modulate(bits)
    assert np.array_equal(modulation.demodulate(symbols), bits)
    assert np.isclose(
        np.mean(np.abs(list(modulation.constellation.values())) ** 2), 1.0
    )


@pytest.mark.parametrize("order", [3, 8, 15, 16.0])
def test_rejects_non_square_power_of_two(order: int | float) -> None:
    with pytest.raises(ValueError):
        QAM(order)
