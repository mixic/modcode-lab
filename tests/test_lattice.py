# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Tests for D4/E8 spherical-cut lattice modulation."""

import numpy as np
import pytest

from modcode.metrics.lattice import (
    qam_per_dimension_points,
    simulate_lattice_ber,
)
from modcode.modulation.lattice import LatticeModulation
from modcode.modulation.qam import QAM
from modcode.modulation.qpsk import QPSK


def _min_distance(points: np.ndarray) -> float:
    distances = np.sum(
        (points[:, None] - points[None, :]) ** 2, axis=2
    )
    np.fill_diagonal(distances, np.inf)
    return float(np.sqrt(np.min(distances)))


def test_d4_and_e8_have_unit_average_symbol_energy() -> None:
    for kind, m_bits in (("D4", 4), ("D4", 8), ("E8", 8)):
        modulation = LatticeModulation(kind, m_bits)
        assert modulation.average_energy() == pytest.approx(1.0)
        assert modulation.points.shape == (modulation.M, modulation.K)
        assert modulation.constellation.shape == (
            modulation.M,
            modulation.K // 2,
        )


def test_d4_and_e8_round_trip() -> None:
    rng = np.random.default_rng(3)
    for kind, m_bits in (("D4", 4), ("D4", 8), ("E8", 8)):
        modulation = LatticeModulation(kind, m_bits)
        bits = rng.integers(
            0, 2, modulation.bits_per_symbol * 100, dtype=np.uint8
        )
        symbols = modulation.modulate(bits)
        assert symbols.shape == (100, modulation.K // 2)
        assert np.array_equal(modulation.demodulate(symbols), bits)


def test_gray_labeler_is_a_permutation() -> None:
    for labeler in ("lex", "gray"):
        modulation = LatticeModulation("D4", 4, labeler=labeler)
        patterns = {tuple(row) for row in modulation.bit_patterns}
        assert len(patterns) == modulation.M


def test_d4_16_matches_qpsk_min_distance_at_equal_dimension() -> None:
    d4 = LatticeModulation("D4", 4)
    qpsk_points, _ = qam_per_dimension_points(QPSK())
    assert d4.min_distance(1.0) == pytest.approx(_min_distance(qpsk_points))


def test_d4_256_beats_16qam_min_distance_at_equal_dimension() -> None:
    d4 = LatticeModulation("D4", 8)
    qam_points, _ = qam_per_dimension_points(QAM(16))
    assert d4.min_distance(1.0) > _min_distance(qam_points)


def test_e8_256_beats_qpsk_min_distance_at_equal_dimension() -> None:
    e8 = LatticeModulation("E8", 8)
    qpsk_points, _ = qam_per_dimension_points(QPSK())
    assert e8.min_distance(1.0) > _min_distance(qpsk_points)


def test_lattice_ber_is_a_fraction() -> None:
    modulation = LatticeModulation("D4", 4)
    ber = simulate_lattice_ber(
        modulation, 10.0, 20_000, np.random.default_rng(10)
    )
    assert 0.0 <= ber <= 1.0
    assert ber < 0.01


def test_lattice_rejects_invalid_inputs() -> None:
    with pytest.raises(ValueError, match="kind"):
        LatticeModulation("D3", 4)
    with pytest.raises(ValueError, match="m_bits"):
        LatticeModulation("D4", 0)
    modulation = LatticeModulation("D4", 4)
    with pytest.raises(ValueError, match="multiple"):
        modulation.modulate(np.zeros(3, dtype=np.uint8))
    with pytest.raises(ValueError, match="shape"):
        modulation.demodulate(np.zeros((2, 3), dtype=complex))
