# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Tests for the sphere-packing gain accounting in the companion note."""

import math

import pytest

from modcode.metrics.sphere_packing import (
    LATTICES,
    SPHERE_SHAPING_GAIN_DB,
    LatticeFacts,
    center_density,
    ebn0_at_ber,
    hermite_parameter,
    nominal_coding_gain_db,
    shannon_ebn0_db,
    shaping_gain_db,
)
from modcode.modulation.qam import QAM
from modcode.modulation.qpsk import QPSK


def _lattice(name: str) -> LatticeFacts:
    for lattice in LATTICES:
        if lattice.name == name:
            return lattice
    raise AssertionError(f"unknown lattice {name}")


def test_nominal_coding_gains_match_the_literature():
    """D4, E8 and Leech give the textbook 1.51, 3.01 and 6.02 dB."""
    assert nominal_coding_gain_db(_lattice("Z^2")) == pytest.approx(0.0, abs=1e-12)
    assert nominal_coding_gain_db(_lattice("D4")) == pytest.approx(1.5051, abs=1e-3)
    assert nominal_coding_gain_db(_lattice("E8")) == pytest.approx(3.0103, abs=1e-3)
    assert nominal_coding_gain_db(_lattice("Leech")) == pytest.approx(6.0206, abs=1e-3)


def test_hermite_parameter_is_scale_free():
    """Doubling a lattice's volume halves d^2 / det^(2/n), not the other way."""
    lattice = _lattice("E8")
    assert hermite_parameter(lattice) == pytest.approx(2.0, abs=1e-12)


def test_center_densities_match_the_literature():
    """The packing densities of D4, E8 and Leech are 1/8, 1/16 and 1."""
    assert center_density(_lattice("D4")) == pytest.approx(0.125, abs=1e-12)
    assert center_density(_lattice("E8")) == pytest.approx(0.0625, abs=1e-12)
    assert center_density(_lattice("Leech")) == pytest.approx(1.0, abs=1e-12)


def test_shaping_gain_ceiling_is_pi_e_over_six():
    """The sphere beats the cube by exactly pi*e/6, which is 1.53 dB."""
    assert SPHERE_SHAPING_GAIN_DB == pytest.approx(1.5329, abs=1e-3)
    assert 10.0 ** (SPHERE_SHAPING_GAIN_DB / 10.0) == pytest.approx(
        math.pi * math.e / 6.0, rel=1e-12
    )


def test_shaping_gain_of_the_cube_is_zero():
    """A region identical to the cube gains nothing."""
    assert shaping_gain_db(1.0 / 12.0) == pytest.approx(0.0, abs=1e-12)


def test_shaping_gain_rejects_regions_denser_than_a_sphere():
    """No region can beat the sphere; asking for one is an error."""
    with pytest.raises(ValueError):
        shaping_gain_db(1.0 / (4.0 * math.pi * math.e))


def test_shannon_limit_of_the_binary_input_is_zero_db():
    """eta = 1/2 bit per real axis is the 0 dB point of the classic curve."""
    assert shannon_ebn0_db(0.5) == pytest.approx(0.0, abs=1e-12)
    assert shannon_ebn0_db(1.0) == pytest.approx(1.7609, abs=1e-3)


def test_uncoded_qpsk_gap_to_capacity_is_about_eight_db():
    """Uncoded QPSK needs ~9.6 dB at 1e-5, which is ~7.8 dB above the limit."""
    assert ebn0_at_ber(QPSK(), 1.0e-5) == pytest.approx(9.588, abs=0.01)
    gap = ebn0_at_ber(QPSK(), 1.0e-5) - shannon_ebn0_db(1.0)
    assert gap == pytest.approx(7.827, abs=0.01)


def test_square_qam_gaps_are_all_between_seven_and_eight_db():
    """The gap is a property of uncoded QAM, not of the constellation size."""
    for order, eta in ((16, 2.0), (64, 3.0), (256, 4.0)):
        gap = ebn0_at_ber(QAM(order), 1.0e-5) - shannon_ebn0_db(eta)
        assert 7.4 < gap < 7.9


def test_ebn0_at_ber_rejects_impossible_targets():
    """A target outside (0, 1), or a bracket that already meets it, is an error."""
    with pytest.raises(ValueError):
        ebn0_at_ber(QPSK(), 0.0)
    with pytest.raises(ValueError):
        ebn0_at_ber(QPSK(), 0.6, bracket=(20.0, 40.0))
    with pytest.raises(ValueError):
        shannon_ebn0_db(0.0)
