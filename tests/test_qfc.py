# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Tests for photon budgets, transduction noise, and transport loss."""

import numpy as np
import pytest

from modcode.qfc import (
    SPEED_OF_LIGHT,
    Transducer,
    covariance_squeezing_db,
    fiber_loss,
    free_space_link_loss,
    photon_energy,
    squeezed_vacuum_covariance,
    thermal_occupation,
    transport_budget,
    wigner_function,
)

RF = 10.0e9
OPTICAL = SPEED_OF_LIGHT / 1550e-9


def test_photon_energy_ratio_is_about_twenty_thousand() -> None:
    ratio = photon_energy(OPTICAL) / photon_energy(RF)
    assert 1.9e4 < ratio < 2.0e4


def test_thermal_occupation_matches_the_bose_einstein_formula() -> None:
    assert thermal_occupation(RF, 290.0) == pytest.approx(603.76, rel=1e-3)
    assert thermal_occupation(RF, 4.0) == pytest.approx(7.8446, rel=1e-3)
    assert thermal_occupation(OPTICAL, 290.0) < 1e-13
    assert thermal_occupation(RF, 0.0) == 0.0


def test_calibrated_and_uncalibrated_excess_noise_differ() -> None:
    device = Transducer(0.5, 0.1)
    assert device.excess_noise(calibrated=True) == pytest.approx(0.2)
    assert device.excess_noise(calibrated=False) == pytest.approx(1.2)
    assert Transducer(1.0, 0.01).excess_noise(calibrated=False) == pytest.approx(
        0.01
    )


def test_loss_only_conversion_preserves_a_vacuum_covariance() -> None:
    assert np.allclose(Transducer(0.5).convert_covariance(np.eye(2)), np.eye(2))


def test_squeezing_is_degraded_by_efficiency_and_added_noise() -> None:
    covariance = squeezed_vacuum_covariance(10.0)
    assert covariance_squeezing_db(covariance) == pytest.approx(-10.0)
    at_half = Transducer(0.5).convert_covariance(covariance)
    assert covariance_squeezing_db(at_half) == pytest.approx(
        10.0 * np.log10(0.55), abs=0.01
    )
    noisy = Transducer(0.5, 0.1).surviving_squeezing_db(10.0)
    assert noisy > covariance_squeezing_db(at_half)


def test_a_four_kelvin_bath_erases_ten_db_of_squeezing() -> None:
    bath = thermal_occupation(RF, 4.0)
    assert Transducer(1.0).surviving_squeezing_db(10.0, bath_occupation=bath) > 0.0
    cold = Transducer(1.0).surviving_squeezing_db(10.0)
    assert cold == pytest.approx(-10.0)


def test_free_space_link_loss_matches_the_friis_arithmetic() -> None:
    rf = free_space_link_loss(100e3, SPEED_OF_LIGHT / RF, 0.3, 0.3)
    optical = free_space_link_loss(100e3, 1550e-9, 0.3, 0.3)
    assert rf == pytest.approx(92.55, abs=0.1)
    assert optical == pytest.approx(6.82, abs=0.1)
    assert fiber_loss(100e3) == pytest.approx(22.0)


def test_pointing_jitter_costs_optical_links_dearly() -> None:
    clean = free_space_link_loss(100e3, 1550e-9, 0.3, 0.3)
    jittered = free_space_link_loss(
        100e3, 1550e-9, 0.3, 0.3, pointing_jitter_rad=10e-6
    )
    assert jittered - clean > 20.0


def test_transport_budget_has_one_row_per_carrier() -> None:
    paths = transport_budget(100e3)
    assert len(paths) == 6
    losses = {path.label: path.loss_db for path in paths}
    assert losses["Fiber 1550 nm (0.2 dB/km)"] == pytest.approx(22.0)


def test_wigner_function_is_normalized() -> None:
    grid = np.linspace(-6.0, 6.0, 601)
    spacing = grid[1] - grid[0]
    wigner = wigner_function(grid, grid, np.eye(2))
    assert wigner.max() == pytest.approx(1.0 / (2.0 * np.pi))
    assert float(np.sum(wigner) * spacing**2) == pytest.approx(1.0, abs=0.01)


def test_invalid_transducer_and_covariance_inputs() -> None:
    with pytest.raises(ValueError):
        Transducer(0.0)
    with pytest.raises(ValueError):
        Transducer(1.1)
    with pytest.raises(ValueError):
        Transducer(0.5, -0.1)
    with pytest.raises(ValueError):
        covariance_squeezing_db(np.zeros((2, 2)))
    with pytest.raises(ValueError):
        wigner_function([0.0], [0.0], np.eye(3))
