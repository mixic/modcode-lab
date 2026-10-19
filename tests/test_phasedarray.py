# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Tests for linear-array steering, cyclic quantisation, and synthesis."""

import numpy as np
import pytest

from modcode.phasedarray import (
    LinearArray,
    grating_lobe_limit,
    quantize_phases,
    squint_angle,
)


@pytest.fixture()
def array() -> LinearArray:
    return LinearArray(elements=16, spacing_over_lambda=0.5, frequency=10.0e9)


def test_geometry_matches_a_half_wavelength_design(array: LinearArray) -> None:
    assert np.isclose(array.wavelength, 0.0299792458)
    assert np.isclose(array.spacing, array.wavelength / 2.0)
    assert np.isclose(array.aperture, (array.elements - 1) * array.spacing)
    positions = array.positions()
    assert np.isclose(positions[0], -positions[-1])


@pytest.mark.parametrize("theta", [0.0, 15.0, 30.0, 45.0])
def test_exact_ramp_steers_to_the_commanded_angle(
    array: LinearArray, theta: float
) -> None:
    phases = array.exact_phases(theta)
    assert abs(array.peak_angle(phases) - theta) < 0.01
    assert array.array_factor(theta, phases)[0] == pytest.approx(1.0, abs=1e-9)


@pytest.mark.parametrize("theta", [17.0, 30.0])
def test_full_field_sum_agrees_with_the_array_factor(
    array: LinearArray, theta: float
) -> None:
    phases = array.exact_phases(theta)
    assert abs(array.far_field_peak(phases) - array.peak_angle(phases)) < 0.03


def test_element_rounding_keeps_pointing_but_costs_gain(array: LinearArray) -> None:
    theta = 17.0
    phases = array.element_quantized_phases(theta, 8)
    assert abs(array.peak_angle(phases) - theta) < 0.5
    assert -0.5 < array.scan_loss_db(theta, phases) < 0.0


def test_step_rounding_keeps_gain_but_moves_the_beam(array: LinearArray) -> None:
    theta = 17.0
    phases = array.step_quantized_phases(theta, 4)
    assert abs(array.peak_angle(phases) - 30.0) < 0.05
    assert array.scan_loss_db(30.0, phases) == pytest.approx(0.0, abs=1e-9)
    assert array.scan_loss_db(theta, phases) < -10.0


@pytest.mark.parametrize("order", [4, 8, 16])
def test_marechal_strehl_predicts_the_aberration_loss(
    array: LinearArray, order: int
) -> None:
    theta = 17.0
    phases = array.element_quantized_phases(theta, order)
    predicted = 10.0 * np.log10(array.strehl_estimate(phases, theta))
    measured = array.scan_loss_db(theta, phases)
    assert abs(predicted - measured) < 0.05


def test_aberration_rms_matches_a_uniform_quantisation_error(
    array: LinearArray,
) -> None:
    quantum = 2.0 * np.pi / 8.0
    phases = array.element_quantized_phases(17.0, 8)
    rms = array.aberration_rms(phases, 17.0)
    assert rms == pytest.approx(quantum / np.sqrt(12.0), rel=0.06)


def test_quantized_phases_wrap_into_the_hardware_range() -> None:
    values = quantize_phases([0.1, 1.0, 3.0], 4, wrap=True)
    assert np.all(values >= 0.0)
    assert np.all(values < 2.0 * np.pi)
    steps = values / (np.pi / 2.0)
    assert np.allclose(steps, np.round(steps))


def test_grating_lobe_limit_and_full_wavelength_images(array: LinearArray) -> None:
    assert grating_lobe_limit(30.0) == pytest.approx(2.0 / 3.0)
    assert grating_lobe_limit(0.0) == pytest.approx(1.0)

    wide = LinearArray(elements=16, spacing_over_lambda=1.0, frequency=10.0e9)
    levels = wide.array_factor([-30.0, 30.0], wide.exact_phases(30.0))
    assert np.allclose(levels, 1.0)


def test_squint_moves_toward_broadside_with_frequency() -> None:
    assert squint_angle(30.0, 1.0) == pytest.approx(30.0)
    narrowed = squint_angle(30.0, 1.1)
    widened = squint_angle(30.0, 0.9)
    assert narrowed < 30.0 < widened
    assert narrowed == pytest.approx(np.degrees(np.arcsin(0.5 / 1.1)))
    assert squint_angle(30.0, 0.4) == pytest.approx(90.0)


def test_half_power_width_matches_the_scanned_array_formula(
    array: LinearArray,
) -> None:
    big = LinearArray(elements=100, frequency=10.0e9)
    width = big.half_power_width_deg(30.0, big.exact_phases(30.0), span_deg=5.0)
    theory = np.degrees(0.886 / (big.elements * big.spacing_over_lambda))
    theory /= np.cos(np.radians(30.0))
    assert width == pytest.approx(theory, rel=0.02)


def test_invalid_inputs_are_rejected(array: LinearArray) -> None:
    with pytest.raises(ValueError):
        LinearArray(elements=1)
    with pytest.raises(ValueError):
        LinearArray(spacing_over_lambda=0.0)
    with pytest.raises(ValueError):
        array.array_factor(0.0, np.zeros(array.elements - 1))
    with pytest.raises(ValueError):
        quantize_phases([0.0], 0)
