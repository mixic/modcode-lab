# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Tests separating distribution shaping from sequence symmetry properties."""

import numpy as np
import pytest

from modcode.shaping import (
    maxwell_boltzmann_probabilities,
    periodic_autocorrelation,
    sample_shaped_symbols,
)


def test_energy_shaping_favors_inner_points_and_has_unit_probability() -> None:
    constellation = np.array([-3, -1, 1, 3], dtype=np.complex128)
    probabilities = maxwell_boltzmann_probabilities(constellation, 0.5)
    assert np.isclose(probabilities.sum(), 1.0)
    assert probabilities[1] > probabilities[0]
    assert probabilities[2] > probabilities[3]
    assert np.allclose(
        maxwell_boltzmann_probabilities(constellation, 0.0),
        np.full(4, 0.25),
    )


def test_shaped_sampler_is_seed_reproducible_and_matches_shape() -> None:
    points = np.array([-1.0, 1.0], dtype=np.complex128)
    probabilities = np.array([0.8, 0.2])
    first = sample_shaped_symbols(
        points, probabilities, 100_000, np.random.default_rng(7)
    )
    second = sample_shaped_symbols(
        points, probabilities, 100_000, np.random.default_rng(7)
    )
    assert np.array_equal(first, second)
    assert np.isclose(np.mean(first == points[0]), 0.8, atol=0.01)


def test_periodic_autocorrelation_is_invariant_under_cyclic_shift() -> None:
    sequence = np.array([1, 1j, -1, -1j, 1], dtype=np.complex128)
    assert np.allclose(
        periodic_autocorrelation(sequence),
        periodic_autocorrelation(np.roll(sequence, 2)),
    )


def test_arbitrary_permutation_preserves_histogram_but_not_autocorrelation() -> None:
    sequence = np.array([0, 0, 1, 1, 1, 0], dtype=np.complex128)
    permutation = np.array([0, 2, 4, 1, 3, 5])
    permuted = sequence[permutation]
    assert np.array_equal(np.sort(sequence), np.sort(permuted))
    assert not np.allclose(
        periodic_autocorrelation(sequence),
        periodic_autocorrelation(permuted),
    )


@pytest.mark.parametrize("parameter", [-1.0, float("nan")])
def test_shaping_parameter_must_be_finite_and_nonnegative(parameter: float) -> None:
    with pytest.raises(ValueError):
        maxwell_boltzmann_probabilities(np.array([1 + 0j, -1 + 0j]), parameter)
