# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Tests for prescribed residual-index densities and finite counts."""

from fractions import Fraction

import numpy as np
import pytest

from modcode.number_theory.residual_index import (
    artin_constant,
    count_residual_indices,
    known_decimal_density_multiple,
    residual_index_convergence,
    residual_index_density,
)


def test_artin_constant_matches_the_known_value() -> None:
    assert artin_constant(10000) == pytest.approx(0.3739558136, abs=5e-7)


def test_known_decimal_density_multiples() -> None:
    assert known_decimal_density_multiple(1) == Fraction(1, 1)
    assert known_decimal_density_multiple(2) == Fraction(3, 4)
    assert known_decimal_density_multiple(4) == Fraction(29, 152)
    assert known_decimal_density_multiple(3) is None


def test_density_matches_the_closed_forms() -> None:
    constant = artin_constant(1000)
    assert residual_index_density(1, terms=80_000) == pytest.approx(
        constant, rel=2e-4
    )
    assert residual_index_density(2, terms=80_000) == pytest.approx(
        0.75 * constant, rel=2e-4
    )
    assert residual_index_density(4, terms=80_000) == pytest.approx(
        Fraction(29, 152) * constant, rel=3e-4
    )
    # The often-quoted 29/52 value is contradicted by the formula and counts.
    assert residual_index_density(4, terms=80_000) < 0.1


def test_finite_counts_are_close_to_the_predicted_densities() -> None:
    counts = count_residual_indices(200_000)
    assert counts.total > 0
    assert counts.ratio(1) == pytest.approx(0.373956, abs=0.01)
    assert counts.ratio(2) == pytest.approx(0.280467, abs=0.01)
    assert counts.ratio(4) == pytest.approx(29 / 152 * 0.373956, abs=0.01)


def test_convergence_series_have_the_expected_shape() -> None:
    series = residual_index_convergence([10_000, 50_000])
    assert set(series) == {1, 2, 4}
    for values in series.values():
        assert values.shape == (2,)
        assert np.all((values >= 0.0) & (values <= 1.0))


def test_invalid_inputs_are_rejected() -> None:
    with pytest.raises(ValueError):
        artin_constant(0)
    with pytest.raises(ValueError):
        residual_index_density(0)
    with pytest.raises(ValueError):
        residual_index_density(2, base=2)
    with pytest.raises(ValueError):
        count_residual_indices(2)
    with pytest.raises(ValueError):
        count_residual_indices(100, base=2)
    with pytest.raises(ValueError):
        residual_index_convergence([])
