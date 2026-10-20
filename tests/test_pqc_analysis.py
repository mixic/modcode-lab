# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Tests for the QC-MDPC security experiments."""

from math import comb, log2

import numpy as np
import pytest

from modcode.pqc import (
    QCMDPC,
    decoder_success_rates,
    dual_weight_distribution,
    syndrome_search_space_bits,
)


def test_search_space_matches_the_binomial_coefficient() -> None:
    assert syndrome_search_space_bits(10, 3) == pytest.approx(
        log2(comb(10, 3))
    )
    assert syndrome_search_space_bits(94, 6) == pytest.approx(
        log2(comb(94, 6))
    )


def test_search_space_grows_with_error_weight() -> None:
    values = [syndrome_search_space_bits(94, weight) for weight in range(1, 7)]
    assert values[0] < values[-1]
    assert all(
        first < second
        for first, second in zip(values, values[1:], strict=False)
    )


def test_dual_weight_distribution_contains_the_sparse_secret_rows() -> None:
    code = QCMDPC(block_length=15, row_weight=3, seed=3)
    counts = dual_weight_distribution(code)
    assert counts.sum() == 1 << 15
    minimum = int(np.flatnonzero(counts[1:])[0] + 1)
    assert minimum <= 2 * code.row_weight
    assert counts[2 * code.row_weight] >= code.block_length


def test_decoder_success_rates_are_probabilities() -> None:
    code = QCMDPC(block_length=31, row_weight=3, seed=5)
    weights, rates = decoder_success_rates(
        code, [1, 2, 3], trials=4, seed=7
    )
    assert np.array_equal(weights, [1, 2, 3])
    assert rates.shape == (3,)
    assert np.all((rates >= 0.0) & (rates <= 1.0))
    assert rates[0] >= 0.5


def test_invalid_analysis_inputs() -> None:
    with pytest.raises(ValueError):
        syndrome_search_space_bits(10, 11)
    with pytest.raises(ValueError):
        syndrome_search_space_bits(0, 0)
    with pytest.raises(ValueError):
        dual_weight_distribution(
            QCMDPC(block_length=31, row_weight=3, seed=5)
        )
    with pytest.raises(ValueError):
        decoder_success_rates(
            QCMDPC(block_length=15, row_weight=3, seed=3),
            [1],
            trials=0,
        )
