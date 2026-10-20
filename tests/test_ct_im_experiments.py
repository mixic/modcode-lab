# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Tests for the CT-IM experimental milestones and baselines."""

import numpy as np
import pytest

from modcode.metrics.index_modulation import (
    block_papr_values,
    papr_ccdf,
    simulate_index_modulation,
)
from modcode.modulation.ct_im import (
    CTIM,
    BCHIndexCode,
    ClassicalOFDMIM,
    ConstantWeightCyclicCode,
    ShortenedConstantWeightCode,
)


def test_shortened_codebook_is_constant_weight_with_distance_four() -> None:
    code = ShortenedConstantWeightCode()
    assert (code.n, code.weight, len(code), code.minimum_distance) == (6, 3, 4, 4)
    assert np.all(code.patterns.sum(axis=1) == 3)
    for label in range(len(code)):
        word = code.encode(label)
        for position in range(code.n):
            corrupted = word.copy()
            corrupted[position] ^= 1
            decoded, distance = code.decode(corrupted)
            assert decoded == label
            assert distance == 1


def test_puncturing_without_shortening_breaks_constant_weight() -> None:
    orbit = ConstantWeightCyclicCode()
    punctured = orbit.patterns[:, 1:]
    weights = punctured.sum(axis=1)
    assert set(weights.tolist()) == {2, 3}
    assert not np.all(weights == weights[0])


def test_shortened_ctim_rate_and_budget() -> None:
    ct = CTIM(M=4, code=ShortenedConstantWeightCode(), detector="ml")
    assert ct.n == 6
    assert ct.active == 3
    assert ct.index_bits == 2
    assert ct.bits_per_block == 8
    assert ct.spectral_efficiency == pytest.approx(8 / 6)


def test_bch_index_code_is_the_weight_seven_orbit() -> None:
    code = BCHIndexCode(weight=7)
    assert (code.n, code.weight, len(code), code.minimum_distance) == (15, 7, 15, 8)
    assert np.all(code.patterns.sum(axis=1) == 7)
    for label in range(len(code)):
        for steps in range(code.n):
            assert code.is_codeword(np.roll(code.encode(label), steps))


def test_bch_ctim_block_budget() -> None:
    ct = CTIM(M=4, code=BCHIndexCode(weight=7), detector="ml")
    assert ct.n == 15
    assert ct.index_bits == 4
    assert ct.bits_per_block == 18
    assert ct.spectral_efficiency == pytest.approx(18 / 15)


def test_classical_ofdm_im_codebook_and_rate() -> None:
    scheme = ClassicalOFDMIM()
    assert len(scheme.index_code) == 32
    assert scheme.index_code.n == 7
    assert np.all(scheme.index_code.patterns.sum(axis=1) == 3)
    assert scheme.bits_per_block == 11
    assert scheme.index_bits == 5
    assert scheme.spectral_efficiency == pytest.approx(11 / 7)


def test_classical_ofdm_im_noiseless_round_trip_and_detection() -> None:
    scheme = ClassicalOFDMIM()
    rng = np.random.default_rng(4)
    for _ in range(20):
        label = int(rng.integers(0, 32))
        bits = np.concatenate(
            (
                np.asarray(
                    [(label >> shift) & 1 for shift in range(4, -1, -1)],
                    dtype=np.uint8,
                ),
                rng.integers(0, 2, 6, dtype=np.uint8),
            )
        )
        np.testing.assert_array_equal(scheme.decode(scheme.encode(bits)), bits)


def test_ml_detector_beats_the_threshold_detector() -> None:
    threshold = CTIM(M=4, detector="threshold")
    matched = CTIM(M=4, detector="ml")
    low = simulate_index_modulation(
        threshold, 0.0, blocks=600, rng=np.random.default_rng(7)
    )
    high = simulate_index_modulation(
        matched, 0.0, blocks=600, rng=np.random.default_rng(7)
    )
    assert high.index_ber < 0.6 * low.index_ber


def test_rayleigh_costs_snr_at_the_same_eb_n0() -> None:
    scheme = CTIM(M=4, detector="ml")
    awgn = simulate_index_modulation(
        scheme, 10.0, blocks=400, rng=np.random.default_rng(9)
    )
    fading = simulate_index_modulation(
        scheme, 10.0, blocks=400, rayleigh=True, rng=np.random.default_rng(9)
    )
    assert fading.ber > 10.0 * awgn.ber


def test_cyclic_shift_is_papr_invariant() -> None:
    rng = np.random.default_rng(5)
    scheme = CTIM(M=4, detector="ml")
    rows = []
    for _ in range(64):
        label = int(rng.integers(0, len(scheme.index_code)))
        bits = np.concatenate(
            (
                np.asarray(
                    [(label >> shift) & 1 for shift in range(2, -1, -1)],
                    dtype=np.uint8,
                ),
                rng.integers(0, 2, 6, dtype=np.uint8),
            )
        )
        rows.append(scheme.encode(bits))
    vectors = np.asarray(rows)
    base = block_papr_values(vectors)
    for shift in (1, 3, 6):
        np.testing.assert_allclose(
            block_papr_values(vectors, offset=shift), base, atol=1e-10
        )
    thresholds = np.asarray([3.0, 5.0, 7.0])
    ccdf = papr_ccdf(vectors, thresholds)
    assert np.all(ccdf >= 0.0) and np.all(ccdf <= 1.0)
    assert np.all(np.diff(ccdf) <= 1e-12)


def test_invalid_reserved_labels_and_bit_lengths() -> None:
    with pytest.raises(ValueError, match="reserved"):
        CTIM(M=4).encode(np.array([1, 1, 1, 0, 0, 0, 0, 0, 0], dtype=np.uint8))
    with pytest.raises(ValueError, match="expected"):
        ClassicalOFDMIM().encode(np.zeros(10, dtype=np.uint8))
    with pytest.raises(ValueError, match="detector"):
        CTIM(M=4, detector="mlr")
    with pytest.raises(ValueError, match="weight"):
        BCHIndexCode(weight=1)
