# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Tests for the short BI-AWGN Polar code and its SC/SCL decoder."""

import numpy as np
import pytest

from modcode.channel.awgn import AWGNChannel
from modcode.fec.crc import CRC8
from modcode.fec.polar import (
    PolarAWGN,
    awgn_synthetic_llr_means,
    gaussian_phi,
    gaussian_phi_inverse,
)
from modcode.modulation.bpsk import BPSK
from modcode.pipeline.chain import Link, run_link


def test_gaussian_phi_matches_direct_integration() -> None:
    """phi(x) = 1 - E[tanh(U/2)] for U ~ N(x, 2x), checked by quadrature."""
    for x in (0.5, 1.0, 2.0, 5.0):
        spread = np.sqrt(2.0 * x)
        nodes = x + spread * np.linspace(-8.0, 8.0, 20001)
        density = np.exp(-0.5 * ((nodes - x) / spread) ** 2)
        density /= density.sum()
        expected = 1.0 - np.sum(np.tanh(0.5 * nodes) * density)
        assert gaussian_phi(x) == pytest.approx(expected, rel=2e-3)


def test_gaussian_phi_inverse_roundtrip() -> None:
    for x in (1e-3, 1e-2, 0.1, 1.0, 10.0, 100.0):
        assert gaussian_phi_inverse(gaussian_phi(x)) == pytest.approx(x, rel=2e-2)


def test_construction_matches_the_bec_reliability_order_at_eight() -> None:
    """The GA ordering for N=8 reproduces the exact BEC order at erasure 1/2."""
    means = awgn_synthetic_llr_means(8, 2.0, 0.5)
    assert means.size == 8
    assert np.all(means > 0.0)
    assert np.argsort(means).tolist() == [0, 1, 2, 4, 3, 5, 6, 7]
    code = PolarAWGN(8, 4, crc=None, list_size=1)
    assert code.information_indices.tolist() == [3, 5, 6, 7]
    assert code.frozen_mask.tolist() == [
        True,
        True,
        True,
        False,
        True,
        False,
        False,
        False,
    ]


def test_noiseless_round_trip_for_sc_and_scl() -> None:
    rng = np.random.default_rng(11)
    for crc in (None, CRC8):
        for list_size in (1, 4):
            code = PolarAWGN(64, 24, crc=crc, list_size=list_size)
            data = rng.integers(0, 2, 24 * 3, dtype=np.uint8)
            encoded = code.encode(data)
            assert encoded.size == 64 * 3
            np.testing.assert_array_equal(code.decode(encoded), data)
            llr = (1.0 - 2.0 * encoded.astype(np.float64)) * 6.0
            np.testing.assert_array_equal(code.decode_soft(llr), data)


def _reference_sc(llr: np.ndarray, frozen: np.ndarray) -> np.ndarray:
    """Independent recursive SC decoder, written directly from the recursion."""
    if llr.size == 1:
        bit = 0 if frozen[0] else int(llr[0] < 0.0)
        return np.asarray([bit], dtype=np.uint8)
    half = llr.size // 2
    first, second = llr[:half], llr[half:]
    products = np.tanh(0.5 * first) * np.tanh(0.5 * second)
    merged = 2.0 * np.arctanh(np.clip(products, -1.0 + 1e-12, 1.0 - 1e-12))
    left = _reference_sc(merged, frozen[:half])
    left_beta = _reference_transform(left)
    right = _reference_sc(second + (1.0 - 2.0 * left_beta) * first, frozen[half:])
    return np.concatenate((left, right))


def _reference_transform(bits: np.ndarray) -> np.ndarray:
    if bits.size == 1:
        return bits.copy()
    half = bits.size // 2
    left = _reference_transform(bits[:half])
    right = _reference_transform(bits[half:])
    return np.concatenate((left ^ right, right))


def test_sc_matches_an_independent_reference_decoder() -> None:
    for block_length in (8, 16, 32):
        code = PolarAWGN(block_length, block_length // 2, crc=None, list_size=1)
        rng = np.random.default_rng(block_length)
        for _ in range(5):
            llr = rng.normal(2.0, 6.0, block_length)
            reference = _reference_sc(llr, code.frozen_mask)
            np.testing.assert_array_equal(
                code.decode_soft(llr), reference[code.information_indices]
            )


def test_work_per_block_is_n_log_n_for_sc_and_scales_with_the_list() -> None:
    code = PolarAWGN(64, 32, crc=None, list_size=1)
    rng = np.random.default_rng(5)
    noisy = rng.normal(0.0, 1.0, 64)
    code.decode_soft(noisy)
    assert code.last_work == 64 * 6

    larger = PolarAWGN(64, 32, crc=None, list_size=4)
    larger.decode_soft(noisy)
    assert 2.5 * 64 * 6 <= larger.last_work <= 4.2 * 64 * 6


def test_list_decoding_reduces_errors_on_a_fixed_sample() -> None:
    rng = np.random.default_rng(2027)
    data = rng.integers(0, 2, 39_936, dtype=np.uint8)
    sc = PolarAWGN(256, 128, crc=None, list_size=1)
    scl = PolarAWGN(256, 128, crc=None, list_size=4)
    llr = (1.0 - 2.0 * sc.encode(data).astype(np.float64)) * 2.0
    noise = rng.normal(0.0, 1.5, llr.size)
    received = llr + noise
    sc_errors = int(np.count_nonzero(sc.decode_soft(received) != data))
    scl_errors = int(np.count_nonzero(scl.decode_soft(received) != data))
    assert scl_errors < 0.5 * sc_errors


def test_decisions_are_invariant_to_a_positive_llr_scale() -> None:
    code = PolarAWGN(64, 32, crc=None, list_size=4)
    rng = np.random.default_rng(12)
    llr = rng.normal(0.0, 4.0, 64)
    np.testing.assert_array_equal(code.decode_soft(llr), code.decode_soft(3.7 * llr))


def test_link_soft_decoding_at_three_db() -> None:
    result = run_link(
        Link(BPSK(), AWGNChannel(), PolarAWGN(256, 128, list_size=4), decoding="soft"),
        n_bits=20_000,
        ebn0_db=3.0,
        rng=np.random.default_rng(9),
    )
    assert result.decoding == "soft"
    assert result.ber < 1e-3


def test_invalid_parameters_and_lengths() -> None:
    with pytest.raises(ValueError, match="power of two"):
        PolarAWGN(60, 30)
    with pytest.raises(ValueError, match="data_bits plus the CRC"):
        PolarAWGN(16, 16, crc=CRC8)
    with pytest.raises(ValueError, match="data_bits must be positive"):
        PolarAWGN(16, 0, crc=None)
    with pytest.raises(ValueError, match="list_size must be positive"):
        PolarAWGN(16, 8, list_size=0)
    code = PolarAWGN(16, 8, crc=None, list_size=2)
    with pytest.raises(ValueError, match="multiple"):
        code.encode(np.zeros(7, dtype=np.uint8))
    with pytest.raises(ValueError, match="multiple"):
        code.decode(np.zeros(15, dtype=np.uint8))
    with pytest.raises(ValueError, match="multiple"):
        code.decode_soft(np.zeros(15))
    with pytest.raises(ValueError, match="finite"):
        code.decode_soft(np.full(16, np.inf))
    with pytest.raises(ValueError, match="real"):
        code.decode_soft(np.ones(16, dtype=np.complex128))
