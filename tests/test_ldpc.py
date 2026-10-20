# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Tests for the short quasi-cyclic LDPC code and min-sum decoder."""

import numpy as np
import pytest

from modcode.channel.awgn import AWGNChannel
from modcode.fec.ldpc import LDPC36
from modcode.modulation.bpsk import BPSK
from modcode.pipeline.chain import Link, run_link


@pytest.fixture(scope="module")
def code() -> LDPC36:
    return LDPC36()


def test_parity_check_is_regular_and_generator_is_orthogonal(
    code: LDPC36,
) -> None:
    assert code.H.shape == (code.m, code.n)
    assert code.H.shape[0] == 504
    assert code.n == 1008
    assert code.rank == 502
    assert code.k == 506
    assert code.rate == pytest.approx(506 / 1008)
    assert np.all(code.H.sum(axis=1) == 6)
    assert np.all(code.H.sum(axis=0) == 3)
    assert np.count_nonzero((code.G @ code.H.T) % 2) == 0


def test_noiseless_encode_decode(code: LDPC36) -> None:
    rng = np.random.default_rng(4)
    data = rng.integers(0, 2, code.k * 4, dtype=np.uint8)
    encoded = code.encode(data)
    assert encoded.size == code.n * 4
    np.testing.assert_array_equal(code.decode(encoded), data)
    llr = (1.0 - 2.0 * encoded.astype(np.float64)) * 8.0
    np.testing.assert_array_equal(code.decode_soft(llr), data)


def test_bp_link_decodes_at_moderate_snr() -> None:
    result = run_link(
        Link(BPSK(), AWGNChannel(), LDPC36(), decoding="soft"),
        n_bits=100_000,
        ebn0_db=3.0,
        rng=np.random.default_rng(6),
    )
    assert result.decoding == "soft"
    assert result.ber < 1e-3


def test_invalid_lengths(code: LDPC36) -> None:
    with pytest.raises(ValueError, match="multiple"):
        code.encode(np.zeros(code.k - 1, dtype=np.uint8))
    with pytest.raises(ValueError, match="multiple"):
        code.decode(np.zeros(code.n - 1, dtype=np.uint8))
    with pytest.raises(ValueError, match="multiple"):
        code.decode_soft(np.zeros(code.n - 1))
