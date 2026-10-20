# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Tests for the pedagogical BCH(15,7) and RS(15,11) codes."""

import numpy as np
import pytest

from modcode.channel.awgn import AWGNChannel
from modcode.fec.bch import BCH15
from modcode.fec.reed_solomon import ReedSolomon
from modcode.modulation.qpsk import QPSK
from modcode.pipeline.chain import Link, run_link


def test_bch_generator_and_parameters() -> None:
    code = BCH15()
    assert code.n == 15
    assert code.k == 7
    assert code.generator.size == 9
    assert code.rate == pytest.approx(7 / 15)


def test_bch_round_trip_and_all_double_errors() -> None:
    code = BCH15()
    rng = np.random.default_rng(4)
    data = rng.integers(0, 2, 7, dtype=np.uint8)
    codeword = code.encode(data)
    assert np.array_equal(code.decode(codeword), data)
    for first in range(15):
        for second in range(first + 1, 15):
            corrupted = codeword.copy()
            corrupted[[first, second]] ^= 1
            assert np.array_equal(code.decode(corrupted), data)


def test_bch_soft_decoding_round_trip() -> None:
    code = BCH15()
    rng = np.random.default_rng(5)
    data = rng.integers(0, 2, 7 * 20, dtype=np.uint8)
    codeword = code.encode(data)
    llr = 1.0 - 2.0 * codeword.astype(np.float64)
    assert np.array_equal(code.decode_soft(llr), data)


def test_bch_link_uses_soft_decoding() -> None:
    result = run_link(
        Link(QPSK(), AWGNChannel(), BCH15(), decoding="soft"),
        n_bits=14_000,
        ebn0_db=5.0,
        rng=np.random.default_rng(6),
    )
    assert result.decoding == "soft"
    assert result.ber < 0.02


def test_rs_round_trip_and_all_double_symbol_errors() -> None:
    code = ReedSolomon()
    rng = np.random.default_rng(7)
    data = rng.integers(0, 2, code.data_bits_per_block, dtype=np.uint8)
    codeword = code.encode(data)
    assert np.array_equal(code.decode(codeword), data)
    for first in range(code.n):
        for second in range(first + 1, code.n):
            corrupted = codeword.copy()
            corrupted[first * code.m : (first + 1) * code.m] ^= 1
            corrupted[second * code.m : (second + 1) * code.m] ^= 1
            assert np.array_equal(code.decode(corrupted), data)


def test_rs_link_falls_back_to_hard_decoding() -> None:
    result = run_link(
        Link(QPSK(), AWGNChannel(), ReedSolomon(), decoding="auto"),
        n_bits=8_000,
        ebn0_db=6.0,
        rng=np.random.default_rng(8),
    )
    assert result.decoding == "hard"
    assert result.ber < 0.1


def test_codes_reject_partial_blocks() -> None:
    with pytest.raises(ValueError, match="multiple"):
        BCH15().encode(np.zeros(6, dtype=np.uint8))
    with pytest.raises(ValueError, match="multiple"):
        ReedSolomon().decode(np.zeros(59, dtype=np.uint8))
