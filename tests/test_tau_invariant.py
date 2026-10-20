# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Tests for the tau-invariant RM(1,3) code."""

import numpy as np
import pytest

from modcode.channel.awgn import AWGNChannel
from modcode.fec.tau_invariant import TauInvariantCode, rm13
from modcode.modulation.bpsk import BPSK
from modcode.pipeline.chain import Link, run_link


def test_rm13_parameters() -> None:
    code = rm13()
    assert code.n == 8
    assert code.k == 4
    assert code.rate == 0.5
    assert code.min_distance() == 4
    assert code.codebook.shape == (16, 8)


def test_rm13_is_tau_invariant() -> None:
    code = rm13()
    for codeword in code.codebook:
        assert tuple(code.tau_action(codeword)) in code._codebook_set


def test_rm13_is_not_cyclic() -> None:
    code = rm13()
    codeword = np.asarray([1, 1, 1, 1, 0, 0, 0, 0], dtype=np.uint8)
    assert tuple(codeword) in code._codebook_set
    assert tuple(np.roll(codeword, 1)) not in code._codebook_set


def test_rm13_round_trip_and_single_error_correction() -> None:
    code = rm13()
    rng = np.random.default_rng(0)
    bits = rng.integers(0, 2, code.k * 40, dtype=np.uint8)
    encoded = code.encode(bits)
    np.testing.assert_array_equal(code.decode(encoded), bits)
    for position in range(encoded.size):
        corrupted = encoded.copy()
        corrupted[position] ^= 1
        np.testing.assert_array_equal(code.decode(corrupted), bits)


def test_rm13_soft_decoding_round_trip() -> None:
    code = rm13()
    rng = np.random.default_rng(1)
    bits = rng.integers(0, 2, code.k * 20, dtype=np.uint8)
    encoded = code.encode(bits)
    llr = 1.0 - 2.0 * encoded.astype(np.float64)
    np.testing.assert_array_equal(code.decode_soft(llr), bits)


def test_rm13_link_uses_soft_decoding() -> None:
    result = run_link(
        Link(BPSK(), AWGNChannel(), rm13(), decoding="soft"),
        n_bits=20_000,
        ebn0_db=5.0,
        rng=np.random.default_rng(2),
    )
    assert result.decoding == "soft"
    assert result.ber < 0.02


def test_non_invariant_code_raises() -> None:
    generator = np.zeros((1, 8), dtype=np.uint8)
    generator[0, 0] = 1
    tau = np.asarray([index ^ 1 for index in range(8)])
    with pytest.raises(ValueError, match="not tau-invariant"):
        TauInvariantCode(generator, tau)


def test_bad_tau_permutation_raises() -> None:
    code = rm13()
    with pytest.raises(ValueError, match="involution"):
        TauInvariantCode(code.generator, np.roll(np.arange(8), 1))
