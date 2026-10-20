# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Tests for the rate-1/2 convolutional code and Viterbi decoding."""

import numpy as np
import pytest

from modcode.channel.awgn import AWGNChannel
from modcode.fec.convolutional import ConvolutionalCode
from modcode.modulation.bpsk import BPSK
from modcode.pipeline.chain import Link, run_link


def test_convolutional_round_trip_hard_and_soft() -> None:
    code = ConvolutionalCode(block_bits=16)
    rng = np.random.default_rng(4)
    bits = rng.integers(0, 2, 16 * 20, dtype=np.uint8)
    encoded = code.encode(bits)
    np.testing.assert_array_equal(code.decode(encoded), bits)
    llr = 1.0 - 2.0 * encoded.astype(np.float64)
    np.testing.assert_array_equal(code.decode_soft(llr), bits)


def test_convolutional_single_error_is_corrected() -> None:
    code = ConvolutionalCode(block_bits=16)
    rng = np.random.default_rng(9)
    bits = rng.integers(0, 2, 16 * 5, dtype=np.uint8)
    encoded = code.encode(bits)
    corrupted = encoded.copy()
    corrupted[7] ^= 1
    np.testing.assert_array_equal(code.decode(corrupted), bits)
    llr = 1.0 - 2.0 * corrupted.astype(np.float64)
    np.testing.assert_array_equal(code.decode_soft(llr), bits)


def test_soft_viterbi_beats_hard_viterbi_under_noise() -> None:
    code = ConvolutionalCode(block_bits=32)
    rng = np.random.default_rng(42)
    bits = rng.integers(0, 2, 20_000, dtype=np.uint8)
    encoded = code.encode(bits)
    bpsk = BPSK()
    symbols = bpsk.modulate(encoded)
    variance = AWGNChannel.noise_variance(2.0, 1, code.rate)
    received = AWGNChannel().transmit(
        symbols, variance, np.random.default_rng(43)
    )
    hard = code.decode((received.real < 0).astype(np.uint8))
    soft = code.decode_soft(bpsk.soft_demodulate(received, variance))
    hard_ber = float(np.mean(hard != bits))
    soft_ber = float(np.mean(soft != bits))
    assert soft_ber < hard_ber
    assert soft_ber < 0.03


def test_convolutional_rate_and_block_lengths() -> None:
    code = ConvolutionalCode(block_bits=64)
    assert code.data_bits_per_block == 64
    assert code.encoded_bits_per_block == 132
    assert code.rate == pytest.approx(64 / 132)


def test_link_uses_soft_viterbi() -> None:
    result = run_link(
        Link(BPSK(), AWGNChannel(), ConvolutionalCode(block_bits=32)),
        n_bits=4_000,
        ebn0_db=3.0,
        rng=np.random.default_rng(11),
    )
    assert result.decoding == "soft"
    assert result.ber < 0.05


def test_convolutional_invalid_inputs() -> None:
    with pytest.raises(ValueError, match="block_bits"):
        ConvolutionalCode(block_bits=0)
    code = ConvolutionalCode(block_bits=8)
    with pytest.raises(ValueError, match="multiple"):
        code.encode(np.zeros(7, dtype=np.uint8))
    with pytest.raises(ValueError, match="whole"):
        code.decode(np.zeros(5, dtype=np.uint8))
    with pytest.raises(ValueError, match="finite"):
        code.decode_soft(np.asarray([np.inf]))
