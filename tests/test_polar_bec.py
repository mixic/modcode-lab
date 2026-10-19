# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Tests for Polar transforms and exact BEC polarization calculations."""

import itertools

import numpy as np
import pytest

from modcode.fec.polar import (
    PolarBEC,
    bec_capacity,
    bec_synthetic_erasure_probabilities,
)


def test_two_channel_polarization_recursion() -> None:
    probabilities = bec_synthetic_erasure_probabilities(2, 0.5)
    assert np.allclose(probabilities, [0.75, 0.25])


@pytest.mark.parametrize("block_length", [1, 2, 4, 8, 32, 128])
def test_synthetic_channels_preserve_average_bec_capacity(block_length: int) -> None:
    erasure_probability = 0.37
    probabilities = bec_synthetic_erasure_probabilities(
        block_length, erasure_probability
    )
    assert np.isclose(np.mean(1.0 - probabilities), bec_capacity(erasure_probability))


def test_larger_blocks_polarize_more_synthetic_channels() -> None:
    erasure_probability = 0.4
    small = bec_synthetic_erasure_probabilities(16, erasure_probability)
    large = bec_synthetic_erasure_probabilities(1024, erasure_probability)
    small_extreme_fraction = np.mean((small < 0.01) | (small > 0.99))
    large_extreme_fraction = np.mean((large < 0.01) | (large > 0.99))
    assert large_extreme_fraction > small_extreme_fraction


def test_polar_encoder_and_sc_decoder_recover_all_short_messages_without_erasures() -> (
    None
):
    code = PolarBEC(block_length=8, information_bits=4, erasure_probability=0.35)
    for message in itertools.product((0, 1), repeat=code.information_bits):
        information = np.asarray(message, dtype=np.uint8)
        encoded = code.encode(information)
        assert np.array_equal(code.decode_bec(encoded), information)


def test_all_erased_bec_decodes_information_as_zero_fallback() -> None:
    code = PolarBEC(block_length=8, information_bits=4, erasure_probability=0.5)
    received = np.full(code.block_length, -1, dtype=np.int8)
    assert np.array_equal(code.decode_bec(received), np.zeros(4, dtype=np.uint8))


def test_sc_decoder_handles_deterministic_random_bec_erasures() -> None:
    rng = np.random.default_rng(22)
    code = PolarBEC(block_length=8, information_bits=4, erasure_probability=0.35)
    frame_errors = 0
    trials = 1000
    for _ in range(trials):
        message = rng.integers(0, 2, code.information_bits, dtype=np.uint8)
        received = code.encode(message).astype(np.int8)
        received[rng.random(code.block_length) < code.erasure_probability] = -1
        decoded = code.decode_bec(received)
        frame_errors += not np.array_equal(decoded, message)
    assert frame_errors / trials < 0.2


def test_polar_code_validates_parameters_and_received_values() -> None:
    with pytest.raises(ValueError):
        PolarBEC(block_length=6, information_bits=3, erasure_probability=0.5)
    with pytest.raises(ValueError):
        PolarBEC(block_length=8, information_bits=9, erasure_probability=0.5)
    code = PolarBEC(block_length=8, information_bits=4, erasure_probability=0.5)
    with pytest.raises(ValueError):
        code.decode_bec(np.full(8, 2, dtype=np.int8))
