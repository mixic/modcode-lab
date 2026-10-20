# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Tests for block Rayleigh fading."""

import numpy as np
import pytest

from modcode.channel.rayleigh import BlockRayleighChannel
from modcode.metrics.theory import theoretical_rayleigh_ber
from modcode.modulation.bpsk import BPSK
from modcode.pipeline.chain import Link, run_link


def test_equalized_fading_is_identity_without_noise() -> None:
    channel = BlockRayleighChannel(block_length=3)
    symbols = np.asarray([1 + 1j, -1 + 0.5j, 0.25 - 0.75j, 2 + 0j])
    received = channel.transmit(
        symbols, 0.0, rng=np.random.default_rng(2)
    )
    np.testing.assert_allclose(received, symbols, atol=1e-12)
    assert channel.last_coefficients is not None
    assert channel.last_coefficients.shape == symbols.shape


def test_coefficients_are_constant_within_a_block() -> None:
    channel = BlockRayleighChannel(block_length=3)
    received = channel.transmit(
        np.ones(9, dtype=complex), 0.0, rng=np.random.default_rng(3)
    )
    coefficients = channel.last_coefficients
    assert coefficients is not None
    assert np.allclose(coefficients[0:3], coefficients[0])
    assert np.allclose(coefficients[3:6], coefficients[3])
    assert np.allclose(coefficients[6:9], coefficients[6])
    assert np.allclose(received, 1.0)


def test_rayleigh_bpsk_ber_matches_known_reference() -> None:
    n_bits = 200_000
    ebn0_db = 10.0
    result = run_link(
        Link(
            BPSK(),
            BlockRayleighChannel(block_length=10),
            decoding="hard",
        ),
        n_bits=n_bits,
        ebn0_db=ebn0_db,
        rng=np.random.default_rng(4),
    )
    expected = theoretical_rayleigh_ber(BPSK(), [ebn0_db])[0]
    standard_error = np.sqrt(expected * (1.0 - expected) / n_bits)
    assert abs(result.ber - expected) < 4.0 * standard_error


def test_rayleigh_theory_value() -> None:
    assert np.isclose(
        theoretical_rayleigh_ber(BPSK(), [0.0])[0],
        0.5 * (1.0 - np.sqrt(0.5)),
    )


def test_rayleigh_validation() -> None:
    with pytest.raises(ValueError, match="block_length"):
        BlockRayleighChannel(block_length=0)
    with pytest.raises(ValueError, match="noise_variance"):
        BlockRayleighChannel().transmit(
            np.ones(2), -1.0, rng=np.random.default_rng(1)
        )
    with pytest.raises(ValueError, match="Rayleigh"):
        theoretical_rayleigh_ber(type("Other", (), {})(), [0.0])
