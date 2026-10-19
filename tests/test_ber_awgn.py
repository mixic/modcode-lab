# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Tests for AWGN normalization and link-level BER behavior."""

import numpy as np

from modcode.channel.awgn import AWGNChannel
from modcode.fec.hamming import Hamming74
from modcode.metrics.ber import ber_sweep
from modcode.modulation.bpsk import BPSK
from modcode.modulation.qam import QAM
from modcode.modulation.qpsk import QPSK
from modcode.pipeline.chain import Link, run_link


def test_awgn_variance_uses_bits_per_symbol_and_code_rate() -> None:
    assert np.isclose(AWGNChannel.noise_variance(0.0, 2), 0.25)
    assert np.isclose(AWGNChannel.noise_variance(0.0, 2, 4 / 7), 7 / 16)


def test_awgn_empirical_variance_matches_requested_variance() -> None:
    channel = AWGNChannel()
    noise = channel.transmit(
        np.zeros(100_000, dtype=np.complex128),
        ebn0_db=3.0,
        bits_per_symbol=2,
        rng=np.random.default_rng(1),
    )
    expected = channel.noise_variance(3.0, 2)
    assert np.isclose(np.var(noise.real), expected, rtol=0.03)
    assert np.isclose(np.var(noise.imag), expected, rtol=0.03)


def test_bpsk_ber_decreases_with_ebn0() -> None:
    results = ber_sweep(
        {"BPSK": BPSK()},
        [0.0, 2.0, 4.0],
        n_bits=100_000,
        seed=123,
    )
    ber = results["BPSK"]
    assert ber[0] > ber[1] > ber[2]


def test_link_with_hamming_returns_original_bit_count() -> None:
    result = run_link(
        Link(BPSK(), AWGNChannel(), Hamming74()),
        n_bits=101,
        ebn0_db=12.0,
        rng=np.random.default_rng(7),
    )
    assert result.source_bits.size == 101
    assert result.received_bits.size == 101
    assert result.ber < 0.1


def test_sweep_supports_all_documented_modulations() -> None:
    results = ber_sweep(
        {"BPSK": BPSK(), "QPSK": QPSK(), "16-QAM": QAM(16)},
        [0.0, 4.0],
        n_bits=2_000,
        seed=4,
    )
    assert set(results) == {"BPSK", "QPSK", "16-QAM"}
    assert all(values.shape == (2,) for values in results.values())
