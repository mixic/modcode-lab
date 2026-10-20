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
from modcode.metrics.ber import ber_sweep, ber_sweep_detailed
from modcode.metrics.theory import theoretical_ber
from modcode.modulation.bpsk import BPSK
from modcode.modulation.qam import QAM
from modcode.modulation.qpsk import QPSK
from modcode.pipeline.chain import Link, run_link


def test_awgn_variance_uses_bits_per_symbol_and_code_rate() -> None:
    assert np.isclose(AWGNChannel.noise_variance(0.0, 2), 0.25)
    assert np.isclose(AWGNChannel.noise_variance(0.0, 2, 4 / 7), 7 / 16)


def test_awgn_empirical_variance_matches_requested_variance() -> None:
    channel = AWGNChannel()
    expected = channel.noise_variance(3.0, 2)
    noise = channel.transmit(
        np.zeros(100_000, dtype=np.complex128),
        expected,
        rng=np.random.default_rng(1),
    )
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


def test_simulated_bpsk_ber_matches_the_q_function() -> None:
    n_bits = 400_000
    ebn0_db = 3.0
    result = run_link(
        Link(BPSK(), AWGNChannel()),
        n_bits=n_bits,
        ebn0_db=ebn0_db,
        rng=np.random.default_rng(2026),
    )
    expected = theoretical_ber(BPSK(), [ebn0_db])[0]
    standard_error = np.sqrt(expected * (1.0 - expected) / n_bits)
    assert abs(result.ber - expected) < 4.0 * standard_error


def test_simulated_qpsk_ber_matches_the_q_function() -> None:
    n_bits = 400_000
    ebn0_db = 3.0
    result = run_link(
        Link(QPSK(), AWGNChannel()),
        n_bits=n_bits,
        ebn0_db=ebn0_db,
        rng=np.random.default_rng(2026),
    )
    expected = theoretical_ber(QPSK(), [ebn0_db])[0]
    standard_error = np.sqrt(expected * (1.0 - expected) / n_bits)
    assert abs(result.ber - expected) < 4.0 * standard_error


def test_ber_sweep_detailed_uses_stopping_rule_and_reports_intervals() -> None:
    curves = ber_sweep_detailed(
        {"BPSK": BPSK()},
        [4.0],
        min_errors=50,
        max_bits=400_000,
        batch_bits=20_000,
        seed=11,
    )
    curve = curves["BPSK"]
    assert curve.errors[0] >= 50 or curve.bits[0] == 400_000
    assert curve.bits[0] <= 400_000
    assert 0.0 <= curve.ci_low[0] <= curve.ber[0] <= curve.ci_high[0] <= 1.0


def test_zero_error_point_has_a_nonzero_upper_bound() -> None:
    curves = ber_sweep_detailed(
        {"BPSK": BPSK()},
        [20.0],
        min_errors=10,
        max_bits=2_000,
        batch_bits=1_000,
        seed=5,
    )
    curve = curves["BPSK"]
    assert curve.errors[0] == 0
    assert curve.ber[0] == 0.0
    assert curve.ci_high[0] > 0.0


def test_link_streams_do_not_depend_on_other_modulations() -> None:
    both = ber_sweep(
        {"BPSK": BPSK(), "QPSK": QPSK()},
        [2.0],
        n_bits=20_000,
        seed=7,
    )
    bpsk_only = ber_sweep(
        {"BPSK": BPSK()},
        [2.0],
        n_bits=20_000,
        seed=7,
    )
    qpsk_only = ber_sweep(
        {"QPSK": QPSK()},
        [2.0],
        n_bits=20_000,
        seed=7,
    )
    assert both["BPSK"][0] == bpsk_only["BPSK"][0]
    assert both["QPSK"][0] == qpsk_only["QPSK"][0]


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
