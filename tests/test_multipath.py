# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Tests for the multipath/Doppler channel and the OFDM estimation path."""

import numpy as np
import pytest

from modcode.channel.multipath import (
    MultipathChannel,
    TapProcess,
    subcarrier_gains,
)
from modcode.fec.hamming import Hamming74
from modcode.modulation.qpsk import QPSK
from modcode.pipeline.ofdm import OFDMLink, interpolate_channel, run_ofdm_link


def test_tap_process_has_unit_average_power() -> None:
    process = TapProcess((0.7, 0.3), doppler_hz=1.0, block_symbols=1, seed=3)
    assert 0.0 <= process.correlation < 0.5
    taps = process.taps(4000)
    assert taps.shape == (4000, 2)
    measured = np.mean(np.abs(taps) ** 2, axis=0)
    np.testing.assert_allclose(measured, [0.7, 0.3], rtol=0.12)


def test_static_channel_has_constant_taps_and_fast_channel_does_not() -> None:
    static = TapProcess((0.5, 0.5), doppler_hz=0.0, block_symbols=1, seed=5)
    assert static.correlation == 1.0
    values = static.taps(20)
    np.testing.assert_allclose(values, np.repeat(values[:1], 20, axis=0))

    fast = TapProcess((0.5, 0.5), doppler_hz=2.0, block_symbols=1, seed=5)
    assert fast.correlation < 0.2
    moved = np.abs(fast.taps(20) - fast.taps(20))
    assert np.all(moved > 0.0)


def test_correlation_decreases_with_doppler() -> None:
    slow = TapProcess((1.0,), doppler_hz=0.1, block_symbols=10, seed=1)
    fast = TapProcess((1.0,), doppler_hz=1.0, block_symbols=10, seed=1)
    assert slow.correlation > fast.correlation >= 0.0


def test_multipath_channel_preserves_delay_structure() -> None:
    channel = MultipathChannel((0.7, 0.3), (0, 2), doppler_hz=0.0, seed=7)
    symbols = np.zeros(8, dtype=complex)
    symbols[2] = 1.0 + 0.0j
    received = channel.transmit(symbols, 0.0, np.random.default_rng(0))
    taps = channel.last_taps
    assert taps.shape == (8, 2)
    np.testing.assert_allclose(received[2], taps[2, 0])
    np.testing.assert_allclose(received[4], taps[4, 1])
    mask = np.ones(8, dtype=bool)
    mask[[2, 4]] = False
    np.testing.assert_allclose(received[mask], 0.0)


def test_multipath_channel_noise_variance_matches_moment() -> None:
    channel = MultipathChannel((1.0,), (0,), doppler_hz=0.0, seed=9)
    variance = 0.25
    received = channel.transmit(
        np.zeros(40_000, dtype=complex), variance, np.random.default_rng(1)
    )
    assert np.var(received.real) == pytest.approx(variance, rel=0.05)
    assert np.var(received.imag) == pytest.approx(variance, rel=0.05)


def test_subcarrier_gains_match_the_manual_dft() -> None:
    taps = np.asarray([1.0 + 0.0j, 0.5j])
    delays = (0, 2)
    n_fft = 8
    gains = subcarrier_gains(taps, delays, n_fft)
    expected = taps[0] + taps[1] * np.exp(
        -2j * np.pi * np.arange(n_fft) * 2 / n_fft
    )
    np.testing.assert_allclose(gains, expected)


def test_interpolation_is_exact_for_a_flat_channel() -> None:
    pilots = np.arange(0, 64, 4)
    gains = np.full(pilots.size, 0.3 - 0.4j)
    estimate = interpolate_channel(gains, pilots, 64)
    np.testing.assert_allclose(estimate, 0.3 - 0.4j)


def test_ofdm_link_noiseless_flat_channel() -> None:
    process = TapProcess((1.0,), doppler_hz=1.0, block_symbols=64, seed=0)
    for estimation in ("genie", "ls"):
        link = OFDMLink(
            QPSK(),
            process,
            (0,),
            n_fft=64,
            pilot_spacing=4,
            estimation=estimation,
        )
        result = run_ofdm_link(
            link, n_bits=20_000, ebn0_db=80.0, rng=np.random.default_rng(2)
        )
        assert result.ber == 0.0
        assert result.channel_nmse < 1e-6
        assert result.decoding == "hard"


def test_ofdm_link_with_fec_and_estimation() -> None:
    process = TapProcess((0.7, 0.3), doppler_hz=1.0, block_symbols=64, seed=4)
    link = OFDMLink(
        QPSK(),
        process,
        (0, 1),
        fec=Hamming74(),
        n_fft=64,
        pilot_spacing=4,
        estimation="ls",
    )
    result = run_ofdm_link(
        link, n_bits=20_000, ebn0_db=20.0, rng=np.random.default_rng(5)
    )
    assert result.decoding == "soft"
    assert result.ber < 1e-2
    assert 0.0 < result.channel_nmse < 1.0
    assert link.effective_rate == pytest.approx(0.75 * 4 / 7)


def test_estimation_costs_snr_and_doppler_decorrelates() -> None:
    def measure(process, estimation: str) -> float:
        link = OFDMLink(
            QPSK(),
            process,
            (0, 1),
            n_fft=64,
            pilot_spacing=4,
            estimation=estimation,
        )
        return run_ofdm_link(
            link, n_bits=40_000, ebn0_db=10.0, rng=np.random.default_rng(6)
        ).ber

    genie = measure(
        TapProcess((0.7, 0.3), doppler_hz=1.0, block_symbols=64, seed=8), "genie"
    )
    estimated = measure(
        TapProcess((0.7, 0.3), doppler_hz=1.0, block_symbols=64, seed=8), "ls"
    )
    assert estimated > genie > 0.0

    slow = TapProcess((1.0,), doppler_hz=0.0, block_symbols=64, seed=1)
    fast = TapProcess((1.0,), doppler_hz=1.0, block_symbols=64, seed=1)
    assert slow.correlation == 1.0
    assert fast.correlation < 0.2


def test_invalid_multipath_configurations() -> None:
    with pytest.raises(ValueError, match="same length"):
        MultipathChannel((0.5, 0.5), (0,))
    with pytest.raises(ValueError, match="non-negative"):
        MultipathChannel((1.0,), (-1,))
    with pytest.raises(ValueError, match="strictly positive"):
        TapProcess((0.0,))
    with pytest.raises(ValueError, match="finite and non-negative"):
        TapProcess((1.0,), doppler_hz=-1.0)

    process = TapProcess((1.0,), doppler_hz=0.0, block_symbols=64, seed=0)
    with pytest.raises(ValueError, match="data subcarrier"):
        OFDMLink(QPSK(), process, (0,), n_fft=64, pilot_spacing=1)
    with pytest.raises(ValueError, match="estimation"):
        OFDMLink(QPSK(), process, (0,), n_fft=64, estimation="perfect")
    with pytest.raises(ValueError, match="delays"):
        OFDMLink(QPSK(), process, (0, 1), n_fft=64)
    link = OFDMLink(QPSK(), process, (0,), n_fft=64)
    with pytest.raises(ValueError, match="positive"):
        run_ofdm_link(link, n_bits=0)
