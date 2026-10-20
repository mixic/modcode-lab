# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Tests for the I/Q phase-imbalance channel and the tau ambiguity."""

import numpy as np
import pytest

from modcode.channel.awgn import AWGNChannel
from modcode.channel.phase_imbalance import PhaseImbalanceChannel
from modcode.fec.crc import CRC8
from modcode.modulation.involution import (
    block_tau_hypothesis,
    correct_tau,
    demodulate_phase_imbalance,
    demodulate_tau_aware,
    pilot_tau_hypothesis,
    tau_log_likelihood_ratio,
)
from modcode.modulation.qpsk import QPSK


def _bits() -> np.ndarray:
    return np.asarray([0, 1, 1, 0, 1, 0, 0, 1], dtype=np.uint8)


def test_phase_imbalance_identity_and_pi_reflection() -> None:
    symbols = np.asarray([1.0 + 2.0j, -0.5 + 0.25j])
    np.testing.assert_allclose(
        PhaseImbalanceChannel(0.0).impair(symbols), symbols
    )
    np.testing.assert_allclose(
        PhaseImbalanceChannel(np.pi).impair(symbols),
        np.conj(symbols),
        atol=1e-12,
    )


def test_phase_imbalance_matrix_at_pi_is_an_involution() -> None:
    matrix = PhaseImbalanceChannel(np.pi).matrix()
    np.testing.assert_allclose(matrix, np.diag([1.0, -1.0]), atol=1e-12)
    np.testing.assert_allclose(matrix @ matrix, np.eye(2), atol=1e-12)
    assert np.isclose(np.linalg.det(matrix), -1.0)


def test_quarter_turn_phase_imbalance_is_rank_deficient() -> None:
    matrix = PhaseImbalanceChannel(np.pi / 2.0).matrix()
    assert np.isclose(np.linalg.det(matrix), 0.0, atol=1e-12)
    assert np.linalg.matrix_rank(matrix) == 1
    with pytest.raises(ValueError, match="singular"):
        demodulate_phase_imbalance(
            QPSK(), np.asarray([1.0 + 1.0j]), np.pi / 2.0
        )


def test_known_phase_correction_recovers_qpsk() -> None:
    qpsk = QPSK()
    bits = _bits()
    transmitted = qpsk.modulate(bits)
    received = PhaseImbalanceChannel(np.pi).impair(transmitted)
    decoded = demodulate_phase_imbalance(qpsk, received, np.pi)
    np.testing.assert_array_equal(decoded, bits)


def test_known_phase_correction_handles_a_non_pi_angle() -> None:
    qpsk = QPSK()
    bits = _bits()
    transmitted = qpsk.modulate(bits)
    received = PhaseImbalanceChannel(0.7).impair(transmitted)
    decoded = demodulate_phase_imbalance(qpsk, received, 0.7)
    np.testing.assert_array_equal(decoded, bits)


def test_standard_detector_flips_only_the_q_bit_at_pi() -> None:
    qpsk = QPSK()
    bits = _bits()
    received = np.conj(qpsk.modulate(bits))
    decoded = qpsk.demodulate(received)
    np.testing.assert_array_equal(decoded[0::2], bits[0::2])
    np.testing.assert_array_equal(decoded[1::2], 1 - bits[1::2])
    assert np.mean(decoded != bits) == 0.5


def test_tau_aware_demodulation_uses_the_hypothesis() -> None:
    qpsk = QPSK()
    bits = _bits()
    transmitted = qpsk.modulate(bits)
    np.testing.assert_array_equal(
        demodulate_tau_aware(qpsk, transmitted, hypothesis=0), bits
    )
    np.testing.assert_array_equal(
        demodulate_tau_aware(qpsk, np.conj(transmitted), hypothesis=1), bits
    )
    with pytest.raises(ValueError, match="hypothesis"):
        demodulate_tau_aware(qpsk, transmitted, hypothesis=2)


def test_tau_likelihood_ratio_is_zero_for_qpsk() -> None:
    rng = np.random.default_rng(3)
    received = rng.standard_normal(100) + 1j * rng.standard_normal(100)
    llr = tau_log_likelihood_ratio(QPSK(), received, noise_variance=0.7)
    np.testing.assert_allclose(llr, 0.0, atol=1e-12)


def test_reflected_qpsk_gives_a_tie_and_chance_level_ber() -> None:
    qpsk = QPSK()
    bits = _bits()
    received = np.conj(qpsk.modulate(bits))
    llr = tau_log_likelihood_ratio(qpsk, received)
    assert abs(float(np.sum(llr))) < 1e-9
    assert block_tau_hypothesis(qpsk, received) == 0
    decoded = demodulate_tau_aware(qpsk, received, hypothesis=0)
    assert np.mean(decoded != bits) == 0.5


def test_pilot_selects_the_reflection() -> None:
    pilot = np.asarray([(1.0 + 1.0j) / np.sqrt(2.0)])
    assert pilot_tau_hypothesis(pilot, pilot) == 0
    assert pilot_tau_hypothesis(np.conj(pilot), pilot) == 1


def test_crc_constraint_breaks_the_blind_tie() -> None:
    payload = np.asarray(
        [1, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1],
        dtype=np.uint8,
    )
    codeword = CRC8.encode(payload)
    qpsk = QPSK()
    received = PhaseImbalanceChannel(np.pi).impair(qpsk.modulate(codeword))
    ideal = qpsk.demodulate(received)
    reflected = qpsk.demodulate(np.conj(received))
    assert not CRC8.check(ideal)
    assert CRC8.check(reflected)
    np.testing.assert_array_equal(reflected, codeword)


def test_correct_tau_and_validation() -> None:
    values = np.asarray([1.0 + 2.0j])
    np.testing.assert_allclose(correct_tau(values, 0), values)
    np.testing.assert_allclose(correct_tau(values, 1), np.conj(values))
    with pytest.raises(ValueError, match="hypothesis"):
        correct_tau(values, 2)
    with pytest.raises(ValueError, match="one-dimensional"):
        correct_tau(np.zeros((2, 2)), 0)


def test_phase_imbalance_channel_adds_awgn() -> None:
    channel = PhaseImbalanceChannel(np.pi)
    expected = AWGNChannel.noise_variance(3.0, 2)
    received = channel.transmit(
        np.zeros(100_000, dtype=np.complex128),
        expected,
        rng=np.random.default_rng(1),
    )
    assert np.isclose(np.var(received.real), expected, rtol=0.03)
    assert np.isclose(np.var(received.imag), expected, rtol=0.03)


def test_phase_imbalance_energy_and_validation() -> None:
    symbols = np.asarray([1.0 + 1.0j])
    assert np.isclose(
        PhaseImbalanceChannel(np.pi).mean_energy_ratio(symbols), 1.0
    )
    assert np.isclose(
        PhaseImbalanceChannel(np.pi / 2.0).mean_energy_ratio(symbols), 0.0
    )
    assert np.isclose(
        PhaseImbalanceChannel.from_degrees(180.0).phi, np.pi
    )
    with pytest.raises(ValueError, match="finite"):
        PhaseImbalanceChannel(np.inf)
    with pytest.raises(ValueError, match="one-dimensional"):
        PhaseImbalanceChannel(0.0).impair(np.zeros((2, 2)))
    with pytest.raises(ValueError, match="all zero"):
        PhaseImbalanceChannel(0.0).mean_energy_ratio(np.zeros(2))
