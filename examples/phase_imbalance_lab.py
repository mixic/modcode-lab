#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""I/Q phase imbalance and the tau involution, end to end.

Run from the repository root:

    python examples/phase_imbalance_lab.py

The script compares four detectors at ``phi = pi``:

* a standard coherent QPSK slicer that assumes ``phi = 0``;
* an oracle that knows ``phi = pi``;
* a pilot-assisted detector that estimates the hypothesis from known symbols;
* a blind detector that tests ``phi = 0`` against ``phi = pi`` from the data.

The blind detector is the control.  For QPSK with a uniform prior, the two
hypotheses have identical marginal likelihoods for every observation, so the
log-likelihood ratio is identically zero and no blind threshold can do better
than chance.  A pilot breaks that symmetry; an oracle assumes it is already
broken.  The figure and the printed report make the distinction explicit.

The output is ``examples/figures/phase_imbalance.png``.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from modcode.channel.awgn import AWGNChannel
from modcode.channel.phase_imbalance import PhaseImbalanceChannel
from modcode.modulation.involution import (
    block_tau_hypothesis,
    correct_tau,
    pilot_tau_hypothesis,
    tau_log_likelihood_ratio,
)
from modcode.modulation.qpsk import QPSK
from modcode.viz.phase_imbalance import (
    plot_channel_determinant,
    plot_phase_imbalance_ber,
    plot_phase_imbalance_constellation,
    plot_tau_llr,
)

plt.switch_backend("Agg")

OUTPUT_DIR = Path(__file__).resolve().parent / "figures"
RULE = "=" * 78
EBN0_DB = np.asarray([0.0, 2.0, 4.0, 6.0, 8.0])
N_BITS = 40_000
PHI = np.pi
PILOT_PERIOD = 32
PILOT = complex((1.0 + 1.0j) / np.sqrt(2.0))


def _ber(reference: np.ndarray, decoded: np.ndarray) -> float:
    if reference.shape != decoded.shape:
        raise ValueError("reference and decoded bits must have the same shape")
    return float(np.mean(reference != decoded))


def _coherent(
    qpsk: QPSK,
    bits: np.ndarray,
    ebn0_db: float,
    rng: np.random.Generator,
) -> np.ndarray:
    """Standard coherent detector: assume the phase reference is exact."""
    transmitted = qpsk.modulate(bits)
    received = PhaseImbalanceChannel(PHI).transmit(
        transmitted, ebn0_db, qpsk.bits_per_symbol, rng=rng
    )
    return qpsk.demodulate(received)


def _oracle(
    qpsk: QPSK,
    bits: np.ndarray,
    ebn0_db: float,
    rng: np.random.Generator,
) -> np.ndarray:
    """Oracle detector: knows that the channel is the phi = pi reflection."""
    transmitted = qpsk.modulate(bits)
    received = PhaseImbalanceChannel(PHI).transmit(
        transmitted, ebn0_db, qpsk.bits_per_symbol, rng=rng
    )
    return qpsk.demodulate_phase_imbalance(received, PHI)


def _pilot_assisted(
    qpsk: QPSK,
    bits: np.ndarray,
    ebn0_db: float,
    rng: np.random.Generator,
    period: int = PILOT_PERIOD,
) -> tuple[np.ndarray, float]:
    """Insert one known pilot per block and decide the hypothesis per block."""
    data = qpsk.modulate(bits)
    frames = []
    blocks = []
    for start in range(0, data.size, period):
        pilot_index = len(frames)
        frames.append(PILOT)
        block_start = len(frames)
        frames.extend(data[start : start + period])
        block_end = len(frames)
        blocks.append((pilot_index, block_start, block_end))
    frames = np.asarray(frames, dtype=np.complex128)
    received = PhaseImbalanceChannel(PHI).transmit(
        frames, ebn0_db, qpsk.bits_per_symbol, rng=rng
    )
    decoded_blocks = []
    wrong_hypotheses = 0
    for pilot_index, block_start, block_end in blocks:
        hypothesis = pilot_tau_hypothesis(
            received[[pilot_index]], np.asarray([PILOT])
        )
        wrong_hypotheses += int(hypothesis != 1)
        corrected = correct_tau(received[block_start:block_end], hypothesis)
        decoded_blocks.append(qpsk.demodulate(corrected))
    return (
        np.concatenate(decoded_blocks),
        wrong_hypotheses / len(blocks),
    )


def _blind(
    qpsk: QPSK,
    bits: np.ndarray,
    ebn0_db: float,
    rng: np.random.Generator,
    block: int = PILOT_PERIOD,
) -> tuple[np.ndarray, np.ndarray]:
    """Test the two hypotheses from the data alone, one block at a time."""
    transmitted = qpsk.modulate(bits)
    received = PhaseImbalanceChannel(PHI).transmit(
        transmitted, ebn0_db, qpsk.bits_per_symbol, rng=rng
    )
    variance = AWGNChannel.noise_variance(ebn0_db, qpsk.bits_per_symbol)
    decoded_blocks = []
    hypotheses = []
    for start in range(0, received.size, block):
        part = received[start : start + block]
        hypothesis = block_tau_hypothesis(
            qpsk, part, noise_variance=variance
        )
        hypotheses.append(hypothesis)
        decoded_blocks.append(
            qpsk.demodulate(correct_tau(part, hypothesis))
        )
    return np.concatenate(decoded_blocks), np.asarray(hypotheses)


def _noiseless_report(qpsk: QPSK, bits: np.ndarray) -> None:
    transmitted = qpsk.modulate(bits)
    received = PhaseImbalanceChannel(PHI).impair(transmitted)
    coherent = qpsk.demodulate(received)
    oracle = qpsk.demodulate_phase_imbalance(received, PHI)
    pilot_hypothesis = pilot_tau_hypothesis(
        np.asarray([np.conj(PILOT)]), np.asarray([PILOT])
    )
    pilot = qpsk.demodulate_tau_aware(received, hypothesis=pilot_hypothesis)
    llr = tau_log_likelihood_ratio(qpsk, received)
    blind_hypothesis = block_tau_hypothesis(qpsk, received)
    blind = qpsk.demodulate_tau_aware(received, hypothesis=blind_hypothesis)

    print()
    print(RULE)
    print("NOISELESS REFLECTION AT PHI = PI")
    print(RULE)
    print("channel matrix R(pi)      : [[1, 0], [0, -1]]")
    print("R(pi) @ R(pi)             : identity (involution)")
    print("R(pi) maps x to           : conj(x)")
    print(f"standard coherent BER     : {_ber(bits, coherent):.4f}")
    print(f"oracle tau-aware BER      : {_ber(bits, oracle):.4f}")
    print(
        "pilot-assisted BER        : "
        f"{_ber(bits, pilot):.4f} (hypothesis {pilot_hypothesis})"
    )
    print(
        "blind tau-aware BER       : "
        f"{_ber(bits, blind):.4f} (hypothesis {blind_hypothesis})"
    )
    print(f"max |tau LLR|             : {np.max(np.abs(llr)):.3e}")
    print()
    print("  The standard slicer inverts the Q bit, so exactly half the bits")
    print("  are wrong. The oracle and the pilot recover. The blind detector")
    print("  has no evidence: the QPSK likelihood is invariant under conj().")


def _finite_snr_report(qpsk: QPSK) -> dict[str, np.ndarray]:
    curves = {
        "coherent (phi=0)": [],
        "oracle (phi=pi)": [],
        "pilot-assisted": [],
        "blind (phi=0/pi)": [],
    }
    rows = []
    for index, ebn0_db in enumerate(EBN0_DB):
        bits = np.random.default_rng(1000 + index).integers(
            0, 2, size=N_BITS, dtype=np.uint8
        )
        coherent = _coherent(
            qpsk, bits, float(ebn0_db), np.random.default_rng(2000 + index)
        )
        oracle = _oracle(
            qpsk, bits, float(ebn0_db), np.random.default_rng(3000 + index)
        )
        pilot, pilot_error = _pilot_assisted(
            qpsk, bits, float(ebn0_db), np.random.default_rng(4000 + index)
        )
        blind, _ = _blind(
            qpsk, bits, float(ebn0_db), np.random.default_rng(5000 + index)
        )
        curves["coherent (phi=0)"].append(_ber(bits, coherent))
        curves["oracle (phi=pi)"].append(_ber(bits, oracle))
        curves["pilot-assisted"].append(_ber(bits, pilot))
        curves["blind (phi=0/pi)"].append(_ber(bits, blind))
        rows.append(
            (
                float(ebn0_db),
                curves["coherent (phi=0)"][-1],
                curves["oracle (phi=pi)"][-1],
                curves["pilot-assisted"][-1],
                curves["blind (phi=0/pi)"][-1],
                float(pilot_error),
            )
        )

    print()
    print(RULE)
    print("BER AT PHI = PI (QPSK, 40,000 bits per point)")
    print(RULE)
    print("Eb/N0   coherent   oracle   pilot    blind   pilot-hyp-err")
    for ebn0_db, coherent, oracle, pilot, blind, pilot_error in rows:
        print(
            f"{ebn0_db:5.1f}   {coherent:8.4f}   {oracle:7.4f}   "
            f"{pilot:7.4f}   {blind:6.4f}   {pilot_error:8.3f}"
        )
    print()
    print("  The coherent curve is pinned near 1/2: the reflection swaps Q,")
    print("  so the Q bit is wrong regardless of SNR. The oracle follows the")
    print("  ideal QPSK waterfall. The pilot approaches the oracle as its")
    print("  hypothesis test becomes reliable; it pays one pilot per block.")
    print("  The blind detector stays at 1/2 because its LLR is zero.")
    return {name: np.asarray(values) for name, values in curves.items()}


def _figure(qpsk: QPSK, curves: dict[str, np.ndarray]) -> None:
    bits = np.random.default_rng(7).integers(
        0, 2, size=20_000, dtype=np.uint8
    )
    transmitted = qpsk.modulate(bits)
    noisy = PhaseImbalanceChannel(PHI).transmit(
        transmitted,
        8.0,
        qpsk.bits_per_symbol,
        rng=np.random.default_rng(8),
    )
    llr = tau_log_likelihood_ratio(
        qpsk,
        noisy,
        noise_variance=AWGNChannel.noise_variance(8.0, qpsk.bits_per_symbol),
    )

    figure, axes = plt.subplots(2, 2, figsize=(11.0, 8.0))
    plot_phase_imbalance_constellation(axes[0, 0], qpsk, PHI)
    plot_phase_imbalance_ber(axes[0, 1], EBN0_DB, curves)
    plot_channel_determinant(axes[1, 0], np.linspace(0.0, 2.0 * np.pi, 721))
    plot_tau_llr(
        axes[1, 1],
        llr[:4000],
        title="Blind tau evidence at 8 dB (QPSK, phi = pi)",
    )
    figure.suptitle(
        "I/Q phase imbalance: a Q-axis reflection, not a quadrature collapse",
        fontsize=11,
    )
    figure.tight_layout()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    path = OUTPUT_DIR / "phase_imbalance.png"
    figure.savefig(path, dpi=140)
    plt.close(figure)
    print(f"  figure -> {path}")


def main() -> None:
    print(RULE)
    print("I/Q PHASE IMBALANCE AND THE TAU INVOLUTION")
    print(
        f"Python {sys.version_info.major}.{sys.version_info.minor}."
        f"{sys.version_info.micro} | numpy {np.__version__}"
    )
    print(RULE)

    qpsk = QPSK()
    bits = np.random.default_rng(42).integers(
        0, 2, size=N_BITS, dtype=np.uint8
    )
    _noiseless_report(qpsk, bits)
    curves = _finite_snr_report(qpsk)
    _figure(qpsk, curves)

    print()
    print(RULE)
    print("WHAT THE EXPERIMENT DOES AND DOES NOT SHOW")
    print(RULE)
    print("  phi = pi is a unitary reflection: invertible, energy preserving,")
    print("  and an involution. phi = pi/2 is the rank-collapse case; it loses")
    print("  the Q component and is not invertible.")
    print("  A blind detector cannot resolve phi=0 versus phi=pi for a")
    print("  symmetric constellation and prior. A pilot, a known preamble, an")
    print("  asymmetric constellation, or a code constraint can break the tie.")
    print("  This is an engineering analogue of the DSIN blindness result, not")
    print("  an identity between the two physical models and not evidence for RH.")
    print()
    print(f"done -- figure in {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
