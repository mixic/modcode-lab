#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Lattice modulation appendix: BER, PAPR diagnostic, and complexity.

Run from the repository root:

    python examples/lattice_lab.py

The BER comparison is at equal energy per real dimension and equal bits per
real dimension, not at equal bits per symbol.  That distinction is essential:
D4-16 carries 4 bits over 4 real dimensions, so it is comparable to QPSK
(2 bits over 2 dimensions), not to 16-QAM (4 bits over 2 dimensions).

The PAPR panel is a component-stream diagnostic, not a full system
comparison; see docs/lattice-modulation.md.
"""

from __future__ import annotations

import time
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from modcode.metrics.lattice import (
    qam_per_dimension_points,
    simulate_lattice_ber,
    simulate_point_ber,
)
from modcode.modulation.base import Modulation
from modcode.modulation.lattice import LatticeModulation
from modcode.modulation.qam import QAM
from modcode.modulation.qpsk import QPSK
from modcode.viz.ber_plot import plot_ber

plt.switch_backend("Agg")

OUTPUT = Path(__file__).resolve().parent / "figures" / "lattice.png"
EBN0_DB = np.asarray([2.0, 4.0, 6.0, 8.0, 10.0, 12.0])
N_BITS = 100_000
FLOOR = 1e-6


def _qam_ber(modulation: Modulation, ebn0_db: float, rng) -> float:
    points, labels = qam_per_dimension_points(modulation)
    return simulate_point_ber(
        points,
        labels,
        modulation.bits_per_symbol,
        ebn0_db,
        N_BITS,
        rng,
    )


def _ber_curves() -> dict[str, np.ndarray]:
    qpsk = QPSK()
    d4_16 = LatticeModulation("D4", 4, labeler="gray")
    e8_256 = LatticeModulation("E8", 8, labeler="gray")
    qam16 = QAM(16)
    d4_256 = LatticeModulation("D4", 8, labeler="gray")

    curves = {
        "QPSK (eta=1)": [],
        "D4-16 (eta=1)": [],
        "E8-256 (eta=1)": [],
        "16-QAM (eta=2)": [],
        "D4-256 (eta=2)": [],
    }
    for index, ebn0_db in enumerate(EBN0_DB):
        curves["QPSK (eta=1)"].append(
            _qam_ber(qpsk, float(ebn0_db), np.random.default_rng(2026 + index))
        )
        curves["D4-16 (eta=1)"].append(
            simulate_lattice_ber(
                d4_16,
                float(ebn0_db),
                N_BITS,
                np.random.default_rng(2026 + index),
            )
        )
        curves["E8-256 (eta=1)"].append(
            simulate_lattice_ber(
                e8_256,
                float(ebn0_db),
                N_BITS,
                np.random.default_rng(2026 + index),
            )
        )
        curves["16-QAM (eta=2)"].append(
            _qam_ber(qam16, float(ebn0_db), np.random.default_rng(2026 + index))
        )
        curves["D4-256 (eta=2)"].append(
            simulate_lattice_ber(
                d4_256,
                float(ebn0_db),
                N_BITS,
                np.random.default_rng(2026 + index),
            )
        )
    return {name: np.asarray(values) for name, values in curves.items()}


def _papr_db(components: np.ndarray) -> float:
    waveform = np.fft.ifft(components)
    power = np.abs(waveform) ** 2
    return float(10.0 * np.log10(power.max() / power.mean()))


def _block_papr(modulation, n_symbols: int, rng) -> float:
    bits = rng.integers(
        0, 2, n_symbols * modulation.bits_per_symbol, dtype=np.uint8
    )
    symbols = modulation.modulate(bits)
    if symbols.ndim == 1:
        symbols = symbols[:, None]
    return max(_papr_db(symbols[:, dimension]) for dimension in range(symbols.shape[1]))


def _papr_curves(blocks: int = 2_000) -> tuple[dict[str, np.ndarray], dict[str, float]]:
    n_symbols = 64
    modulations = {
        "256-QAM": QAM(256),
        "D4-256": LatticeModulation("D4", 8),
        "E8-256": LatticeModulation("E8", 8),
    }
    thresholds = np.linspace(0.0, 18.0, 73)
    ccdf: dict[str, np.ndarray] = {}
    percentiles: dict[str, float] = {}
    for index, (name, modulation) in enumerate(modulations.items()):
        rng = np.random.default_rng(7000 + index)
        values = np.asarray(
            [
                _block_papr(modulation, n_symbols, rng)
                for _ in range(blocks)
            ]
        )
        ccdf[name] = np.asarray(
            [np.mean(values > threshold) for threshold in thresholds]
        )
        percentiles[name] = float(np.percentile(values, 99.9))
    return ccdf, percentiles


def _throughput(modulation, n_symbols: int = 2_000) -> float:
    rng = np.random.default_rng(9000)
    bits = rng.integers(
        0, 2, n_symbols * modulation.bits_per_symbol, dtype=np.uint8
    )
    symbols = modulation.modulate(bits)
    modulation.demodulate(symbols[:50])
    start = time.perf_counter()
    modulation.demodulate(symbols)
    elapsed = time.perf_counter() - start
    return n_symbols / elapsed


def _eta(modulation) -> float:
    if isinstance(modulation, LatticeModulation):
        return modulation.bits_per_real_dimension
    return modulation.bits_per_symbol / 2.0


def main() -> None:
    curves = _ber_curves()

    print("minimum distances at unit energy per real dimension")
    for name, modulation in (
        ("QPSK", QPSK()),
        ("D4-16", LatticeModulation("D4", 4)),
        ("E8-256", LatticeModulation("E8", 8)),
        ("16-QAM", QAM(16)),
        ("D4-256", LatticeModulation("D4", 8)),
    ):
        if not isinstance(modulation, LatticeModulation):
            points, _ = qam_per_dimension_points(modulation)
            distances = np.sum(
                (points[:, None] - points[None, :]) ** 2, axis=2
            )
            np.fill_diagonal(distances, np.inf)
            dmin = float(np.sqrt(np.min(distances)))
        else:
            dmin = modulation.min_distance(1.0)
        print(
            f"  {name:8s}  eta={_eta(modulation):.2f}  d_min={dmin:.4f}"
        )

    print("BER")
    for index, ebn0_db in enumerate(EBN0_DB):
        row = "  ".join(
            f"{name}={curves[name][index]:.3e}" for name in curves
        )
        print(f"  {ebn0_db:5.1f} dB  {row}")

    ccdf, percentiles = _papr_curves()
    print("PAPR diagnostic, 64 symbols/block, max across component streams")
    for name, value in percentiles.items():
        print(f"  {name:8s}  99.9% = {value:.2f} dB")

    throughput = {
        "256-QAM": _throughput(QAM(256)),
        "D4-256": _throughput(LatticeModulation("D4", 8)),
        "E8-256": _throughput(LatticeModulation("E8", 8)),
    }
    print("naive nearest-neighbour demodulation throughput")
    for name, value in throughput.items():
        print(f"  {name:8s}  {value:.0f} symbols/s")

    figure, axes = plt.subplots(1, 3, figsize=(15.0, 4.4))
    plot_ber(
        EBN0_DB,
        {name: np.maximum(values, FLOOR) for name, values in curves.items()},
        axes[0],
    )
    axes[0].set_title("BER at matched bits and energy per dimension")
    axes[0].set_ylim(FLOOR / 2.0, 1.0)

    thresholds = np.linspace(0.0, 18.0, 73)
    for name, values in ccdf.items():
        axes[1].semilogy(thresholds, values, label=name)
    axes[1].set_xlabel("PAPR threshold (dB)")
    axes[1].set_ylabel("CCDF")
    axes[1].set_title("Component-stream PAPR diagnostic")
    axes[1].grid(True, which="both", alpha=0.3)
    axes[1].legend(fontsize=7)

    names = list(throughput)
    axes[2].bar(names, [throughput[name] for name in names], color="tab:gray")
    axes[2].set_yscale("log")
    axes[2].set_ylabel("symbols / s (log scale)")
    axes[2].set_title("Software demodulation cost")
    axes[2].tick_params(axis="x", rotation=20)

    figure.suptitle("Lattice appendix: coding gain, PAPR, and complexity")
    figure.tight_layout()
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(OUTPUT, dpi=140)
    plt.close(figure)
    print(f"figure -> {OUTPUT}")


if __name__ == "__main__":
    main()
