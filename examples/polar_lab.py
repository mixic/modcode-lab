#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Compare a short BI-AWGN Polar code with LDPC and convolutional coding.

Run from the repository root:

    python examples/polar_lab.py

The polar code is ``modcode.fec.polar.PolarAWGN``: Gaussian-approximation
construction, successive cancellation, list decoding, and CRC-aided list
decoding.  The second panel measures decoder work and wall-clock time, so the
list-size trade-off is visible rather than asserted.
"""

from __future__ import annotations

import time
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import NullFormatter

from modcode.channel.awgn import AWGNChannel
from modcode.fec.convolutional import ConvolutionalCode
from modcode.fec.crc import CRC8
from modcode.fec.ldpc import LDPC36
from modcode.fec.polar import PolarAWGN
from modcode.metrics.theory import theoretical_ber
from modcode.modulation.bpsk import BPSK
from modcode.pipeline.chain import Link, run_link
from modcode.viz.ber_plot import plot_ber

plt.switch_backend("Agg")

OUTPUT = Path(__file__).resolve().parent / "figures" / "polar.png"
EBN0_DB = np.asarray([0.0, 0.5, 1.0, 1.5, 2.0, 2.5])
N_BITS = 80_000
FLOOR = 1e-6
POLAR_N = 1024
POLAR_K = 512
DEMO_DB = 2.0
DESIGN_DB = 2.0


def _sweep(fec) -> np.ndarray:
    values = []
    for index, ebn0_db in enumerate(EBN0_DB):
        result = run_link(
            Link(BPSK(), AWGNChannel(), fec, decoding="soft"),
            n_bits=N_BITS,
            ebn0_db=float(ebn0_db),
            rng=np.random.default_rng(2026 + index),
        )
        values.append(result.ber)
    return np.asarray(values)


def _ber_at(fec, ebn0_db: float = DEMO_DB) -> float:
    """Measure one sweep point with the same seed the sweep would use."""
    index = int(np.searchsorted(EBN0_DB, ebn0_db))
    result = run_link(
        Link(BPSK(), AWGNChannel(), fec, decoding="soft"),
        n_bits=N_BITS,
        ebn0_db=ebn0_db,
        rng=np.random.default_rng(2026 + index),
    )
    return result.ber


def _benchmark(fec, blocks: int = 16, ebn0_db: float = DEMO_DB) -> tuple[float, int]:
    """Return (milliseconds per block, element-level decoder work per block)."""
    rng = np.random.default_rng(31)
    variance = AWGNChannel.noise_variance(ebn0_db, 1, fec.rate)
    data = rng.integers(0, 2, fec.data_bits_per_block * blocks, dtype=np.uint8)
    encoded = fec.encode(data)
    symbols = BPSK().modulate(encoded)
    received = symbols + np.sqrt(variance) * rng.standard_normal(symbols.size)
    llrs = BPSK().soft_demodulate(received, variance)
    blocks_of_llr = llrs.reshape(blocks, fec.encoded_bits_per_block)
    started = time.perf_counter()
    for block in blocks_of_llr:
        fec.decode_soft(block)
    milliseconds = (time.perf_counter() - started) / blocks * 1e3
    return milliseconds, int(getattr(fec, "last_work", 0))


def _spot_check(label: str, fec, ebn0_db: float, n_bits: int, seed: int) -> None:
    result = run_link(
        Link(BPSK(), AWGNChannel(), fec, decoding="soft"),
        n_bits=n_bits,
        ebn0_db=ebn0_db,
        rng=np.random.default_rng(seed),
    )
    print(f"  {label:<22} {result.ber:.3e}")


def main() -> None:
    curves = {
        "polar SC (L=1)": PolarAWGN(
            POLAR_N, POLAR_K, crc=None, list_size=1, design_ebn0_db=DESIGN_DB
        ),
        "polar SCL (L=4)": PolarAWGN(
            POLAR_N, POLAR_K, crc=None, list_size=4, design_ebn0_db=DESIGN_DB
        ),
        "polar CA-SCL (L=4, CRC-8)": PolarAWGN(
            POLAR_N, POLAR_K, crc=CRC8, list_size=4, design_ebn0_db=DESIGN_DB
        ),
        "LDPC(1008,506), min-sum": LDPC36(),
        "convolutional K=3, soft": ConvolutionalCode(block_bits=64),
    }
    results = {name: _sweep(fec) for name, fec in curves.items()}

    print(
        f"polar N={POLAR_N}, k={POLAR_K}, rate={POLAR_K / POLAR_N:.3f}; "
        f"{N_BITS} bits per point, one fixed seed per point\n"
    )
    print(f"{'Eb/N0':>5}" + "".join(f"{name:>27}" for name in results))
    for index, ebn0_db in enumerate(EBN0_DB):
        row = "".join(f"{values[index]:27.3e}" for values in results.values())
        print(f"{ebn0_db:5.1f}{row}")

    print("\ndecoder work and measured time at 2 dB")
    print(
        f"{'decoder':<28}{'list':>5}{'work/block':>12}"
        f"{'ms/block':>10}{'BER@2dB':>11}"
    )
    demo_index = int(np.searchsorted(EBN0_DB, DEMO_DB))
    benchmarks = [
        (
            "polar SC (L=1)",
            PolarAWGN(POLAR_N, POLAR_K, crc=None, list_size=1),
            results["polar SC (L=1)"][demo_index],
        ),
        (
            "polar SCL (L=2)",
            PolarAWGN(POLAR_N, POLAR_K, crc=None, list_size=2),
            _ber_at(PolarAWGN(POLAR_N, POLAR_K, crc=None, list_size=2)),
        ),
        (
            "polar SCL (L=4)",
            PolarAWGN(POLAR_N, POLAR_K, crc=None, list_size=4),
            results["polar SCL (L=4)"][demo_index],
        ),
        (
            "polar SCL (L=8)",
            PolarAWGN(POLAR_N, POLAR_K, crc=None, list_size=8),
            _ber_at(PolarAWGN(POLAR_N, POLAR_K, crc=None, list_size=8)),
        ),
        (
            "polar CA-SCL (L=2, CRC-8)",
            PolarAWGN(POLAR_N, POLAR_K, crc=CRC8, list_size=2),
            _ber_at(PolarAWGN(POLAR_N, POLAR_K, crc=CRC8, list_size=2)),
        ),
        (
            "polar CA-SCL (L=4, CRC-8)",
            PolarAWGN(POLAR_N, POLAR_K, crc=CRC8, list_size=4),
            results["polar CA-SCL (L=4, CRC-8)"][demo_index],
        ),
        (
            "polar CA-SCL (L=8, CRC-8)",
            PolarAWGN(POLAR_N, POLAR_K, crc=CRC8, list_size=8),
            _ber_at(PolarAWGN(POLAR_N, POLAR_K, crc=CRC8, list_size=8)),
        ),
        (
            "LDPC min-sum",
            LDPC36(),
            results["LDPC(1008,506), min-sum"][demo_index],
        ),
        (
            "convolutional",
            ConvolutionalCode(block_bits=64),
            results["convolutional K=3, soft"][demo_index],
        ),
    ]
    cost_points = {}
    for label, fec, ber in benchmarks:
        milliseconds, work = _benchmark(fec)
        list_size = getattr(fec, "list_size", "-")
        print(
            f"{label:<28}{list_size!s:>5}{work:>12,}{milliseconds:>10.1f}"
            f"{ber:>11.2e}"
        )
        cost_points[label] = (milliseconds, ber)

    print("\nblock length, SCL L=4, 50k bits at 1.5 dB")
    for block_length in (256, 512, 1024, 2048):
        _spot_check(
            f"N={block_length}",
            PolarAWGN(block_length, block_length // 2, crc=None, list_size=4),
            1.5,
            50_000,
            2031,
        )

    print("\nlist size, N=1024, 80k bits at 1.5 dB")
    for list_size in (1, 2, 4, 8, 16):
        _spot_check(
            f"L={list_size}",
            PolarAWGN(POLAR_N, POLAR_K, crc=None, list_size=list_size),
            1.5,
            80_000,
            2029,
        )

    print("\nconstruction design point, N=1024, SCL L=4, 50k bits at 1.5 dB")
    for design in (1.0, 2.0, 3.0):
        _spot_check(
            f"design {design:3.1f} dB",
            PolarAWGN(POLAR_N, POLAR_K, crc=None, list_size=4, design_ebn0_db=design),
            1.5,
            50_000,
            2030,
        )

    figure, (ber_axis, cost_axis) = plt.subplots(1, 2, figsize=(12.8, 5.0))
    plot_ber(
        EBN0_DB,
        {name: np.maximum(values, FLOOR) for name, values in results.items()},
        ber_axis,
        theory={"uncoded BPSK": theoretical_ber(BPSK(), EBN0_DB)},
    )
    ber_axis.set_ylim(FLOOR / 2.0, 1.0)
    ber_axis.set_title("Polar(1024,512), LDPC(1008,506), and convolutional")
    ber_axis.legend(fontsize=8)

    for label, (milliseconds, ber) in cost_points.items():
        cost_axis.scatter(milliseconds, max(ber, FLOOR), label=label, s=42)
    cost_axis.axhline(
        FLOOR, color="0.5", linestyle=":", linewidth=1.0,
        label="no errors in 80k bits",
    )
    cost_axis.set_xscale("log")
    cost_axis.set_yscale("log")
    cost_axis.set_xticks([2, 5, 10, 20, 50], labels=["2", "5", "10", "20", "50"])
    cost_axis.xaxis.set_minor_formatter(NullFormatter())
    cost_axis.set_xlabel("measured decoder time per block (ms)")
    cost_axis.set_ylabel(f"BER at {DEMO_DB:.1f} dB")
    cost_axis.set_ylim(FLOOR / 3.0, 0.3)
    cost_axis.set_title("Decoding cost versus error rate")
    cost_axis.grid(True, which="both", alpha=0.3)
    cost_axis.legend(fontsize=8)

    figure.tight_layout()
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(OUTPUT, dpi=140)
    plt.close(figure)
    print(f"\nfigure -> {OUTPUT}")


if __name__ == "__main__":
    main()
