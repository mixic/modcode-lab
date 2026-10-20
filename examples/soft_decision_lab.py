#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Compare hard and soft decoding for Hamming and convolutional codes.

Run from the repository root:

    python examples/soft_decision_lab.py

The figure separates three questions: whether uncoded BPSK matches the
Q-function, whether Hamming(7,4) benefits from LLRs, and how much a
rate-1/2 convolutional code gains from soft Viterbi decoding.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from modcode.channel.awgn import AWGNChannel
from modcode.fec.convolutional import ConvolutionalCode
from modcode.fec.hamming import Hamming74
from modcode.fec.tau_invariant import rm13
from modcode.metrics.theory import theoretical_ber
from modcode.modulation.bpsk import BPSK
from modcode.pipeline.chain import Link, run_link
from modcode.viz.ber_plot import plot_ber

plt.switch_backend("Agg")

OUTPUT = Path(__file__).resolve().parent / "figures" / "soft_decision.png"
EBN0_DB = np.asarray([0.0, 2.0, 4.0, 6.0])
N_BITS = 80_000


def _sweep(fec, decoding: str) -> np.ndarray:
    values = []
    for index, ebn0_db in enumerate(EBN0_DB):
        result = run_link(
            Link(
                BPSK(),
                AWGNChannel(),
                fec,
                decoding=decoding,
            ),
            n_bits=N_BITS,
            ebn0_db=float(ebn0_db),
            rng=np.random.default_rng(2026 + 10 * index + len(decoding)),
        )
        values.append(result.ber)
    return np.asarray(values)


def main() -> None:
    hamming = Hamming74()
    rm = rm13()
    convolutional = ConvolutionalCode(block_bits=32)

    uncoded = _sweep(None, "auto")
    hamming_hard = _sweep(hamming, "hard")
    hamming_soft = _sweep(hamming, "soft")
    rm_hard = _sweep(rm, "hard")
    rm_soft = _sweep(rm, "soft")
    convolutional_hard = _sweep(convolutional, "hard")
    convolutional_soft = _sweep(convolutional, "soft")

    print(
        "Eb/N0   uncoded   Ham-hard  Ham-soft  RM-hard  RM-soft  "
        "Conv-hard  Conv-soft"
    )
    for index, ebn0_db in enumerate(EBN0_DB):
        print(
            f"{ebn0_db:5.1f}   {uncoded[index]:8.5f}  "
            f"{hamming_hard[index]:8.5f}  {hamming_soft[index]:8.5f}  "
            f"{rm_hard[index]:7.5f}  {rm_soft[index]:7.5f}  "
            f"{convolutional_hard[index]:9.5f}  "
            f"{convolutional_soft[index]:9.5f}"
        )

    figure, axes = plt.subplots(2, 2, figsize=(11.0, 8.0))
    floor = 1e-6
    plot_ber(
        EBN0_DB,
        {"uncoded BPSK": np.maximum(uncoded, floor)},
        axes[0, 0],
        theory={"BPSK": theoretical_ber(BPSK(), EBN0_DB)},
    )
    axes[0, 0].set_title("Uncoded BPSK vs Q-function")

    plot_ber(
        EBN0_DB,
        {
            "Hamming hard": np.maximum(hamming_hard, floor),
            "Hamming soft": np.maximum(hamming_soft, floor),
        },
        axes[0, 1],
    )
    axes[0, 1].set_title("Cyclic Hamming(7,4): hard vs soft")

    plot_ber(
        EBN0_DB,
        {
            "RM(1,3) hard": np.maximum(rm_hard, floor),
            "RM(1,3) soft FHT": np.maximum(rm_soft, floor),
        },
        axes[1, 0],
    )
    axes[1, 0].set_title("Tau-invariant RM(1,3): hard vs FHT")

    plot_ber(
        EBN0_DB,
        {
            "convolutional hard": np.maximum(convolutional_hard, floor),
            "convolutional soft": np.maximum(convolutional_soft, floor),
        },
        axes[1, 1],
    )
    axes[1, 1].set_title("Rate-1/2 convolutional: hard vs soft")
    axes[1, 1].set_yscale("log")

    for axis in axes.flat:
        axis.set_ylim(floor / 2.0, 1.0)

    figure.suptitle("Soft-decision decoding and the theoretical BER sanity check")
    figure.tight_layout()
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(OUTPUT, dpi=140)
    plt.close(figure)
    print(f"figure -> {OUTPUT}")
    print("  Values at the plotting floor are zero-error runs, not measured BER.")


if __name__ == "__main__":
    main()
