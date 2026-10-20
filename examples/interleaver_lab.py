#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Show what a bit interleaver buys for Hamming + QPSK.

Run from the repository root:

    python examples/interleaver_lab.py

The AWGN panel is a control: with independent symbol errors the gain is
small.  The burst panel is the reason to interleave: several adjacent QPSK
symbols are erased, so an un-interleaved Hamming codeword collects multiple
errors while an interleaved codeword gets at most one bit from each burst.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from modcode.channel.awgn import AWGNChannel
from modcode.channel.burst import BurstChannel
from modcode.fec.hamming import Hamming74
from modcode.modulation.qpsk import QPSK
from modcode.pipeline.chain import Link, run_link
from modcode.pipeline.interleaver import PermutationInterleaver
from modcode.viz.ber_plot import plot_ber

plt.switch_backend("Agg")

OUTPUT = Path(__file__).resolve().parent / "figures" / "interleaver.png"
EBN0_DB = np.asarray([4.0, 6.0, 8.0, 10.0])
N_BITS = 100_000
FLOOR = 1e-6


def _sweep(channel: object, use_interleaver: bool) -> np.ndarray:
    values = []
    for index, ebn0_db in enumerate(EBN0_DB):
        link = Link(
            QPSK(),
            channel,
            Hamming74(),
            decoding="soft",
            interleaver=(
                PermutationInterleaver(seed=17) if use_interleaver else None
            ),
        )
        result = run_link(
            link,
            n_bits=N_BITS,
            ebn0_db=float(ebn0_db),
            rng=np.random.default_rng(2026 + index),
        )
        values.append(result.ber)
    return np.asarray(values)


def main() -> None:
    awgn_plain = _sweep(AWGNChannel(), False)
    awgn_interleaved = _sweep(AWGNChannel(), True)
    burst_plain = _sweep(
        BurstChannel(burst_length=8, burst_probability=0.03), False
    )
    burst_interleaved = _sweep(
        BurstChannel(burst_length=8, burst_probability=0.03), True
    )

    print("Eb/N0   AWGN plain  AWGN IL   burst plain  burst IL")
    for index, ebn0_db in enumerate(EBN0_DB):
        print(
            f"{ebn0_db:5.1f}   {awgn_plain[index]:10.5f}  "
            f"{awgn_interleaved[index]:8.5f}  "
            f"{burst_plain[index]:11.5f}  "
            f"{burst_interleaved[index]:8.5f}"
        )

    figure, axes = plt.subplots(1, 2, figsize=(11.0, 4.2))
    plot_ber(
        EBN0_DB,
        {
            "AWGN, no IL": np.maximum(awgn_plain, FLOOR),
            "AWGN, interleaved": np.maximum(awgn_interleaved, FLOOR),
        },
        axes[0],
    )
    axes[0].set_title("Hamming(7,4) + QPSK, AWGN")

    plot_ber(
        EBN0_DB,
        {
            "bursts, no IL": np.maximum(burst_plain, FLOOR),
            "bursts, interleaved": np.maximum(burst_interleaved, FLOOR),
        },
        axes[1],
    )
    axes[1].set_title("Hamming(7,4) + QPSK, 8-symbol bursts")

    for axis in axes:
        axis.set_ylim(FLOOR / 2.0, 1.0)

    figure.suptitle("Bit interleaving spreads a burst across codewords")
    figure.tight_layout()
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(OUTPUT, dpi=140)
    plt.close(figure)
    print(f"figure -> {OUTPUT}")


if __name__ == "__main__":
    main()
