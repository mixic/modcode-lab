#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Four-state 8-PSK TCM versus uncoded QPSK at 2 bit/symbol.

Run from the repository root:

    python examples/tcm_lab.py

The TCM decoder is a soft Viterbi search over the four-state trellis.  The
comparison is at equal spectral efficiency: both schemes carry 2 information
bits per channel use.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from modcode.channel.awgn import AWGNChannel
from modcode.metrics.theory import theoretical_ber
from modcode.modulation.qpsk import QPSK
from modcode.modulation.tcm import TCM8PSK
from modcode.pipeline.chain import Link, run_link
from modcode.viz.ber_plot import plot_ber

plt.switch_backend("Agg")

OUTPUT = Path(__file__).resolve().parent / "figures" / "tcm.png"
EBN0_DB = np.asarray([2.0, 4.0, 6.0, 8.0, 10.0])
N_BITS = 100_000
FLOOR = 1e-6


def main() -> None:
    tcm = TCM8PSK()
    tcm_ber = []
    qpsk_ber = []
    for index, ebn0_db in enumerate(EBN0_DB):
        rng = np.random.default_rng(2026 + index)
        bits = rng.integers(0, 2, N_BITS, dtype=np.uint8)
        symbols = tcm.encode(bits)
        variance = AWGNChannel.noise_variance(ebn0_db, 2)
        received = AWGNChannel().transmit(symbols, variance, rng)
        tcm_ber.append(float(np.mean(tcm.decode(received) != bits)))
        qpsk_ber.append(
            run_link(
                Link(QPSK(), AWGNChannel()),
                n_bits=N_BITS,
                ebn0_db=float(ebn0_db),
                rng=np.random.default_rng(2026 + index),
            ).ber
        )

    print("Eb/N0   TCM 8PSK   uncoded QPSK")
    for index, ebn0_db in enumerate(EBN0_DB):
        print(
            f"{ebn0_db:5.1f}   {tcm_ber[index]:9.5f}   {qpsk_ber[index]:12.5f}"
        )

    figure, axis = plt.subplots(figsize=(7.2, 4.8))
    plot_ber(
        EBN0_DB,
        {
            "TCM 8PSK": np.maximum(tcm_ber, FLOOR),
            "uncoded QPSK": np.maximum(qpsk_ber, FLOOR),
        },
        axis,
        theory={"QPSK": theoretical_ber(QPSK(), EBN0_DB)},
    )
    axis.set_ylim(FLOOR / 2.0, 1.0)
    axis.set_title("TCM 8-PSK vs uncoded QPSK at 2 bit/symbol")
    figure.tight_layout()
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(OUTPUT, dpi=140)
    plt.close(figure)
    print(f"figure -> {OUTPUT}")


if __name__ == "__main__":
    main()
