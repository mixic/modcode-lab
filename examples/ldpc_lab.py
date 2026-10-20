#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Compare a short LDPC code with the convolutional code.

Run from the repository root:

    python examples/ldpc_lab.py

The LDPC code is the pedagogical quasi-cyclic (3,6) construction in
``modcode.fec.ldpc``, decoded with normalized min-sum.  It is not a 5G LDPC
implementation and the curve is a finite-length result, not an asymptotic
threshold claim.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from modcode.channel.awgn import AWGNChannel
from modcode.fec.convolutional import ConvolutionalCode
from modcode.fec.ldpc import LDPC36
from modcode.metrics.theory import theoretical_ber
from modcode.modulation.bpsk import BPSK
from modcode.pipeline.chain import Link, run_link
from modcode.viz.ber_plot import plot_ber

plt.switch_backend("Agg")

OUTPUT = Path(__file__).resolve().parent / "figures" / "ldpc.png"
EBN0_DB = np.asarray([0.0, 1.0, 2.0, 3.0, 4.0, 5.0])
N_BITS = 100_000
FLOOR = 1e-6


def _sweep(fec, decoding: str) -> np.ndarray:
    values = []
    for index, ebn0_db in enumerate(EBN0_DB):
        result = run_link(
            Link(BPSK(), AWGNChannel(), fec, decoding=decoding),
            n_bits=N_BITS,
            ebn0_db=float(ebn0_db),
            rng=np.random.default_rng(2026 + index),
        )
        values.append(result.ber)
    return np.asarray(values)


def main() -> None:
    convolutional = ConvolutionalCode(block_bits=64)
    ldpc = LDPC36()
    conv_ber = _sweep(convolutional, "soft")
    ldpc_ber = _sweep(ldpc, "soft")

    print(
        f"rates: convolutional {convolutional.rate:.4f}, "
        f"LDPC {ldpc.rate:.4f} (N={ldpc.n}, k={ldpc.k})"
    )
    print("Eb/N0   convolutional   LDPC")
    for index, ebn0_db in enumerate(EBN0_DB):
        print(
            f"{ebn0_db:5.1f}   {conv_ber[index]:13.5f}   "
            f"{ldpc_ber[index]:.5f}"
        )

    figure, axis = plt.subplots(figsize=(7.2, 4.8))
    plot_ber(
        EBN0_DB,
        {
            "convolutional K=3, soft": np.maximum(conv_ber, FLOOR),
            "LDPC (3,6), min-sum": np.maximum(ldpc_ber, FLOOR),
        },
        axis,
        theory={"uncoded BPSK": theoretical_ber(BPSK(), EBN0_DB)},
    )
    axis.set_ylim(FLOOR / 2.0, 1.0)
    axis.set_title("Short LDPC versus convolutional coding")
    figure.tight_layout()
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(OUTPUT, dpi=140)
    plt.close(figure)
    print(f"figure -> {OUTPUT}")


if __name__ == "__main__":
    main()
