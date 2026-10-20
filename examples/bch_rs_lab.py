#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Compare Hamming, BCH(15,7), RS(15,11), and convolutional coding.

Run from the repository root:

    python examples/bch_rs_lab.py

All curves use BPSK and AWGN.  The BCH and Hamming codes support soft ML
decoding at these lengths; RS(15,11) uses hard Berlekamp-Massey/Chien
decoding, and the convolutional code uses soft Viterbi.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from modcode.channel.awgn import AWGNChannel
from modcode.fec.bch import BCH15
from modcode.fec.convolutional import ConvolutionalCode
from modcode.fec.hamming import Hamming74
from modcode.fec.reed_solomon import ReedSolomon
from modcode.metrics.theory import theoretical_ber
from modcode.modulation.bpsk import BPSK
from modcode.pipeline.chain import Link, run_link
from modcode.viz.ber_plot import plot_ber

plt.switch_backend("Agg")

OUTPUT = Path(__file__).resolve().parent / "figures" / "bch_rs.png"
EBN0_DB = np.asarray([0.0, 2.0, 4.0, 6.0, 8.0, 10.0])
N_BITS = 100_000
FLOOR = 1e-6


def _sweep(fec, decoding: str | None) -> np.ndarray:
    values = []
    for index, ebn0_db in enumerate(EBN0_DB):
        link = Link(
            BPSK(),
            AWGNChannel(),
            fec,
            decoding=decoding if decoding is not None else "auto",
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
    hamming_hard = _sweep(Hamming74(), "hard")
    hamming_soft = _sweep(Hamming74(), "soft")
    bch_hard = _sweep(BCH15(), "hard")
    bch_soft = _sweep(BCH15(), "soft")
    rs_hard = _sweep(ReedSolomon(), "hard")
    convolutional_soft = _sweep(ConvolutionalCode(64), "soft")

    print("Eb/N0   Ham-hard Ham-soft BCH-hard BCH-soft RS-hard Conv-soft")
    for index, ebn0_db in enumerate(EBN0_DB):
        print(
            f"{ebn0_db:5.1f}   {hamming_hard[index]:8.5f} "
            f"{hamming_soft[index]:8.5f} {bch_hard[index]:8.5f} "
            f"{bch_soft[index]:8.5f} {rs_hard[index]:7.5f} "
            f"{convolutional_soft[index]:8.5f}"
        )
    print()
    print("rate: Hamming 4/7, BCH 7/15, RS 11/15, convolutional 64/132")

    figure, axis = plt.subplots(figsize=(7.2, 4.8))
    plot_ber(
        EBN0_DB,
        {
            "Hamming hard": np.maximum(hamming_hard, FLOOR),
            "Hamming soft": np.maximum(hamming_soft, FLOOR),
            "BCH hard": np.maximum(bch_hard, FLOOR),
            "BCH soft": np.maximum(bch_soft, FLOOR),
            "RS hard": np.maximum(rs_hard, FLOOR),
            "convolutional soft": np.maximum(convolutional_soft, FLOOR),
        },
        axis,
        theory={"uncoded BPSK": theoretical_ber(BPSK(), EBN0_DB)},
    )
    axis.set_ylim(FLOOR / 2.0, 1.0)
    axis.set_title("BCH, Reed-Solomon, Hamming, and convolutional codes")
    figure.tight_layout()
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(OUTPUT, dpi=140)
    plt.close(figure)
    print(f"figure -> {OUTPUT}")


if __name__ == "__main__":
    main()
