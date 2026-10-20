#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Block Rayleigh fading: uncoded references and hard-decoded FEC.

Run from the repository root:

    python examples/fading_lab.py

The receiver knows the fading coefficient and equalizes before slicing.
Hard decisions are used deliberately: after zero-forcing, the noise is
``sigma^2 / |h|^2`` per symbol, so a single scalar variance would be wrong
for soft LLRs.  A per-symbol-variance soft path is a later extension.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from modcode.channel.rayleigh import BlockRayleighChannel
from modcode.fec.bch import BCH15
from modcode.fec.convolutional import ConvolutionalCode
from modcode.fec.hamming import Hamming74
from modcode.fec.reed_solomon import ReedSolomon
from modcode.metrics.theory import theoretical_rayleigh_ber
from modcode.modulation.bpsk import BPSK
from modcode.modulation.qpsk import QPSK
from modcode.pipeline.chain import Link, run_link
from modcode.pipeline.interleaver import PermutationInterleaver
from modcode.viz.ber_plot import plot_ber

plt.switch_backend("Agg")

OUTPUT = Path(__file__).resolve().parent / "figures" / "rayleigh.png"
EBN0_DB = np.asarray([0.0, 4.0, 8.0, 12.0, 16.0, 20.0])
N_BITS = 100_000
FLOOR = 1e-6


def _run(modulation, fec=None, *, interleaver=None) -> np.ndarray:
    values = []
    for index, ebn0_db in enumerate(EBN0_DB):
        result = run_link(
            Link(
                modulation,
                BlockRayleighChannel(block_length=10),
                fec,
                decoding="hard" if fec is not None else "auto",
                interleaver=interleaver,
            ),
            n_bits=N_BITS,
            ebn0_db=float(ebn0_db),
            rng=np.random.default_rng(2026 + index),
        )
        values.append(result.ber)
    return np.asarray(values)


def main() -> None:
    uncoded_bpsk = _run(BPSK())
    uncoded_qpsk = _run(QPSK())
    hamming = _run(BPSK(), Hamming74())
    hamming_il = _run(
        BPSK(), Hamming74(), interleaver=PermutationInterleaver(seed=9)
    )
    bch = _run(BPSK(), BCH15())
    rs = _run(BPSK(), ReedSolomon())
    convolutional = _run(BPSK(), ConvolutionalCode(block_bits=64))

    print("Eb/N0   BPSK     QPSK     Ham      Ham+IL   BCH      RS       Conv")
    for index, ebn0_db in enumerate(EBN0_DB):
        print(
            f"{ebn0_db:5.1f}   {uncoded_bpsk[index]:7.5f}  "
            f"{uncoded_qpsk[index]:7.5f}  {hamming[index]:7.5f}  "
            f"{hamming_il[index]:7.5f}  {bch[index]:7.5f}  "
            f"{rs[index]:7.5f}  {convolutional[index]:7.5f}"
        )

    figure, axes = plt.subplots(1, 2, figsize=(11.0, 4.2))
    plot_ber(
        EBN0_DB,
        {
            "BPSK": np.maximum(uncoded_bpsk, FLOOR),
            "QPSK": np.maximum(uncoded_qpsk, FLOOR),
        },
        axes[0],
        theory={
            "BPSK": theoretical_rayleigh_ber(BPSK(), EBN0_DB),
            "QPSK": theoretical_rayleigh_ber(QPSK(), EBN0_DB),
        },
    )
    axes[0].set_title("Uncoded block Rayleigh, known CSI")

    plot_ber(
        EBN0_DB,
        {
            "Hamming hard": np.maximum(hamming, FLOOR),
            "Hamming + interleaver": np.maximum(hamming_il, FLOOR),
            "BCH hard": np.maximum(bch, FLOOR),
            "RS hard": np.maximum(rs, FLOOR),
            "convolutional hard": np.maximum(convolutional, FLOOR),
        },
        axes[1],
    )
    axes[1].set_title("Coded BPSK over block Rayleigh")

    for axis in axes:
        axis.set_ylim(FLOOR / 2.0, 1.0)

    figure.suptitle("Block Rayleigh fading: theory and hard-decoded FEC")
    figure.tight_layout()
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(OUTPUT, dpi=140)
    plt.close(figure)
    print(f"figure -> {OUTPUT}")


if __name__ == "__main__":
    main()
