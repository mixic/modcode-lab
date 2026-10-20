#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""RM(1,3), a tau-invariant code, versus cyclic Hamming(7,4).

Run from the repository root:

    python examples/tau_invariant_lab.py

The point is structural, not a gain claim.  RM(1,3) is invariant under the
coordinate flip i -> i XOR 1, which is not a cyclic shift.  With ML decoding
the two codes have similar soft-decision BER because both are short codes
with one-error-correcting distance; the difference is the symmetry class.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from modcode.channel.awgn import AWGNChannel
from modcode.fec.hamming import Hamming74
from modcode.fec.tau_invariant import rm13
from modcode.metrics.theory import theoretical_ber
from modcode.modulation.bpsk import BPSK
from modcode.pipeline.chain import Link, run_link
from modcode.viz.ber_plot import plot_ber

plt.switch_backend("Agg")

OUTPUT = Path(__file__).resolve().parent / "figures" / "tau_invariant.png"
EBN0_DB = np.asarray([0.0, 2.0, 4.0, 6.0, 8.0])
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
    hamming = Hamming74()
    rm = rm13()
    hamming_hard = _sweep(hamming, "hard")
    hamming_soft = _sweep(hamming, "soft")
    rm_hard = _sweep(rm, "hard")
    rm_soft = _sweep(rm, "soft")

    print(
        f"Hamming(7,4): rate {hamming.rate:.4f}, d_min 3\n"
        f"RM(1,3):      rate {rm.rate:.4f}, d_min {rm.min_distance()}\n"
    )
    print("Eb/N0   Ham-hard  Ham-soft  RM-hard  RM-soft")
    for index, ebn0_db in enumerate(EBN0_DB):
        print(
            f"{ebn0_db:5.1f}   {hamming_hard[index]:8.5f}  "
            f"{hamming_soft[index]:8.5f}  {rm_hard[index]:7.5f}  "
            f"{rm_soft[index]:7.5f}"
        )

    figure, axis = plt.subplots(figsize=(7.2, 4.8))
    plot_ber(
        EBN0_DB,
        {
            "Hamming hard (cyclic)": np.maximum(hamming_hard, FLOOR),
            "Hamming soft (cyclic)": np.maximum(hamming_soft, FLOOR),
            "RM(1,3) hard (tau-invariant)": np.maximum(rm_hard, FLOOR),
            "RM(1,3) soft (tau-invariant)": np.maximum(rm_soft, FLOOR),
        },
        axis,
        theory={"uncoded BPSK": theoretical_ber(BPSK(), EBN0_DB)},
    )
    axis.set_ylim(FLOOR / 2.0, 1.0)
    axis.set_title("Cyclic Hamming vs tau-invariant RM(1,3)")
    figure.tight_layout()
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(OUTPUT, dpi=140)
    plt.close(figure)
    print(f"figure -> {OUTPUT}")


if __name__ == "__main__":
    main()
