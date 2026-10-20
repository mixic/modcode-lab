#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Visualise the CT-IM constant-weight index code.

Run from the repository root:

    python examples/ct_im_lab.py

This is the algebraic/index-code demo. BER and PAPR sweeps are future work;
the important object here is the seven-word weight-3 Hamming orbit.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from modcode.modulation.ct_im import CTIM, ConstantWeightCyclicCode
from modcode.viz.ct_im import (
    plot_codebook,
    plot_distance_matrix,
    plot_rate_budget,
)

plt.switch_backend("Agg")

OUTPUT_DIR = Path(__file__).resolve().parent / "figures"
RULE = "=" * 78


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    code = ConstantWeightCyclicCode()
    ct = CTIM(M=4)
    classical_bits = 5 + 3 * 2  # 5 index bits + 3 QPSK symbols

    figure, axes = plt.subplots(1, 3, figsize=(16.0, 4.6))
    figure.suptitle(
        "CT-IM constant-weight index code: seven weight-3 Hamming words, d=4",
        fontsize=12,
    )
    plot_codebook(axes[0], code)
    plot_distance_matrix(axes[1], code)
    plot_rate_budget(axes[2], ct, classical_bits=classical_bits)
    figure.tight_layout(rect=(0, 0, 1, 0.93))
    path = OUTPUT_DIR / "ct_im_codebook.png"
    figure.savefig(path, dpi=140)
    plt.close(figure)

    # Demonstrate that every single index-bit error is corrected.
    corrected = 0
    total = 0
    for label in range(len(code)):
        word = code.encode(label)
        for position in range(code.n):
            corrupted = word.copy()
            corrupted[position] ^= 1
            decoded, distance = code.decode(corrupted)
            corrected += int(decoded == label and distance <= 1)
            total += 1

    rng = np.random.default_rng(2026)
    bits = rng.integers(0, 2, size=ct.bits_per_block, dtype=np.uint8)
    bits[:3] = [0, 0, 1]
    vector = ct.encode(bits)
    recovered = ct.decode(vector)

    print(RULE)
    print("CT-IM CONSTANT-WEIGHT INDEX CODE")
    print(RULE)
    print(f"subcarriers n          : {code.n}")
    print(f"active subcarriers k   : {code.weight}")
    print(f"index words            : {len(code)}")
    print(f"minimum distance       : {code.minimum_distance}")
    print(f"index entropy          : log2(7) = {np.log2(7):.4f} bits")
    print("labels used            : 0..6 (label 7 reserved)")
    print(f"symbol bits (QPSK)     : {ct.symbol_bits}")
    print(f"bits per block         : {ct.bits_per_block}")
    print(f"bits per subcarrier    : {ct.spectral_efficiency:.4f}")
    print(f"classical OFDM-IM rate : {classical_bits / ct.n:.4f}")
    print(f"single-error fixes     : {corrected}/{total}")
    print(f"noiseless round trip   : {np.array_equal(recovered, bits)}")
    print()
    print("  The codebook is one cyclic orbit of size 7, not a linear")
    print("  subspace. It corrects one index bit and costs rate relative")
    print("  to the unconstrained C(7,3)=35-pattern OFDM-IM alphabet.")
    print(f"  figure -> {path}")
    print(RULE)


if __name__ == "__main__":
    main()
