# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Visualize APSK phase orbits, shaped symbol probabilities, and a cyclic TCM."""

import matplotlib.pyplot as plt
import numpy as np

from modcode.modulation.apsk import APSK
from modcode.modulation.cyclic_tcm import TailBitingTCM8PSK
from modcode.shaping import (
    maxwell_boltzmann_probabilities,
    sample_shaped_symbols,
)
from modcode.viz.constellation import plot_constellation


def main() -> None:
    apsk = APSK((4, 12), (1.0, 2.0))
    points = np.asarray(list(apsk.constellation.values()))
    probabilities = maxwell_boltzmann_probabilities(points, shaping_parameter=0.8)
    samples = sample_shaped_symbols(
        points,
        probabilities,
        count=50_000,
        rng=np.random.default_rng(2026),
    )

    figure, axes = plt.subplots(1, 2, figsize=(12, 5))
    plot_constellation(apsk, ax=axes[0])
    axes[0].set_title("16-APSK: each ring has a cyclic phase action")
    axes[1].scatter(
        samples.real,
        samples.imag,
        s=4,
        alpha=0.12,
        label="Maxwell-Boltzmann shaped samples",
    )
    axes[1].scatter(points.real, points.imag, marker="+", color="black")
    axes[1].set_aspect("equal", adjustable="datalim")
    axes[1].set_xlabel("In-phase")
    axes[1].set_ylabel("Quadrature")
    axes[1].set_title("Shaping changes symbol frequencies, not ring geometry")
    axes[1].grid(True, alpha=0.3)

    bits = np.random.default_rng(9).integers(0, 2, 32, dtype=np.uint8)
    tcm = TailBitingTCM8PSK()
    symbols = tcm.encode(bits)
    decoded = tcm.decode(symbols)
    assert np.array_equal(decoded, bits)
    print(
        "Tail-biting 8-PSK TCM noiseless decode:",
        np.array_equal(decoded, bits),
    )
    print("APSK ring sizes:", apsk.ring_sizes, "(cyclic phase groups C4 and C12)")
    figure.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
