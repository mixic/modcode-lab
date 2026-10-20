# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Compare BPSK, QPSK, and 16-QAM over AWGN, with theory overlays."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from modcode.metrics.ber import ber_sweep
from modcode.metrics.theory import theoretical_ber
from modcode.modulation.bpsk import BPSK
from modcode.modulation.qam import QAM
from modcode.modulation.qpsk import QPSK
from modcode.viz.ber_plot import plot_ber

plt.switch_backend("Agg")

OUTPUT = Path(__file__).resolve().parent / "figures" / "ber_sweep.png"


def main() -> None:
    modulations = {"BPSK": BPSK(), "QPSK": QPSK(), "16-QAM": QAM(16)}
    ebn0_db = np.arange(0.0, 11.0, 2.0)
    results = ber_sweep(
        modulations,
        ebn0_db,
        n_bits=200_000,
        seed=2026,
    )
    theory = {
        name: theoretical_ber(modulation, ebn0_db)
        for name, modulation in modulations.items()
    }

    print("Eb/N0   scheme    simulated   theory")
    for name in modulations:
        for point, measured, expected in zip(
            ebn0_db, results[name], theory[name], strict=True
        ):
            print(f"{point:5.1f}   {name:8s}  {measured:9.5f}  {expected:9.5f}")

    figure, _ = plot_ber(ebn0_db, results, theory=theory)
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    figure.tight_layout()
    figure.savefig(OUTPUT, dpi=140)
    plt.close(figure)
    print(f"figure -> {OUTPUT}")


if __name__ == "__main__":
    main()
