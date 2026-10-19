# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Compare BPSK, QPSK, and 16-QAM over an AWGN channel."""

import matplotlib.pyplot as plt
import numpy as np

from modcode.metrics.ber import ber_sweep
from modcode.modulation.bpsk import BPSK
from modcode.modulation.qam import QAM
from modcode.modulation.qpsk import QPSK
from modcode.viz.ber_plot import plot_ber


def main() -> None:
    ebn0_db = np.arange(0.0, 11.0, 2.0)
    results = ber_sweep(
        {"BPSK": BPSK(), "QPSK": QPSK(), "16-QAM": QAM(16)},
        ebn0_db,
        n_bits=100_000,
        seed=2026,
    )
    figure, _ = plot_ber(ebn0_db, results)
    figure.savefig("ber_sweep.png", dpi=150, bbox_inches="tight")
    plt.show()


if __name__ == "__main__":
    main()
