# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""BER curve plotting."""

from collections.abc import Mapping, Sequence

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from numpy.typing import NDArray


def plot_ber(
    ebn0_db: NDArray[np.floating] | Sequence[float],
    results: Mapping[str, NDArray[np.floating] | Sequence[float]],
    ax: Axes | None = None,
    *,
    theory: Mapping[str, NDArray[np.floating] | Sequence[float]] | None = None,
) -> tuple[Figure, Axes]:
    """Plot named BER curves on a semilogarithmic y-axis.

    ``theory`` optionally overlays closed-form references as dashed gray
    curves.  Simulation and theory are separate arguments so a caller cannot
    accidentally label a measured curve as exact.
    """
    points = np.asarray(ebn0_db, dtype=np.float64)
    if points.ndim != 1 or points.size == 0:
        raise ValueError("ebn0_db must be a non-empty one-dimensional sequence")
    if not results:
        raise ValueError("at least one BER result is required")
    if ax is None:
        figure, axes = plt.subplots()
    else:
        axes = ax
        figure = axes.figure
    for name, values in results.items():
        ber = np.asarray(values, dtype=np.float64)
        if ber.shape != points.shape:
            raise ValueError(f"BER values for {name!r} must match ebn0_db")
        if not np.all(np.isfinite(ber)) or np.any((ber < 0) | (ber > 1)):
            raise ValueError(f"BER values for {name!r} must be between 0 and 1")
        axes.semilogy(
            points,
            np.maximum(ber, np.finfo(np.float64).tiny),
            marker="o",
            label=name,
        )
    if theory is not None:
        for name, values in theory.items():
            ber = np.asarray(values, dtype=np.float64)
            if ber.shape != points.shape:
                raise ValueError(
                    f"theoretical values for {name!r} must match ebn0_db"
                )
            if not np.all(np.isfinite(ber)) or np.any((ber < 0) | (ber > 1)):
                raise ValueError(
                    f"theoretical values for {name!r} must be between 0 and 1"
                )
            axes.semilogy(
                points,
                np.maximum(ber, np.finfo(np.float64).tiny),
                linestyle="--",
                color="0.45",
                linewidth=1.2,
                label=f"{name} theory",
            )
    axes.set_xlabel("Eb/N0 (dB)")
    axes.set_ylabel("Bit error rate")
    axes.set_title("BER vs. Eb/N0")
    axes.grid(True, which="both", alpha=0.3)
    axes.legend()
    return figure, axes
