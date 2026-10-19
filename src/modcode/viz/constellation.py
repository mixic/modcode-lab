# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Constellation diagram plotting."""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from numpy.typing import NDArray

from modcode.modulation.base import Modulation


def plot_constellation(
    modulation: Modulation,
    received: NDArray[np.complexfloating] | None = None,
    ax: Axes | None = None,
) -> tuple[Figure, Axes]:
    """Plot a modulation's ideal points and optional received symbols."""
    if ax is None:
        figure, axes = plt.subplots()
    else:
        axes = ax
        figure = axes.figure
    points = np.asarray(list(modulation.constellation.values()), dtype=np.complex128)
    axes.scatter(points.real, points.imag, marker="o", label="Constellation")
    for label, point in modulation.constellation.items():
        axes.annotate("".join(map(str, label)), (point.real, point.imag))
    if received is not None:
        samples = np.asarray(received, dtype=np.complex128)
        if samples.ndim != 1:
            raise ValueError("received symbols must be a one-dimensional array")
        if not np.all(np.isfinite(samples)):
            raise ValueError("received symbols must contain only finite values")
        axes.scatter(
            samples.real,
            samples.imag,
            marker=".",
            alpha=0.25,
            label="Received",
        )
        axes.legend()
    axes.set_xlabel("In-phase")
    axes.set_ylabel("Quadrature")
    axes.set_title(f"{type(modulation).__name__} constellation")
    axes.set_aspect("equal", adjustable="datalim")
    axes.grid(True, alpha=0.3)
    return figure, axes
