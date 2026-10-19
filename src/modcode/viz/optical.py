# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Matplotlib visualizations for optical intensity formats and experiments."""

from collections.abc import Mapping

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from numpy.typing import NDArray

from modcode.optical.channel import IMDDChannel
from modcode.optical.modulation import IntensityPAM


def plot_intensity_constellations(
    modulations: Mapping[str, IntensityPAM],
    ax: Axes | None = None,
) -> tuple[Figure, Axes]:
    """Compare optical intensity levels as labeled points on a power axis."""
    if not modulations:
        raise ValueError("at least one modulation is required")
    if ax is None:
        figure, axes = plt.subplots()
    else:
        axes = ax
        figure = axes.figure
    for index, (name, modulation) in enumerate(modulations.items()):
        levels = modulation.constellation
        vertical = np.full(levels.size, index, dtype=np.float64)
        axes.scatter(levels, vertical, s=60, label=name)
        for level_index, (level, y) in enumerate(zip(levels, vertical, strict=True)):
            label = format(level_index, f"0{modulation.bits_per_symbol}b")
            axes.annotate(
                label,
                (level, y),
                xytext=(0, 8),
                textcoords="offset points",
                ha="center",
            )
    axes.set_yticks(range(len(modulations)), list(modulations))
    axes.set_xlabel("Normalized optical intensity")
    axes.set_title("Optical intensity constellations")
    axes.grid(True, axis="x", alpha=0.3)
    axes.legend()
    return figure, axes


def plot_received_levels(
    modulation: IntensityPAM,
    channel: IMDDChannel,
    n_symbols: int = 200,
    seed: int = 0,
    ax: Axes | None = None,
) -> tuple[Figure, Axes]:
    """Show ideal and impaired symbol levels for a reproducible random sequence."""
    if n_symbols <= 0:
        raise ValueError("n_symbols must be positive")
    generator = np.random.default_rng(seed)
    bits = generator.integers(
        0,
        2,
        n_symbols * modulation.bits_per_symbol,
        dtype=np.uint8,
    )
    transmitted = modulation.modulate(bits)
    received = channel.transmit(transmitted, rng=generator)
    if ax is None:
        figure, axes = plt.subplots()
    else:
        axes = ax
        figure = axes.figure
    indexes = np.arange(transmitted.size)
    axes.plot(indexes, transmitted, ".", label="Transmitted")
    axes.plot(indexes, received, ".", alpha=0.7, label="Received")
    axes.set_xlabel("Symbol index")
    axes.set_ylabel("Normalized optical intensity")
    axes.set_title(f"Received levels: {type(modulation).__name__}")
    axes.grid(True, alpha=0.3)
    axes.legend()
    return figure, axes


def plot_robustness(
    results: Mapping[str, Mapping[str, NDArray[np.floating]]],
    ax: Axes | None = None,
) -> tuple[Figure, Axes]:
    """Plot BER comparison across each configured channel condition."""
    if not results:
        raise ValueError("at least one result is required")
    if ax is None:
        figure, axes = plt.subplots()
    else:
        axes = ax
        figure = axes.figure
    for modulation_name, conditions in results.items():
        for condition_name, values in conditions.items():
            ber = np.asarray(values, dtype=np.float64)
            if ber.ndim != 1 or not np.all(np.isfinite(ber)):
                raise ValueError("BER values must be finite one-dimensional arrays")
            if np.any((ber < 0) | (ber > 1)):
                raise ValueError("BER values must be between 0 and 1")
            axes.plot(
                np.arange(ber.size),
                np.maximum(ber, np.finfo(np.float64).tiny),
                marker="o",
                label=f"{modulation_name} — {condition_name}",
            )
    axes.set_xlabel("Configured channel case")
    axes.set_ylabel("Bit error rate")
    axes.set_yscale("log")
    axes.set_title("Optical IM/DD robustness comparison")
    axes.grid(True, which="both", alpha=0.3)
    axes.legend()
    return figure, axes
