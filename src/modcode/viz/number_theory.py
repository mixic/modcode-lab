# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Plots for residual-index densities and finite prime counts."""

from __future__ import annotations

from collections.abc import Mapping

import numpy as np
from matplotlib.axes import Axes
from numpy.typing import ArrayLike

_COLORS = {
    1: "tab:blue",
    2: "tab:orange",
    4: "tab:green",
}


def plot_residual_index_density(
    ax: Axes,
    m_values: ArrayLike,
    densities: ArrayLike,
    *,
    artin: float,
    title: str = "Density by residual index",
) -> Axes:
    """Bar chart of predicted densities, with Artin's constant marked."""
    indexes = np.asarray(m_values, dtype=int)
    values = np.asarray(densities, dtype=float)
    colors = [_COLORS.get(int(index), "tab:gray") for index in indexes]
    bars = ax.bar(indexes.astype(str), values, color=colors, alpha=0.85)
    ax.axhline(artin, color="black", linestyle="--", linewidth=1.0)
    ax.text(
        0.02,
        artin,
        "Artin A",
        transform=ax.get_yaxis_transform(),
        fontsize=7,
        va="bottom",
    )
    for bar, value in zip(bars, values, strict=True):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            value,
            f"{value:.3f}",
            ha="center",
            va="bottom",
            fontsize=7,
        )
    ax.set_xlabel("residual index m", fontsize=8)
    ax.set_ylabel("predicted density", fontsize=8)
    ax.set_title(title, fontsize=9)
    ax.tick_params(labelsize=7)
    ax.grid(axis="y", alpha=0.3)
    return ax


def plot_residual_index_convergence(
    ax: Axes,
    limits: ArrayLike,
    series: Mapping[int, ArrayLike],
    predicted: Mapping[int, float],
    *,
    title: str = "Finite-count convergence",
) -> Axes:
    """Plot observed ratios against X, with predicted densities as guides."""
    x = np.asarray(limits, dtype=float)
    for index, values in series.items():
        color = _COLORS.get(index, "tab:gray")
        ax.semilogx(
            x,
            np.asarray(values, dtype=float),
            marker="o",
            markersize=3,
            color=color,
            label=f"m={index}",
        )
        if index in predicted:
            ax.axhline(
                predicted[index],
                color=color,
                linestyle=":",
                linewidth=1.0,
                alpha=0.8,
            )
    ax.set_xlabel("prime limit X", fontsize=8)
    ax.set_ylabel("observed fraction", fontsize=8)
    ax.set_title(title, fontsize=9)
    ax.tick_params(labelsize=7)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7)
    return ax


def plot_residual_index_histogram(
    ax: Axes,
    counts: Mapping[int, int],
    *,
    total: int,
    max_index: int = 12,
    title: str = "Residual-index histogram",
) -> Axes:
    """Plot the observed frequency of each small residual index."""
    indexes = np.arange(1, max_index + 1)
    values = np.asarray(
        [counts.get(int(index), 0) / total for index in indexes],
        dtype=float,
    )
    ax.bar(indexes, values, color="tab:purple", alpha=0.85)
    ax.set_xlabel("residual index m", fontsize=8)
    ax.set_ylabel("observed fraction", fontsize=8)
    ax.set_title(title, fontsize=9)
    ax.tick_params(labelsize=7)
    ax.grid(axis="y", alpha=0.3)
    return ax


__all__ = [
    "plot_residual_index_convergence",
    "plot_residual_index_density",
    "plot_residual_index_histogram",
]
