# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Plots for the CT-IM constant-weight index code."""

from __future__ import annotations

from matplotlib.axes import Axes

from modcode.modulation.ct_im import CTIM, ConstantWeightCyclicCode


def plot_codebook(ax: Axes, code: ConstantWeightCyclicCode) -> Axes:
    """Show the seven weight-3 index patterns."""
    ax.imshow(
        code.patterns,
        cmap="Greys",
        vmin=0.0,
        vmax=1.0,
        aspect="auto",
        interpolation="nearest",
    )
    ax.set_xticks(range(code.n))
    ax.set_yticks(range(len(code)))
    ax.set_xlabel("subcarrier", fontsize=8)
    ax.set_ylabel("index label", fontsize=8)
    ax.set_title("weight-3 CT-IM codebook", fontsize=9)
    ax.tick_params(labelsize=7)
    return ax


def plot_distance_matrix(ax: Axes, code: ConstantWeightCyclicCode) -> Axes:
    """Show the pairwise Hamming distances of the index codebook."""
    matrix = code.distance_matrix()
    ax.imshow(matrix, cmap="Blues", aspect="equal", interpolation="nearest")
    for row in range(matrix.shape[0]):
        for column in range(matrix.shape[1]):
            ax.text(
                column,
                row,
                str(int(matrix[row, column])),
                ha="center",
                va="center",
                fontsize=7,
                color="white" if matrix[row, column] > 3 else "black",
            )
    ax.set_xticks(range(len(code)))
    ax.set_yticks(range(len(code)))
    ax.set_xlabel("label", fontsize=8)
    ax.set_ylabel("label", fontsize=8)
    ax.set_title("pairwise Hamming distance", fontsize=9)
    ax.tick_params(labelsize=7)
    return ax


def plot_rate_budget(ax: Axes, ct: CTIM, *, classical_bits: float) -> Axes:
    """Compare CT-IM and classical OFDM-IM bits per subcarrier."""
    labels = ["CT-IM", "OFDM-IM"]
    values = [ct.spectral_efficiency, classical_bits / ct.n]
    bars = ax.bar(labels, values, color=["tab:blue", "tab:gray"], alpha=0.85)
    for bar, value in zip(bars, values, strict=True):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            value,
            f"{value:.3f}",
            ha="center",
            va="bottom",
            fontsize=8,
        )
    ax.set_ylabel("bits per subcarrier", fontsize=8)
    ax.set_title("rate cost of the index code", fontsize=9)
    ax.tick_params(labelsize=7)
    ax.grid(axis="y", alpha=0.3)
    return ax


__all__ = ["plot_codebook", "plot_distance_matrix", "plot_rate_budget"]
