# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Visualisations for quasi-cyclic code-based cryptography."""

from __future__ import annotations

import numpy as np
from matplotlib.axes import Axes
from matplotlib.colors import ListedColormap
from numpy.typing import ArrayLike

from modcode.pqc import QCMDPC


def plot_parity_check(
    ax: Axes,
    code: QCMDPC,
    *,
    title: str = "QC parity check",
) -> Axes:
    """Show ``H = [H0 | H1]`` with the circulant-block boundary."""
    matrix = code.parity_check
    ax.imshow(
        matrix,
        cmap="Greys",
        vmin=0.0,
        vmax=1.0,
        aspect="equal",
        interpolation="nearest",
    )
    ax.axvline(code.block_length - 0.5, color="tab:red", linewidth=1.5)
    ax.text(
        (code.block_length - 1) / 2,
        -1.5,
        "H0",
        ha="center",
        fontsize=9,
        color="tab:red",
    )
    ax.text(
        code.block_length + (code.block_length - 1) / 2,
        -1.5,
        "H1",
        ha="center",
        fontsize=9,
        color="tab:red",
    )
    ax.set_title(f"{title}: {matrix.shape[0]} x {matrix.shape[1]}", fontsize=9)
    ax.set_xlabel("code bit", fontsize=8)
    ax.set_ylabel("check", fontsize=8)
    ax.tick_params(labelsize=7)
    return ax


def plot_circulant_block(
    ax: Axes,
    block: ArrayLike,
    *,
    title: str = "Circulant block",
) -> Axes:
    """Show one circulant block as a shifted copy of its first row."""
    matrix = np.asarray(block, dtype=np.uint8)
    ax.imshow(
        matrix,
        cmap="Greys",
        vmin=0.0,
        vmax=1.0,
        aspect="equal",
        interpolation="nearest",
    )
    ax.set_title(f"{title}: {matrix.shape[0]} rows", fontsize=9)
    ax.set_xlabel("column", fontsize=8)
    ax.set_ylabel("row", fontsize=8)
    ax.tick_params(labelsize=7)
    return ax


def plot_key_scaling(
    ax: Axes,
    block_lengths: ArrayLike,
    dense_bits: ArrayLike,
    qc_bits: ArrayLike,
    *,
    title: str = "Key-size scaling",
) -> Axes:
    """Compare dense and quasi-cyclic public-key representations."""
    sizes = np.asarray(block_lengths, dtype=float)
    dense = np.asarray(dense_bits, dtype=float)
    qc = np.asarray(qc_bits, dtype=float)
    ax.loglog(sizes, dense, marker="o", label="dense Q matrix")
    ax.loglog(sizes, qc, marker="s", label="circulant first row")
    ax.set_xlabel("block length n", fontsize=8)
    ax.set_ylabel("public-key bits", fontsize=8)
    ax.set_title(title, fontsize=9)
    ax.tick_params(labelsize=7)
    ax.grid(which="both", alpha=0.3)
    ax.legend(fontsize=7)
    return ax


def plot_bit_vector(
    ax: Axes,
    bits: ArrayLike,
    *,
    title: str,
    color: str = "tab:blue",
) -> Axes:
    """Show a binary vector as a single-row image."""
    values = np.asarray(bits, dtype=np.uint8)
    ax.imshow(
        values.reshape(1, -1),
        cmap=ListedColormap(["white", color]),
        aspect="auto",
        interpolation="nearest",
    )
    ax.set_yticks([])
    ax.set_xlim(-0.5, values.size - 0.5)
    ax.set_xlabel("bit position", fontsize=8)
    ax.set_title(title, fontsize=9)
    ax.tick_params(labelsize=7)
    return ax


def plot_syndrome_history(
    ax: Axes,
    history: ArrayLike,
    *,
    title: str = "Bit-flipping progress",
) -> Axes:
    """Plot the syndrome weight against decoder iteration."""
    values = np.asarray(history, dtype=float)
    ax.plot(np.arange(values.size), values, marker="o", markersize=3)
    ax.set_xlabel("iteration", fontsize=8)
    ax.set_ylabel("syndrome weight", fontsize=8)
    ax.set_title(title, fontsize=9)
    ax.tick_params(labelsize=7)
    ax.grid(alpha=0.3)
    return ax


__all__ = [
    "plot_bit_vector",
    "plot_circulant_block",
    "plot_key_scaling",
    "plot_parity_check",
    "plot_syndrome_history",
]
