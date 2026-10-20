# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""CRC frame and division-register plots."""

from __future__ import annotations

import numpy as np
from matplotlib.axes import Axes
from matplotlib.patches import Rectangle
from numpy.typing import ArrayLike


def _bit_row(ax: Axes, bits: np.ndarray, *, title: str) -> None:
    ax.imshow(
        bits.reshape(1, -1),
        cmap="Greys",
        vmin=0.0,
        vmax=1.0,
        aspect="auto",
        interpolation="nearest",
    )
    if bits.size <= 64:
        for index, bit in enumerate(bits):
            ax.text(
                index,
                0,
                str(int(bit)),
                ha="center",
                va="center",
                fontsize=7,
                color="white" if bit else "black",
            )
    ax.set_yticks([])
    ax.set_xlim(-0.5, bits.size - 0.5)
    ax.set_title(title, fontsize=9)
    ax.tick_params(labelsize=7)


def plot_crc_frame(
    ax: Axes,
    message: ArrayLike,
    parity: ArrayLike,
    *,
    error_positions: ArrayLike = (),
    title: str = "CRC frame",
) -> Axes:
    """Show the systematic message-plus-CRC frame, optionally marking errors."""
    message_bits = np.asarray(message, dtype=np.uint8)
    parity_bits = np.asarray(parity, dtype=np.uint8)
    frame = np.concatenate([message_bits, parity_bits])
    _bit_row(ax, frame, title=title)
    ax.axvline(message_bits.size - 0.5, color="tab:blue", linewidth=1.5)
    ax.text(
        (message_bits.size - 1) / 2,
        -0.8,
        "message",
        ha="center",
        va="bottom",
        fontsize=8,
        color="tab:blue",
    )
    ax.text(
        message_bits.size + (parity_bits.size - 1) / 2,
        -0.8,
        "CRC",
        ha="center",
        va="bottom",
        fontsize=8,
        color="tab:blue",
    )
    for position in np.atleast_1d(np.asarray(error_positions, dtype=int)):
        ax.add_patch(
            Rectangle(
                (position - 0.5, -0.5),
                1.0,
                1.0,
                fill=False,
                edgecolor="tab:red",
                linewidth=2.0,
            )
        )
    ax.set_xlabel("bit position", fontsize=8)
    return ax


def plot_crc_division(
    ax: Axes,
    trace: ArrayLike,
    *,
    title: str = "CRC division",
) -> Axes:
    """Show the CRC shift-register contents after each message bit."""
    values = np.asarray(trace, dtype=np.uint8)
    ax.imshow(
        values.T,
        cmap="Greys",
        vmin=0.0,
        vmax=1.0,
        aspect="auto",
        interpolation="nearest",
    )
    ax.set_xlabel("message step", fontsize=8)
    ax.set_ylabel("remainder bit", fontsize=8)
    ax.set_title(title, fontsize=9)
    ax.tick_params(labelsize=7)
    return ax


def plot_crc_syndromes(
    ax: Axes,
    positions: ArrayLike,
    syndromes: ArrayLike,
    *,
    title: str = "Single-bit error syndromes",
) -> Axes:
    """Bar chart of the CRC remainder for each single-bit error position."""
    x = np.asarray(positions, dtype=int)
    y = np.asarray(syndromes, dtype=int)
    ax.bar(x, y, color="tab:red", alpha=0.8)
    ax.set_xlabel("error position", fontsize=8)
    ax.set_ylabel("CRC remainder", fontsize=8)
    ax.set_title(title, fontsize=9)
    ax.tick_params(labelsize=7)
    ax.grid(axis="y", alpha=0.3)
    return ax


__all__ = ["plot_crc_division", "plot_crc_frame", "plot_crc_syndromes"]
