# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Plots for the I/Q phase-imbalance and tau-involution experiment."""

from __future__ import annotations

from collections.abc import Mapping

import numpy as np
from matplotlib.axes import Axes
from numpy.typing import ArrayLike

from modcode.channel.phase_imbalance import PhaseImbalanceChannel
from modcode.modulation.base import Modulation


def plot_phase_imbalance_constellation(
    ax: Axes,
    modulation: Modulation,
    phi: float,
    *,
    title: str | None = None,
) -> Axes:
    """Show ideal and phase-impaired points with the I/Q decision axes."""
    points = np.asarray(
        list(modulation.constellation.values()), dtype=np.complex128
    )
    labels = list(modulation.constellation.keys())
    impaired = PhaseImbalanceChannel(phi).impair(points)

    ax.scatter(
        points.real,
        points.imag,
        s=95,
        facecolors="none",
        edgecolors="tab:blue",
        linewidths=1.3,
        label="ideal",
    )
    ax.scatter(
        impaired.real,
        impaired.imag,
        s=55,
        marker="x",
        color="tab:orange",
        linewidths=1.6,
        label=f"impaired, phi={phi:.3f}",
    )
    for point, shifted in zip(points, impaired, strict=True):
        if abs(shifted - point) > 1e-12:
            ax.annotate(
                "",
                xy=(shifted.real, shifted.imag),
                xytext=(point.real, point.imag),
                arrowprops={
                    "arrowstyle": "->",
                    "color": "tab:orange",
                    "linewidth": 1.1,
                    "alpha": 0.75,
                },
            )
    for label, point in zip(labels, points, strict=True):
        text = "".join(str(int(bit)) for bit in label)
        ax.annotate(
            text,
            (point.real, point.imag),
            textcoords="offset points",
            xytext=(6, 6),
            fontsize=7,
            color="tab:blue",
        )

    ax.axhline(0.0, color="0.55", linewidth=0.8, linestyle="--")
    ax.axvline(0.0, color="0.55", linewidth=0.8, linestyle="--")
    ax.set_xlabel("in-phase")
    ax.set_ylabel("quadrature")
    ax.set_title(
        title or "QPSK: Q-axis reflection at phi = pi", fontsize=9
    )
    ax.set_aspect("equal", adjustable="datalim")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=7, loc="lower left")
    ax.tick_params(labelsize=7)
    return ax


def plot_phase_imbalance_ber(
    ax: Axes,
    ebn0_db: ArrayLike,
    curves: Mapping[str, ArrayLike],
    *,
    title: str = "BER under the phi = pi reflection",
) -> Axes:
    """Plot BER curves for several detectors at a fixed phase mismatch."""
    x = np.asarray(ebn0_db, dtype=np.float64)
    if x.ndim != 1 or x.size == 0:
        raise ValueError("ebn0_db must be a non-empty 1-D sequence")
    for name, values in curves.items():
        y = np.asarray(values, dtype=np.float64)
        if y.shape != x.shape:
            raise ValueError(f"curve {name!r} must match ebn0_db")
        ax.semilogy(x, y, marker="o", markersize=3.5, linewidth=1.4, label=name)
    ax.set_xlabel("Eb/N0 (dB)")
    ax.set_ylabel("bit error rate")
    ax.set_title(title, fontsize=9)
    ax.set_ylim(1e-5, 1.0)
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=7)
    ax.tick_params(labelsize=7)
    return ax


def plot_channel_determinant(
    ax: Axes,
    phis: ArrayLike,
    *,
    title: str = "Channel geometry",
) -> Axes:
    """Plot ``det R(phi) = cos(phi)`` against the phase mismatch."""
    values = np.asarray(phis, dtype=np.float64)
    if values.ndim != 1 or values.size == 0:
        raise ValueError("phis must be a non-empty 1-D sequence")
    determinant = np.cos(values)
    ax.plot(values, determinant, color="tab:purple", linewidth=1.6)
    ax.axhline(0.0, color="0.5", linewidth=0.8)
    ax.axvline(np.pi / 2.0, color="tab:red", linestyle=":", linewidth=1.2)
    ax.axvline(np.pi, color="tab:green", linestyle="--", linewidth=1.2)
    ax.text(
        np.pi / 2.0,
        0.08,
        "rank collapse",
        ha="center",
        va="bottom",
        rotation=90,
        fontsize=7,
        color="tab:red",
    )
    ax.text(
        np.pi,
        -0.85,
        "reflection",
        ha="center",
        va="bottom",
        fontsize=7,
        color="tab:green",
    )
    ax.set_xticks([0.0, np.pi / 2.0, np.pi, 1.5 * np.pi, 2.0 * np.pi])
    ax.set_xticklabels(["0", "pi/2", "pi", "3pi/2", "2pi"])
    ax.set_xlabel("phase mismatch phi (rad)")
    ax.set_ylabel("det R(phi) = cos(phi)")
    ax.set_title(title, fontsize=9)
    ax.set_ylim(-1.15, 1.15)
    ax.grid(alpha=0.25)
    ax.tick_params(labelsize=7)
    return ax


def plot_tau_llr(
    ax: Axes,
    llr: ArrayLike,
    *,
    title: str = "Blind tau evidence",
) -> Axes:
    """Plot the per-symbol log-likelihood ratio for the two hypotheses."""
    values = np.asarray(llr, dtype=np.float64)
    if values.ndim != 1 or values.size == 0:
        raise ValueError("llr must be a non-empty 1-D sequence")
    if not np.all(np.isfinite(values)):
        raise ValueError("llr must be finite")
    ax.plot(
        np.arange(values.size),
        values,
        ".",
        markersize=3.0,
        alpha=0.55,
        color="tab:blue",
    )
    ax.axhline(0.0, color="tab:red", linewidth=1.0)
    limit = max(1e-3, 1.25 * float(np.max(np.abs(values))))
    ax.set_ylim(-limit, limit)
    ax.set_xlabel("symbol index")
    ax.set_ylabel("log-likelihood ratio")
    ax.set_title(title, fontsize=9)
    if np.max(np.abs(values)) < 1e-9:
        ax.text(
            0.02,
            0.95,
            "LLR = 0 for every symbol: no threshold separates H0 from H1",
            transform=ax.transAxes,
            va="top",
            fontsize=7,
            color="tab:red",
        )
    ax.grid(alpha=0.25)
    ax.tick_params(labelsize=7)
    return ax


__all__ = [
    "plot_channel_determinant",
    "plot_phase_imbalance_ber",
    "plot_phase_imbalance_constellation",
    "plot_tau_llr",
]
