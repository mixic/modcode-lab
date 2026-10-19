# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Plots for transport budgets and transduction noise."""

from __future__ import annotations

from collections.abc import Mapping, Sequence

import numpy as np
from matplotlib.axes import Axes
from numpy.typing import ArrayLike, NDArray

from modcode.qfc.link import TransportPath


def plot_transport_budget(
    ax: Axes,
    paths: Sequence[TransportPath],
    *,
    title: str | None = None,
) -> Axes:
    """Horizontal bar chart of net loss, on a log scale."""
    labels = [path.label for path in paths]
    values = [max(path.loss_db, 1.0) for path in paths]
    colors = [
        "tab:green",
        "tab:blue",
        "tab:purple",
        "tab:orange",
        "tab:brown",
        "tab:red",
    ]
    bars = ax.barh(labels, values, color=colors[: len(paths)], alpha=0.85)
    ax.set_xscale("log")
    ax.set_xlim(1.0, max(values) * 20.0)
    ax.set_xlabel("net loss over the route (dB, log scale)", fontsize=8)
    ax.tick_params(labelsize=7)
    ax.grid(axis="x", alpha=0.3)
    for bar, value in zip(bars, values, strict=True):
        text = f"{value:,.0f} dB" if value >= 100.0 else f"{value:,.1f} dB"
        ax.text(
            value * 1.2,
            bar.get_y() + bar.get_height() / 2.0,
            text,
            va="center",
            fontsize=7,
        )
    if title:
        ax.set_title(title, fontsize=10)
    return ax


def plot_noise_budget(
    ax: Axes,
    efficiency: ArrayLike,
    calibrated: ArrayLike,
    uncalibrated: ArrayLike,
    *,
    budget: float = 0.1,
    title: str | None = None,
) -> Axes:
    """Input-referred excess noise against conversion efficiency."""
    eta = np.asarray(efficiency, dtype=np.float64)
    ax.loglog(
        eta,
        np.maximum(calibrated, 1e-10),
        color="tab:green",
        linewidth=1.8,
        label="calibrated loss",
    )
    ax.loglog(
        eta,
        np.maximum(uncalibrated, 1e-10),
        color="tab:red",
        linewidth=1.8,
        label="uncalibrated loss",
    )
    ax.axhline(budget, color="black", linestyle="--", linewidth=1.1)
    ax.text(
        eta[0] * 1.2,
        budget * 1.4,
        f"{budget:g} SNU QKD budget",
        fontsize=7,
        color="black",
    )
    ax.set_xlabel("conversion efficiency", fontsize=8)
    ax.set_ylabel("input-referred excess noise (SNU)", fontsize=8)
    ax.set_ylim(1e-6, 1e4)
    ax.tick_params(labelsize=7)
    ax.grid(which="both", alpha=0.3)
    ax.legend(fontsize=7)
    if title:
        ax.set_title(title, fontsize=10)
    return ax


def plot_thermal_occupation(
    ax: Axes,
    temperature: ArrayLike,
    microwave: ArrayLike,
    optical: ArrayLike,
    *,
    budget: float = 0.1,
    title: str | None = None,
) -> Axes:
    """Thermal photons per mode against bath temperature."""
    temps = np.asarray(temperature, dtype=np.float64)
    ax.loglog(
        temps,
        np.maximum(microwave, 1e-30),
        color="tab:blue",
        linewidth=1.8,
        label="10 GHz microwave",
    )
    ax.loglog(
        temps,
        np.maximum(optical, 1e-30),
        color="tab:orange",
        linewidth=1.8,
        label="1550 nm optical",
    )
    ax.axhline(budget, color="tab:red", linestyle="--", linewidth=1.1)
    ax.text(
        temps[0] * 1.3,
        budget * 1.5,
        f"{budget:g} SNU budget",
        fontsize=7,
        color="tab:red",
    )
    ax.set_xlabel("bath temperature (K)", fontsize=8)
    ax.set_ylabel("thermal photons per mode", fontsize=8)
    ax.tick_params(labelsize=7)
    ax.grid(which="both", alpha=0.3)
    ax.legend(fontsize=7)
    if title:
        ax.set_title(title, fontsize=10)
    return ax


def plot_wigner(
    ax: Axes,
    q: ArrayLike,
    p: ArrayLike,
    wigner: NDArray[np.float64],
    *,
    title: str | None = None,
) -> Axes:
    """Contour plot of a Gaussian Wigner function."""
    Q, P = np.meshgrid(np.asarray(q, dtype=np.float64), np.asarray(p, dtype=np.float64))
    peak = float(np.max(wigner))
    levels = np.linspace(0.08, 0.9, 5) * peak
    ax.contour(Q, P, wigner, levels=levels, colors="tab:blue", linewidths=1.0)
    ax.axhline(0.0, color="0.6", linewidth=0.6)
    ax.axvline(0.0, color="0.6", linewidth=0.6)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("q (SNU)", fontsize=8)
    ax.set_ylabel("p (SNU)", fontsize=8)
    ax.tick_params(labelsize=7)
    if title:
        ax.set_title(title, fontsize=9)
    return ax


def plot_squeezing_survival(
    ax: Axes,
    efficiency: ArrayLike,
    curves: Mapping[str, ArrayLike],
    *,
    title: str | None = None,
) -> Axes:
    """Output squeezing against conversion efficiency for several noise floors."""
    eta = np.asarray(efficiency, dtype=np.float64)
    for name, values in curves.items():
        ax.plot(eta, np.asarray(values, dtype=np.float64), linewidth=1.6, label=name)
    ax.axhline(0.0, color="0.35", linestyle=":", linewidth=1.0)
    ax.set_xlabel("conversion efficiency", fontsize=8)
    ax.set_ylabel("output squeezing (dB)", fontsize=8)
    ax.tick_params(labelsize=7)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7)
    if title:
        ax.set_title(title, fontsize=10)
    return ax


__all__ = [
    "plot_noise_budget",
    "plot_squeezing_survival",
    "plot_thermal_occupation",
    "plot_transport_budget",
    "plot_wigner",
]
