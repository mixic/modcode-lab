# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Wavefront, phase-ramp, and beam-pattern plotting for linear arrays."""

from __future__ import annotations

import numpy as np
from matplotlib.axes import Axes
from numpy.typing import ArrayLike

from modcode.phasedarray import LinearArray, squint_angle


def plot_wavefront(
    ax: Axes,
    array: LinearArray,
    x: ArrayLike,
    z: ArrayLike,
    phases: ArrayLike,
    theta_deg: float,
    *,
    component: str = "real",
    title: str | None = None,
    show_guide: bool = True,
) -> Axes:
    """Draw the synthesised field over an ``(x, z)`` cut.

    ``component="real"`` shows the sign of the real field (the fringes that
    read as wavefronts); ``component="phase"`` shows the wrapped phase.
    """
    x_values = np.asarray(x, dtype=np.float64)
    z_values = np.asarray(z, dtype=np.float64)
    field = array.field(x_values, z_values, phases)
    extent = [x_values[0], x_values[-1], z_values[0], z_values[-1]]

    if component == "real":
        data = field.real
        limit = float(np.percentile(np.abs(data), 99.0)) or 1.0
        ax.imshow(
            data,
            extent=extent,
            origin="lower",
            cmap="RdBu_r",
            vmin=-limit,
            vmax=limit,
            aspect="auto",
            interpolation="bilinear",
        )
    elif component == "phase":
        ax.imshow(
            np.angle(field),
            extent=extent,
            origin="lower",
            cmap="twilight",
            vmin=-np.pi,
            vmax=np.pi,
            aspect="auto",
            interpolation="bilinear",
        )
    else:
        raise ValueError("component must be 'real' or 'phase'")

    if show_guide and abs(theta_deg) < 89.0:
        slope = np.tan(np.radians(theta_deg))
        x_end = np.clip(
            slope * (z_values[-1] - z_values[0]), x_values[0], x_values[-1]
        )
        ax.plot(
            [0.0, x_end],
            [z_values[0], z_values[-1]],
            linestyle="--",
            linewidth=1.0,
            color="black",
            alpha=0.55,
        )

    ax.scatter(
        array.positions(),
        np.full(array.elements, z_values[0]),
        marker="^",
        s=10,
        color="black",
        clip_on=True,
    )
    if title:
        ax.set_title(title, fontsize=10)
    ax.set_xlabel("x (m)", fontsize=8)
    ax.tick_params(labelsize=7)
    return ax


def plot_phase_ramp(
    ax: Axes,
    array: LinearArray,
    theta_deg: float,
    phases: ArrayLike,
    *,
    label: str,
    title: str | None = None,
) -> Axes:
    """Plot the applied element phases against the ideal ramp."""
    index = np.arange(array.elements)
    ideal = np.degrees(np.unwrap(array.exact_phases(theta_deg)))
    applied = np.degrees(np.unwrap(np.asarray(phases, dtype=np.float64)))
    ax.plot(index, ideal, color="0.45", linestyle="--", linewidth=1.0, label="ideal")
    ax.plot(index, applied, marker="o", markersize=3.5, linewidth=1.1, label=label)
    ax.set_xlabel("element index", fontsize=8)
    ax.set_ylabel("phase (deg)", fontsize=8)
    ax.tick_params(labelsize=7)
    ax.grid(alpha=0.3)
    if title:
        ax.set_title(title, fontsize=10)
    ax.legend(fontsize=7, loc="upper left")
    return ax


def plot_pattern(
    ax: Axes,
    array: LinearArray,
    phases: ArrayLike,
    *,
    commanded_deg: float | None = None,
    label: str | None = None,
    floor_db: float = -40.0,
    mark_peak: bool = True,
    title: str | None = None,
) -> Axes:
    """Draw the normalised array-factor cut in dB."""
    grid = np.arange(-90.0, 90.0 + 0.05, 0.05)
    pattern = array.pattern_db(grid, phases)
    ax.plot(grid, pattern, linewidth=1.2, label=label)
    if commanded_deg is not None:
        ax.axvline(commanded_deg, color="0.35", linestyle=":", linewidth=1.0)
    if mark_peak:
        peak = array.peak_angle(phases)
        level = float(array.pattern_db([peak], phases)[0])
        ax.plot(peak, level, marker="v", color="tab:red", markersize=6)
        ax.annotate(
            f"{peak:.1f} deg",
            (peak, level),
            textcoords="offset points",
            xytext=(4, 6),
            fontsize=7,
            color="tab:red",
        )
    ax.set_xlim(-90.0, 90.0)
    ax.set_ylim(floor_db, 2.0)
    ax.set_xlabel("theta (deg)", fontsize=8)
    ax.set_ylabel("array factor (dB)", fontsize=8)
    ax.tick_params(labelsize=7)
    ax.grid(alpha=0.3)
    if title:
        ax.set_title(title, fontsize=10)
    if label:
        ax.legend(fontsize=7, loc="upper left")
    return ax


def plot_pattern_family(
    ax: Axes,
    array: LinearArray,
    patterns: dict[str, ArrayLike],
    *,
    commanded_deg: float | None = None,
    floor_db: float = -40.0,
    title: str | None = None,
) -> Axes:
    """Overlay several array-factor cuts for comparison."""
    grid = np.arange(-90.0, 90.0 + 0.05, 0.05)
    for name, phases in patterns.items():
        ax.plot(grid, array.pattern_db(grid, phases), linewidth=1.2, label=name)
    if commanded_deg is not None:
        ax.axvline(commanded_deg, color="0.35", linestyle=":", linewidth=1.0)
    ax.set_xlim(-90.0, 90.0)
    ax.set_ylim(floor_db, 2.0)
    ax.set_xlabel("theta (deg)", fontsize=8)
    ax.set_ylabel("array factor (dB)", fontsize=8)
    ax.tick_params(labelsize=7)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7, loc="upper left")
    if title:
        ax.set_title(title, fontsize=10)
    return ax


def plot_squint(
    ax: Axes,
    theta0_deg: ArrayLike,
    ratios: ArrayLike = (1.0,),
    *,
    title: str | None = None,
) -> Axes:
    """Plot realised pointing against drive frequency for fixed phase ramps."""
    ratio_values = np.asarray(ratios, dtype=np.float64)
    for theta in np.atleast_1d(np.asarray(theta0_deg, dtype=np.float64)):
        angles = [squint_angle(theta, ratio) for ratio in ratio_values]
        ax.plot(ratio_values, angles, linewidth=1.5, label=f"{theta:.0f} deg")
    ax.axvline(1.0, color="0.35", linestyle=":", linewidth=1.0)
    ax.set_xlabel("f / f0", fontsize=8)
    ax.set_ylabel("realised pointing (deg)", fontsize=8)
    ax.tick_params(labelsize=7)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7, loc="upper left")
    if title:
        ax.set_title(title, fontsize=10)
    return ax


__all__ = [
    "plot_pattern",
    "plot_pattern_family",
    "plot_phase_ramp",
    "plot_squint",
    "plot_wavefront",
]
