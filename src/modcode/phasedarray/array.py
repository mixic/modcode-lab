# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Linear phased-array geometry, cyclic phase resolution, and field synthesis.

The array sits on the ``x`` axis, centred on the origin, and radiates towards
``+z``. Element ``n`` is at ``x_n = (n - (N - 1) / 2) * d``. A phase ramp
``psi_n = n * dphi`` with ``dphi = 2*pi*(d/lambda)*sin(theta0)`` steers the
main beam to ``+theta0``.

A digital phase shifter cannot apply an arbitrary phase. If it is built from
the cyclic group ``C_order``, every phase must be a multiple of
``2*pi/order``. There are two distinct ways to spend that budget:

* ``element_quantized_phases`` rounds each element's phase to the nearest
  available value. The average ramp slope survives, so the beam stays near the
  commanded angle, but the aperture picks up a staircase aberration.
* ``step_quantized_phases`` quantises the inter-element step once and applies
  that step uniformly. The wavefront stays perfectly flat, so the beam keeps
  full gain, but it points at the quantised angle instead of the commanded one.

Keeping both is important: they trade pointing error against wavefront error,
and only the second one produces the "told to point at 17 degrees, points at
30" failure that is easy to misattribute to per-element rounding.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

SPEED_OF_LIGHT = 299_792_458.0
"""Speed of light in vacuum, in metres per second."""

TWO_PI = 2.0 * np.pi


def _finite_vector(values: ArrayLike, name: str) -> NDArray[np.float64]:
    array = np.asarray(values, dtype=np.float64)
    if array.ndim != 1:
        raise ValueError(f"{name} must be a one-dimensional array")
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must contain only finite values")
    return array


def _phase_order(order: int) -> int:
    if isinstance(order, bool) or not isinstance(order, (int, np.integer)):
        raise ValueError("order must be an integer")
    if order < 1:
        raise ValueError("order must be a positive integer")
    return int(order)


def quantize_phases(
    phases: ArrayLike, order: int, *, wrap: bool = False
) -> NDArray[np.float64]:
    """Round phases to the nearest multiple of ``2*pi/order``.

    Parameters
    ----------
    phases
        Applied element phases in radians.
    order
        Size of the cyclic phase group. ``order=4`` gives 0/90/180/270 degree
        shifters, ``order=8`` gives 45 degree steps, and so on.
    wrap
        If true, return the phases modulo ``2*pi``, which is what the shifter
        actually stores.
    """
    values = _finite_vector(phases, "phases")
    quantum = TWO_PI / _phase_order(order)
    quantized = np.round(values / quantum) * quantum
    return np.mod(quantized, TWO_PI) if wrap else quantized


def grating_lobe_limit(theta_deg: float) -> float:
    """Largest ``d/lambda`` that keeps grating lobes out of visible space.

    A second main lobe appears once ``d/lambda > 1 / (1 + |sin(theta0)|)``.
    The bound is the standard uniform-array no-grating-lobe condition.
    """
    return 1.0 / (1.0 + abs(np.sin(np.radians(theta_deg))))


def squint_angle(theta_deg: float, frequency_ratio: float) -> float:
    """Pointing angle when a fixed phase ramp is driven at another frequency.

    The ramp is designed at ``f0``. Because the phase is fixed while the
    electrical spacing ``d/lambda`` scales with frequency,
    ``sin(theta) = sin(theta0) / (f / f0)``: a phase-steered beam squints
    towards broadside as the drive frequency rises. Angles whose sine leaves
    ``[-1, 1]`` are reported at the horizon.
    """
    if not np.isfinite(frequency_ratio) or frequency_ratio <= 0:
        raise ValueError("frequency_ratio must be finite and positive")
    sine = np.sin(np.radians(theta_deg)) / frequency_ratio
    if abs(sine) >= 1.0:
        return float(np.sign(sine) * 90.0)
    return float(np.degrees(np.arcsin(sine)))


@dataclass(frozen=True)
class LinearArray:
    """A uniform linear array with a nominal half-wavelength design.

    Parameters
    ----------
    elements
        Number of radiating elements (at least two).
    spacing_over_lambda
        Element spacing in wavelengths. The physical spacing is frozen through
        ``frequency``; changing the drive frequency later is a squint, not a
        resize.
    frequency
        Design frequency in hertz.
    """

    elements: int = 16
    spacing_over_lambda: float = 0.5
    frequency: float = 10.0e9

    def __post_init__(self) -> None:
        if isinstance(self.elements, bool) or not isinstance(
            self.elements, (int, np.integer)
        ):
            raise ValueError("elements must be an integer")
        if self.elements < 2:
            raise ValueError("an array needs at least two elements")
        if (
            not np.isfinite(self.spacing_over_lambda)
            or not 0.0 < self.spacing_over_lambda <= 2.0
        ):
            raise ValueError("spacing_over_lambda must lie in (0, 2]")
        if not np.isfinite(self.frequency) or self.frequency <= 0.0:
            raise ValueError("frequency must be finite and positive")

    @property
    def wavelength(self) -> float:
        return SPEED_OF_LIGHT / self.frequency

    @property
    def spacing(self) -> float:
        return self.spacing_over_lambda * self.wavelength

    @property
    def aperture(self) -> float:
        return (self.elements - 1) * self.spacing

    @property
    def fresnel_distance(self) -> float:
        return self.aperture**2 / self.wavelength

    def positions(self) -> NDArray[np.float64]:
        """Element positions on the x axis, in metres, centred on zero."""
        offsets = np.arange(self.elements) - (self.elements - 1) / 2.0
        return offsets * self.spacing

    def steering_step(self, theta_deg: float) -> float:
        """Inter-element phase step in radians for a beam at ``theta_deg``."""
        return TWO_PI * self.spacing_over_lambda * np.sin(np.radians(theta_deg))

    def exact_phases(self, theta_deg: float) -> NDArray[np.float64]:
        """The ideal continuous phase ramp for ``theta_deg``, in radians."""
        return self.steering_step(theta_deg) * np.arange(self.elements)

    def element_quantized_phases(
        self, theta_deg: float, order: int
    ) -> NDArray[np.float64]:
        """Round every element of the ideal ramp to the ``C_order`` grid."""
        return quantize_phases(self.exact_phases(theta_deg), order)

    def step_quantized_phases(
        self, theta_deg: float, order: int
    ) -> NDArray[np.float64]:
        """Quantise the ramp step once and apply it uniformly."""
        quantum = TWO_PI / _phase_order(order)
        step = np.round(self.steering_step(theta_deg) / quantum) * quantum
        return step * np.arange(self.elements)

    def wrapped_phases(self, phases: ArrayLike) -> NDArray[np.float64]:
        """Return applied phases in ``[0, 2*pi)``, as a shifter would store."""
        return np.mod(_finite_vector(phases, "phases"), TWO_PI)

    def array_factor(
        self,
        theta_deg: ArrayLike,
        phases: ArrayLike,
        *,
        frequency_ratio: float = 1.0,
    ) -> NDArray[np.float64]:
        """Normalised far-field array factor at one or more angles.

        The result is ``|sum_n exp(j*(psi_n - k*x_n*sin(theta)))| / N``. It is
        one at the peak of an ideal, unquantised ramp and never exceeds one.
        ``frequency_ratio`` scales the wave number, which is what reveals
        wideband squint for a fixed phase ramp.
        """
        if not np.isfinite(frequency_ratio) or frequency_ratio <= 0.0:
            raise ValueError("frequency_ratio must be finite and positive")
        angles = np.atleast_1d(np.asarray(theta_deg, dtype=np.float64))
        if not np.all(np.isfinite(angles)):
            raise ValueError("theta_deg must contain only finite values")
        psi = _finite_vector(phases, "phases")
        if psi.size != self.elements:
            raise ValueError("phases must contain one value per element")
        wave_number = TWO_PI * self.spacing_over_lambda * frequency_ratio
        index = np.arange(self.elements)
        phase = psi[None, :] - wave_number * index[None, :] * np.sin(
            np.radians(angles)
        )[:, None]
        return np.abs(np.exp(1j * phase).mean(axis=1))

    def pattern_db(
        self,
        theta_deg: ArrayLike,
        phases: ArrayLike,
        *,
        frequency_ratio: float = 1.0,
    ) -> NDArray[np.float64]:
        """Array factor in decibels, floored at -120 dB."""
        magnitude = self.array_factor(
            theta_deg, phases, frequency_ratio=frequency_ratio
        )
        return 20.0 * np.log10(np.maximum(magnitude, 1e-6))

    def peak_angle(
        self, phases: ArrayLike, resolution_deg: float = 0.002
    ) -> float:
        """Realised main-beam angle, found by scanning the array factor."""
        grid = np.arange(-89.0, 89.0 + resolution_deg, resolution_deg)
        magnitude = self.array_factor(grid, phases)
        value = float(grid[int(np.argmax(magnitude))])
        return 0.0 if abs(value) <= resolution_deg / 2.0 else value

    def scan_loss_db(self, theta_deg: float, phases: ArrayLike) -> float:
        """Gain at the commanded angle relative to an ideal ramp (dB)."""
        value = float(self.array_factor(theta_deg, phases)[0])
        return float(20.0 * np.log10(max(value, 1e-6)))

    def half_power_width_deg(
        self,
        theta_deg: float,
        phases: ArrayLike,
        *,
        span_deg: float = 15.0,
        resolution_deg: float = 0.001,
    ) -> float:
        """Numerical -3 dB beamwidth around the realised peak."""
        peak = self.peak_angle(phases)
        lower = max(peak - span_deg, -89.0)
        upper = min(peak + span_deg, 89.0)
        grid = np.arange(lower, upper + resolution_deg, resolution_deg)
        magnitude = self.array_factor(grid, phases)
        below = np.flatnonzero(magnitude < 1.0 / np.sqrt(2.0))
        peak_index = int(np.argmin(np.abs(grid - peak)))
        left = below[below < peak_index]
        right = below[below > peak_index]
        if left.size == 0 or right.size == 0:
            return float("nan")
        return float(grid[right[0]] - grid[left[-1]])

    def phase_residual(
        self, phases: ArrayLike, theta_deg: float
    ) -> NDArray[np.float64]:
        """Wrapped difference between applied and ideal element phases."""
        applied = _finite_vector(phases, "phases")
        if applied.size != self.elements:
            raise ValueError("phases must contain one value per element")
        ideal = self.exact_phases(theta_deg)
        return np.angle(np.exp(1j * (applied - ideal)))

    def aberration_rms(self, phases: ArrayLike, theta_deg: float) -> float:
        """RMS phase aberration after removing piston and tilt, in radians.

        A linear phase is a pointing change, not a distortion, so the best-fit
        ramp is removed before the root-mean-square is taken. The residual is
        unwrapped first so a perfectly linear but modulo-``2*pi`` staircase
        does not masquerade as aberration.
        """
        residual = np.unwrap(self.phase_residual(phases, theta_deg))
        index = np.arange(self.elements)
        trend = np.polyval(np.polyfit(index, residual, 1), index)
        return float(np.sqrt(np.mean((residual - trend) ** 2)))

    def strehl_estimate(self, phases: ArrayLike, theta_deg: float) -> float:
        """Marechal estimate ``exp(-sigma^2)`` of the aberration loss."""
        rms = self.aberration_rms(phases, theta_deg)
        return float(np.exp(-(rms**2)))

    def field(
        self,
        x: ArrayLike,
        z: ArrayLike,
        phases: ArrayLike,
        *,
        amplitude_decay: bool = True,
    ) -> NDArray[np.complex128]:
        """Synthesise the scalar field on an ``(z, x)`` grid.

        Each element contributes a cylindrical wave
        ``exp(j*(k*r + psi_n)) / sqrt(1 + r/lambda)``. This is the physical
        near field, not the array-factor approximation, so the wavefront
        curvature out to the Fresnel distance is visible.
        """
        x_values = _finite_vector(x, "x")
        z_values = _finite_vector(z, "z")
        psi = _finite_vector(phases, "phases")
        if psi.size != self.elements:
            raise ValueError("phases must contain one value per element")
        X, Z = np.meshgrid(x_values, z_values)
        wave_number = TWO_PI / self.wavelength
        total = np.zeros_like(X, dtype=np.complex128)
        for position, phase in zip(self.positions(), psi, strict=True):
            radius = np.hypot(X - position, Z)
            wave = np.exp(1j * (wave_number * radius + phase))
            if amplitude_decay:
                wave = wave / np.sqrt(1.0 + radius / self.wavelength)
            total += wave
        return total

    def far_field_cut(
        self,
        angles_deg: ArrayLike,
        phases: ArrayLike,
        *,
        radius: float = 10.0,
    ) -> NDArray[np.complex128]:
        """Sum real cylindrical waves on a large circle (independent check)."""
        angles = _finite_vector(angles_deg, "angles_deg")
        psi = _finite_vector(phases, "phases")
        if psi.size != self.elements:
            raise ValueError("phases must contain one value per element")
        if not np.isfinite(radius) or radius <= 0.0:
            raise ValueError("radius must be finite and positive")
        X = radius * np.sin(np.radians(angles))
        Z = radius * np.cos(np.radians(angles))
        wave_number = TWO_PI / self.wavelength
        total = np.zeros(angles.size, dtype=np.complex128)
        for position, phase in zip(self.positions(), psi, strict=True):
            distance = np.hypot(X - position, Z)
            total += np.exp(1j * (wave_number * distance + phase)) / np.sqrt(
                1.0 + distance / self.wavelength
            )
        return total

    def far_field_peak(self, phases: ArrayLike, resolution_deg: float = 0.01) -> float:
        """Peak angle of the full field sum, with no array-factor formula."""
        grid = np.arange(-89.0, 89.0 + resolution_deg, resolution_deg)
        cut = self.far_field_cut(grid, phases)
        return float(grid[int(np.argmax(np.abs(cut)))])
