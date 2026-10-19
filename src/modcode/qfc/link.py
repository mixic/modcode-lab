# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Transport loss for microwave, free-space optical, and fibre carriers.

The point of these numbers is not precision link engineering. It is the order
of magnitude that makes quantum frequency conversion worth its noise: a
10 GHz state cannot be moved 100 km over any room-temperature conductor, while
1550 nm light crosses the same distance on fibre for about 22 dB.

Free-space optical figures assume diffraction-limited apertures and a Gaussian
beam. They deliberately exclude turbulence, cloud, and weather, all of which
can add tens of dB and availability limits that no modulation experiment can
wish away.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

SPEED_OF_LIGHT = 299_792_458.0


@dataclass(frozen=True)
class TransportPath:
    """One carrier/medium option and its loss over a fixed distance."""

    label: str
    loss_db: float


def _positive(value: float, name: str) -> float:
    if not np.isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be finite and positive")
    return float(value)


def free_space_path_loss(distance_m: float, wavelength_m: float) -> float:
    """Friis free-space path loss in dB."""
    distance = _positive(distance_m, "distance_m")
    wavelength = _positive(wavelength_m, "wavelength_m")
    return float(20.0 * np.log10(4.0 * np.pi * distance / wavelength))


def aperture_gain_db(diameter_m: float, wavelength_m: float) -> float:
    """On-axis gain of a uniformly illuminated circular aperture, in dB."""
    diameter = _positive(diameter_m, "diameter_m")
    wavelength = _positive(wavelength_m, "wavelength_m")
    return float(20.0 * np.log10(np.pi * diameter / wavelength))


def pointing_loss_db(
    distance_m: float,
    wavelength_m: float,
    aperture_m: float,
    jitter_rad: float,
) -> float:
    """Gaussian-beam coupling loss caused by angular pointing jitter.

    A pointing error of ``jitter_rad`` displaces the beam by
    ``jitter * distance`` at the receiver. The Gaussian mode radius there is
    ``w0 * sqrt(1 + (distance / z_R)^2)``, and the power coupling falls as
    ``exp(-(offset / w)^2)``. This is a single-axis, small-angle estimate.
    """
    distance = _positive(distance_m, "distance_m")
    wavelength = _positive(wavelength_m, "wavelength_m")
    aperture = _positive(aperture_m, "aperture_m")
    if not np.isfinite(jitter_rad) or jitter_rad < 0.0:
        raise ValueError("jitter_rad must be finite and non-negative")
    if jitter_rad == 0.0:
        return 0.0
    waist = aperture / 2.0
    rayleigh = np.pi * waist**2 / wavelength
    beam_radius = waist * np.sqrt(1.0 + (distance / rayleigh) ** 2)
    offset = jitter_rad * distance
    coupling = float(np.exp(-((offset / beam_radius) ** 2)))
    return float(-10.0 * np.log10(max(coupling, 1e-300)))


def free_space_link_loss(
    distance_m: float,
    wavelength_m: float,
    tx_aperture_m: float,
    rx_aperture_m: float,
    *,
    pointing_jitter_rad: float = 0.0,
) -> float:
    """Two-aperture free-space link loss in dB.

    The Friis term already contains the diffraction spreading; the optional
    jitter term adds a Gaussian pointing penalty on top.
    """
    loss = (
        free_space_path_loss(distance_m, wavelength_m)
        - aperture_gain_db(tx_aperture_m, wavelength_m)
        - aperture_gain_db(rx_aperture_m, wavelength_m)
        + pointing_loss_db(
            distance_m, wavelength_m, rx_aperture_m, pointing_jitter_rad
        )
    )
    return float(loss)


def fiber_loss(
    distance_m: float,
    *,
    attenuation_db_per_km: float = 0.2,
    connector_loss_db: float = 2.0,
) -> float:
    """Standard single-mode fibre loss at 1550 nm."""
    distance = _positive(distance_m, "distance_m")
    if attenuation_db_per_km < 0.0 or connector_loss_db < 0.0:
        raise ValueError("attenuation and connector loss must be non-negative")
    return float(attenuation_db_per_km * distance / 1e3 + connector_loss_db)


def coax_loss(
    distance_m: float, *, attenuation_db_per_m: float = 0.5
) -> float:
    """Room-temperature coaxial cable loss."""
    distance = _positive(distance_m, "distance_m")
    if attenuation_db_per_m < 0.0:
        raise ValueError("attenuation_db_per_m must be non-negative")
    return float(attenuation_db_per_m * distance)


def transport_budget(
    distance_m: float = 100e3,
    *,
    frequency: float = 10.0e9,
    optical_wavelength_m: float = 1550e-9,
    aperture_m: float = 0.3,
    pointing_jitter_rad: float = 1e-6,
) -> list[TransportPath]:
    """Compare plausible 100 km routes for the same information."""
    distance = _positive(distance_m, "distance_m")
    rf_wavelength = SPEED_OF_LIGHT / _positive(frequency, "frequency")
    rows = [
        TransportPath(
            "RF coax, 10 GHz (0.5 dB/m)", coax_loss(distance, attenuation_db_per_m=0.5)
        ),
        TransportPath(
            "RF waveguide (0.05 dB/m)",
            coax_loss(distance, attenuation_db_per_m=0.05),
        ),
        TransportPath(
            "Free-space RF, 0.3 m dishes",
            free_space_link_loss(
                distance, rf_wavelength, aperture_m, aperture_m
            ),
        ),
        TransportPath(
            "Free-space optical, no jitter",
            free_space_link_loss(
                distance, optical_wavelength_m, aperture_m, aperture_m
            ),
        ),
        TransportPath(
            "Free-space optical, 1 urad jitter",
            free_space_link_loss(
                distance,
                optical_wavelength_m,
                aperture_m,
                aperture_m,
                pointing_jitter_rad=pointing_jitter_rad,
            ),
        ),
        TransportPath(
            "Fiber 1550 nm (0.2 dB/km)",
            fiber_loss(distance),
        ),
    ]
    return rows
