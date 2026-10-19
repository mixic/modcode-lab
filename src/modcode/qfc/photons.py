# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Photon budgets and phase-insensitive frequency conversion.

Everything here uses shot-noise units (SNU): the vacuum variance of one
quadrature is one, a coherent state has unit covariance, and a thermal state
with mean occupation ``n`` has covariance ``(2*n + 1) * I``.

A phase-insensitive transducer of efficiency ``eta`` and output-referred added
noise ``n_add`` maps an input covariance ``V`` to

    V_out = eta * V + (1 - eta + n_add) * I

The ``(1 - eta)`` term is the vacuum that leaks in through the loss. If that
loss is calibrated, it is not excess noise, and the input-referred excess noise
is ``n_add / eta``. If the receiver cannot trust ``eta`` -- because it drifts,
or because someone could have inserted a lossy component -- the same loss shows
up as ``(1 - eta + n_add) / eta`` SNU, indistinguishable from an eavesdropper.
That distinction, not the efficiency alone, decides whether a converted
microwave state is still useful.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

PLANCK = 6.626_070_15e-34
"""Planck constant, in joule-seconds."""

BOLTZMANN = 1.380_649e-23
"""Boltzmann constant, in joules per kelvin."""

SPEED_OF_LIGHT = 299_792_458.0
"""Speed of light in vacuum, in metres per second."""


def photon_energy(frequency: float) -> float:
    """Energy of one photon at ``frequency``, in joules."""
    if not np.isfinite(frequency) or frequency <= 0.0:
        raise ValueError("frequency must be finite and positive")
    return PLANCK * frequency


def thermal_occupation(frequency: float, temperature: float) -> float:
    """Mean photon number of a single mode in thermal equilibrium.

    ``n = 1 / (exp(h*f / (k*T)) - 1)``. The zero-temperature limit is the
    vacuum, and the expression underflows gracefully to zero for ``h*f >> k*T``.
    """
    if not np.isfinite(frequency) or frequency <= 0.0:
        raise ValueError("frequency must be finite and positive")
    if not np.isfinite(temperature) or temperature < 0.0:
        raise ValueError("temperature must be finite and non-negative")
    if temperature == 0.0:
        return 0.0
    try:
        ratio = PLANCK * frequency / (BOLTZMANN * temperature)
        return 1.0 / math.expm1(ratio)
    except OverflowError:
        return 0.0


def squeezed_vacuum_covariance(squeezing_db: float) -> NDArray[np.float64]:
    """Covariance of a squeezed vacuum in SNU.

    ``squeezing_db`` is the reduction of the smaller quadrature variance below
    the vacuum, so 10 dB gives ``diag(0.1, 10)``.
    """
    if not np.isfinite(squeezing_db) or squeezing_db < 0.0:
        raise ValueError("squeezing_db must be finite and non-negative")
    variance = 10.0 ** (-squeezing_db / 10.0)
    return np.diag([variance, 1.0 / variance])


def covariance_squeezing_db(covariance: ArrayLike) -> float:
    """Squeezing of the quieter quadrature, in dB relative to vacuum.

    Negative values mean squeezing; a positive value means both quadratures
    are noisier than vacuum and no squeezing survives.
    """
    matrix = _validated_covariance(covariance)
    smaller = float(np.linalg.eigvalsh(matrix)[0])
    return float(10.0 * np.log10(max(smaller, 1e-12)))


def thermal_squeezed_covariance(
    squeezing_db: float, occupation: float
) -> NDArray[np.float64]:
    """Squeezed vacuum embedded in a thermal bath of ``occupation`` photons."""
    if not np.isfinite(occupation) or occupation < 0.0:
        raise ValueError("occupation must be finite and non-negative")
    return (2.0 * occupation + 1.0) * squeezed_vacuum_covariance(squeezing_db)


def _validated_covariance(covariance: ArrayLike) -> NDArray[np.float64]:
    matrix = np.asarray(covariance, dtype=np.float64)
    if matrix.shape != (2, 2):
        raise ValueError("covariance must be a 2x2 matrix")
    if not np.all(np.isfinite(matrix)):
        raise ValueError("covariance must contain only finite values")
    if not np.allclose(matrix, matrix.T, atol=1e-12):
        raise ValueError("covariance must be symmetric")
    if np.linalg.eigvalsh(matrix)[0] <= 0.0:
        raise ValueError("covariance must be positive definite")
    return matrix


@dataclass(frozen=True)
class Transducer:
    """A phase-insensitive microwave-to-optical converter.

    Parameters
    ----------
    efficiency
        Photon-number conversion efficiency in ``(0, 1]``.
    added_noise
        Output-referred added noise in SNU (photons per mode). Zero is an
        ideal loss-only converter.
    """

    efficiency: float
    added_noise: float = 0.0

    def __post_init__(self) -> None:
        if not np.isfinite(self.efficiency) or not 0.0 < self.efficiency <= 1.0:
            raise ValueError("efficiency must lie in (0, 1]")
        if not np.isfinite(self.added_noise) or self.added_noise < 0.0:
            raise ValueError("added_noise must be finite and non-negative")

    def excess_noise(self, *, calibrated: bool = True) -> float:
        """Input-referred excess noise in SNU.

        With ``calibrated=True`` the loss is trusted and only the added noise
        counts. With ``calibrated=False`` the loss itself is treated as
        unknown noise, which is the conservative (Eve-compatible) case.
        """
        if calibrated:
            return self.added_noise / self.efficiency
        return (1.0 - self.efficiency + self.added_noise) / self.efficiency

    def convert_covariance(self, covariance: ArrayLike) -> NDArray[np.float64]:
        """Apply the converter to a Gaussian state covariance."""
        matrix = _validated_covariance(covariance)
        return self.efficiency * matrix + (
            1.0 - self.efficiency + self.added_noise
        ) * np.eye(2)

    def surviving_squeezing_db(
        self, squeezing_db: float, *, bath_occupation: float = 0.0
    ) -> float:
        """Squeezing left after the state passes through the converter."""
        input_covariance = thermal_squeezed_covariance(squeezing_db, bath_occupation)
        return covariance_squeezing_db(self.convert_covariance(input_covariance))

    def coherent_output_occupation(self, input_occupation: float) -> float:
        """Mean photons of a coherent state after conversion (added noise only).

        Loss attenuates the mean but does not add photons at zero temperature;
        the transducer's added noise does.
        """
        if not np.isfinite(input_occupation) or input_occupation < 0.0:
            raise ValueError("input_occupation must be finite and non-negative")
        return self.efficiency * input_occupation + self.added_noise


def wigner_function(
    q: ArrayLike,
    p: ArrayLike,
    covariance: ArrayLike,
    *,
    mean: ArrayLike = (0.0, 0.0),
) -> NDArray[np.float64]:
    """Wigner function of a zero- or displaced-Gaussian state on a grid.

    ``covariance`` is in SNU and ``mean`` is the displacement
    ``(q, p)``. The normalisation is ``integral W dq dp = 1``.
    """
    matrix = _validated_covariance(covariance)
    center = np.asarray(mean, dtype=np.float64)
    if center.shape != (2,) or not np.all(np.isfinite(center)):
        raise ValueError("mean must contain two finite values")
    q_values = np.asarray(q, dtype=np.float64)
    p_values = np.asarray(p, dtype=np.float64)
    if q_values.ndim != 1 or p_values.ndim != 1:
        raise ValueError("q and p must be one-dimensional arrays")
    Q, P = np.meshgrid(q_values, p_values)
    delta_q = Q - center[0]
    delta_p = P - center[1]
    inverse = np.linalg.inv(matrix)
    exponent = 0.5 * (
        inverse[0, 0] * delta_q**2
        + 2.0 * inverse[0, 1] * delta_q * delta_p
        + inverse[1, 1] * delta_p**2
    )
    determinant = float(np.linalg.det(matrix))
    return np.exp(-exponent) / (2.0 * np.pi * math.sqrt(determinant))
