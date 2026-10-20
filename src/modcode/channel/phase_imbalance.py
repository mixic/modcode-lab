# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""I/Q phase-imbalance channel.

A quadrature modulator with a phase error ``phi`` in the Q branch produces,
for the complex baseband symbol ``x = I + jQ``,

    y = I + j Q exp(j phi)
      = (I - Q sin(phi)) + j (Q cos(phi)),

before additive noise.  The real 2-by-2 matrix acting on ``(I, Q)`` is

    R(phi) = [[1, -sin(phi)],
              [0,  cos(phi)]].

The determinant is ``cos(phi)``.  At ``phi = 0`` the channel is the identity;
at ``phi = pi`` it is the Q-axis reflection ``y = conj(x)``, an involution;
at ``phi = pi/2`` the matrix is rank one, so the Q component is lost.

The common shorthand ``I + exp(j phi) Q`` omits the ``j`` that makes ``Q`` the
imaginary axis.  It gives ``I - Q`` at ``phi = pi``, a rank-one projection,
not a reflection and not an involution.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from modcode.channel.awgn import AWGNChannel
from modcode.channel.base import Channel


class PhaseImbalanceChannel(Channel):
    """Quadrature phase mismatch, followed by complex AWGN.

    The noise variance is normalized to the input symbol energy, exactly as in
    :class:`~modcode.channel.awgn.AWGNChannel`.  The impairment is not
    renormalized afterwards: at ``phi=pi`` the map preserves energy, while at
    other phases it can change the mean symbol energy.  That change is part of
    the modelled impairment, not a separate gain control.

    Parameters
    ----------
    phi
        Phase mismatch in radians.  ``0`` is ideal; ``pi`` reflects the Q
        axis; ``pi/2`` (mod ``pi``) is singular.
    """

    def __init__(self, phi: float = 0.0) -> None:
        value = float(phi)
        if not np.isfinite(value):
            raise ValueError("phi must be finite")
        self.phi = value

    @classmethod
    def from_degrees(cls, degrees: float) -> PhaseImbalanceChannel:
        """Build a channel from a phase mismatch in degrees."""
        value = float(degrees)
        if not np.isfinite(value):
            raise ValueError("degrees must be finite")
        return cls(np.deg2rad(value))

    def matrix(self) -> NDArray[np.float64]:
        """Return the real matrix acting on the column vector ``(I, Q)``."""
        cosine = np.cos(self.phi)
        sine = np.sin(self.phi)
        return np.asarray([[1.0, -sine], [0.0, cosine]], dtype=np.float64)

    def impair(self, symbols: NDArray[np.complexfloating]) -> NDArray[np.complex128]:
        """Apply the phase mismatch without adding noise."""
        values = np.asarray(symbols, dtype=np.complex128)
        if values.ndim != 1:
            raise ValueError("symbols must be a one-dimensional array")
        if not np.all(np.isfinite(values)):
            raise ValueError("symbols must contain only finite values")
        return values.real + 1j * values.imag * np.exp(1j * self.phi)

    def mean_energy_ratio(
        self, symbols: NDArray[np.complexfloating]
    ) -> float:
        """Return the impaired-to-input mean symbol-energy ratio."""
        values = np.asarray(symbols, dtype=np.complex128)
        if values.ndim != 1:
            raise ValueError("symbols must be a one-dimensional array")
        if not np.all(np.isfinite(values)):
            raise ValueError("symbols must contain only finite values")
        input_energy = float(np.mean(np.abs(values) ** 2))
        if input_energy == 0.0:
            raise ValueError("symbols must not be all zero")
        impaired = values.real + 1j * values.imag * np.exp(1j * self.phi)
        return float(np.mean(np.abs(impaired) ** 2) / input_energy)

    def transmit(
        self,
        symbols: NDArray[np.complexfloating],
        ebn0_db: float,
        bits_per_symbol: int,
        code_rate: float = 1.0,
        rng: np.random.Generator | None = None,
    ) -> NDArray[np.complex128]:
        """Apply the mismatch and then add AWGN at the requested Eb/N0."""
        impaired = self.impair(symbols)
        variance = AWGNChannel.noise_variance(
            ebn0_db, bits_per_symbol, code_rate
        )
        generator = rng if rng is not None else np.random.default_rng()
        noise = np.sqrt(variance) * (
            generator.standard_normal(impaired.size)
            + 1j * generator.standard_normal(impaired.size)
        )
        return impaired + noise


__all__ = ["PhaseImbalanceChannel"]
