# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""BER simulation for multi-dimensional lattice constellations.

The comparison is at equal **energy per real dimension**, not merely equal
bits per symbol.  For a constellation with ``b`` bits per symbol and ``K``
real dimensions, let ``eta = b/K``.  With unit energy per real dimension,
``Eb = 1/eta`` and the noise variance per real dimension at a given Eb/N0 is

    sigma^2 = 1 / (2 * eta * Eb/N0).

This makes a D4-256 symbol and a 16-QAM symbol directly comparable: both
carry 2 bits per real dimension and see the same per-dimension noise.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from modcode.modulation.base import Modulation
from modcode.modulation.lattice import LatticeModulation


def qam_per_dimension_points(
    modulation: Modulation,
) -> tuple[NDArray[np.float64], NDArray[np.uint8]]:
    """Return 2D QAM/QPSK points at unit energy per real dimension."""
    labels = np.asarray(
        list(modulation.constellation.keys()), dtype=np.uint8
    )
    complex_points = np.asarray(
        list(modulation.constellation.values()), dtype=np.complex128
    )
    points = np.column_stack([complex_points.real, complex_points.imag])
    points = points * np.sqrt(2.0)
    return points, labels


def simulate_point_ber(
    points: NDArray[np.float64],
    labels: NDArray[np.uint8],
    bits_per_symbol: int,
    ebn0_db: float,
    n_bits: int,
    rng: np.random.Generator,
    *,
    batch: int = 1_000,
) -> float:
    """Simulate hard-decision BER for real points at unit per-dim energy."""
    if points.ndim != 2 or labels.ndim != 2:
        raise ValueError("points and labels must be two-dimensional")
    if points.shape[0] != labels.shape[0]:
        raise ValueError("points and labels must have the same length")
    if labels.shape[1] != bits_per_symbol:
        raise ValueError("labels do not match bits_per_symbol")
    if n_bits <= 0 or not np.isfinite(ebn0_db):
        raise ValueError("n_bits must be positive and eb/n0 finite")
    if batch <= 0:
        raise ValueError("batch must be positive")

    dimensions = points.shape[1]
    eta = bits_per_symbol / dimensions
    variance = 1.0 / (2.0 * eta * 10.0 ** (ebn0_db / 10.0))
    symbols = n_bits // bits_per_symbol
    label_map = {
        tuple(pattern): index for index, pattern in enumerate(labels)
    }
    errors = 0
    used = 0
    for start in range(0, symbols, batch):
        stop = min(symbols, start + batch)
        bits = rng.integers(
            0, 2, size=(stop - start, bits_per_symbol), dtype=np.uint8
        )
        indices = np.asarray(
            [label_map[tuple(row)] for row in bits], dtype=np.intp
        )
        transmitted = points[indices]
        received = transmitted + np.sqrt(variance) * rng.standard_normal(
            transmitted.shape
        )
        distances = np.sum(
            (received[:, None, :] - points[None, :, :]) ** 2, axis=2
        )
        decoded = labels[np.argmin(distances, axis=1)]
        errors += int(np.count_nonzero(decoded != bits))
        used += bits.size
    return errors / used if used else 0.0


def simulate_lattice_ber(
    modulation: LatticeModulation,
    ebn0_db: float,
    n_bits: int,
    rng: np.random.Generator,
) -> float:
    """Simulate hard-decision BER for a lattice constellation."""
    return simulate_point_ber(
        modulation.per_dimension_points(1.0),
        modulation.bit_patterns,
        modulation.bits_per_symbol,
        ebn0_db,
        n_bits,
        rng,
    )


__all__ = [
    "qam_per_dimension_points",
    "simulate_lattice_ber",
    "simulate_point_ber",
]
