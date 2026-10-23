# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Gain accounting for sphere packing and multidimensional modulation.

Three different gains separate a lattice constellation from square QAM, and
they are routinely added together by mistake:

* the *nominal coding gain*, ``10 log10(gamma)`` with the Hermite parameter
  ``gamma = d**2 / det**(2/n)`` of the underlying lattice;
* the *shaping gain*, ``(1/12) / G``, where ``G`` is the normalized second
  moment of the constellation region (a cube has ``G = 1/12``);
* the *gap to capacity*, which is the quantity the first two are measured
  against and which neither of them closes.

The helpers here are deliberately pure functions of published invariants, so
the note in ``docs/sphere-packing-modulation.md`` can show its arithmetic and
``tests/test_sphere_packing.py`` can pin the literature values.
"""

from __future__ import annotations

import math
from typing import NamedTuple

import numpy as np

from modcode.metrics.theory import theoretical_ber
from modcode.modulation.base import Modulation


class LatticeFacts(NamedTuple):
    """Packing invariants of a lattice, in its own dimension ``n``."""

    name: str
    n: int
    min_squared_norm: float
    determinant: float


#: Densest known lattices in their dimensions; ``determinant`` is the volume
#: of the fundamental cell and ``min_squared_norm`` the squared length of a
#: shortest non-zero vector.
LATTICES: tuple[LatticeFacts, ...] = (
    LatticeFacts("Z^2", 2, 1.0, 1.0),
    LatticeFacts("D4", 4, 2.0, 2.0),
    LatticeFacts("E8", 8, 2.0, 1.0),
    LatticeFacts("Leech", 24, 4.0, 1.0),
)

#: Normalized second moments ``G`` of Voronoi cells, Conway and Sloane,
#: *Sphere Packings, Lattices and Groups*, 3rd ed., Table 2.3.  The cube has
#: ``G = 1/12``; the lower bound over all regions is the sphere, ``1/(2 pi e)``.
NORMALIZED_SECOND_MOMENTS: dict[str, float] = {
    "cube": 1.0 / 12.0,
    "D4": 0.0766032,
    "E8": 0.0716821,
    "Leech": 0.0657710,
}

CUBE_NORMALIZED_SECOND_MOMENT = 1.0 / 12.0
SPHERE_NORMALIZED_SECOND_MOMENT = 1.0 / (2.0 * math.pi * math.e)


def hermite_parameter(lattice: LatticeFacts) -> float:
    """Return ``d**2 / det**(2/n)``, the dimensionless packing figure."""
    exponent = 2.0 / lattice.n
    return lattice.min_squared_norm / lattice.determinant**exponent


def nominal_coding_gain_db(lattice: LatticeFacts) -> float:
    """Return the asymptotic gain of the lattice over ``Z**n`` in dB."""
    return 10.0 * math.log10(hermite_parameter(lattice))


def center_density(lattice: LatticeFacts) -> float:
    """Return the fraction of space covered by the packing spheres."""
    radius = 0.5 * math.sqrt(lattice.min_squared_norm)
    return radius**lattice.n / lattice.determinant


def shaping_gain_db(normalized_second_moment: float) -> float:
    """Return the gain of a region over the cube with the same volume, in dB."""
    if not (
        SPHERE_NORMALIZED_SECOND_MOMENT
        <= normalized_second_moment
        <= CUBE_NORMALIZED_SECOND_MOMENT
    ):
        raise ValueError("G must lie between the sphere's 1/(2 pi e) and 1/12")
    return 10.0 * math.log10(CUBE_NORMALIZED_SECOND_MOMENT / normalized_second_moment)


#: The shaping gain ceiling, ``pi e / 6`` in linear terms: 1.53 dB.
SPHERE_SHAPING_GAIN_DB = shaping_gain_db(SPHERE_NORMALIZED_SECOND_MOMENT)


def shannon_ebn0_db(eta: float) -> float:
    """Return the Shannon ``Eb/N0`` limit in dB for ``eta`` bits per real axis.

    The real-input AWGN capacity is ``0.5 log2(1 + 2 eta Eb/N0)``, so the
    limit for spectral efficiency ``eta`` sits at ``(2**(2 eta) - 1) / (2 eta)``.
    """
    if eta <= 0.0:
        raise ValueError("eta must be positive")
    return 10.0 * math.log10((2.0 ** (2.0 * eta) - 1.0) / (2.0 * eta))


def ebn0_at_ber(
    modulation: Modulation,
    target_ber: float,
    bracket: tuple[float, float] = (-5.0, 40.0),
    iterations: int = 200,
) -> float:
    """Bisect the closed-form uncoded BER reference for ``target_ber``."""
    if not 0.0 < target_ber < 1.0:
        raise ValueError("target_ber must lie strictly between 0 and 1")
    low, high = bracket
    if not float(theoretical_ber(modulation, np.array([low]))[0]) > target_ber:
        raise ValueError("low end of the bracket is already below the target")
    for _ in range(iterations):
        middle = 0.5 * (low + high)
        reached = float(theoretical_ber(modulation, np.array([middle]))[0])
        if reached > target_ber:
            low = middle
        else:
            high = middle
    return 0.5 * (low + high)


def gap_to_capacity_db(modulation: Modulation, eta: float, target_ber: float) -> float:
    """Return how far uncoded operation sits above the Shannon limit, in dB."""
    return ebn0_at_ber(modulation, target_ber) - shannon_ebn0_db(eta)


__all__ = [
    "CUBE_NORMALIZED_SECOND_MOMENT",
    "LATTICES",
    "NORMALIZED_SECOND_MOMENTS",
    "SPHERE_NORMALIZED_SECOND_MOMENT",
    "SPHERE_SHAPING_GAIN_DB",
    "LatticeFacts",
    "center_density",
    "ebn0_at_ber",
    "gap_to_capacity_db",
    "hermite_parameter",
    "nominal_coding_gain_db",
    "shaping_gain_db",
    "shannon_ebn0_db",
]
