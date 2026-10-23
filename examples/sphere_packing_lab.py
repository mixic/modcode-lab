#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Sphere-packing gain accounting: the three tables of the companion note.

Run from the repository root:

    python examples/sphere_packing_lab.py

The tables answer three separate questions, and the note in
``docs/sphere-packing-modulation.md`` quotes them:

1. how far uncoded square QAM sits above the Shannon limit at a fixed BER;
2. what the dense lattices would give if they were available at that size;
3. how much of the answer is shaping and how much is packing.

Only the first table is a measured property of a real decoder; the other two
are properties of the asymptotic objects, which is precisely the point.
"""

from __future__ import annotations

import numpy as np

from modcode.metrics.sphere_packing import (
    LATTICES,
    NORMALIZED_SECOND_MOMENTS,
    SPHERE_SHAPING_GAIN_DB,
    center_density,
    ebn0_at_ber,
    hermite_parameter,
    nominal_coding_gain_db,
    shannon_ebn0_db,
    shaping_gain_db,
)
from modcode.modulation.qam import QAM
from modcode.modulation.qpsk import QPSK

TARGET_BER = 1.0e-5


def gap_table() -> list[tuple[float, str, float, float, float]]:
    """Return eta, scheme, Shannon limit, uncoded Eb/N0, and the gap in dB."""
    schemes = (
        (1.0, "QPSK", QPSK()),
        (2.0, "16-QAM", QAM(16)),
        (3.0, "64-QAM", QAM(64)),
        (4.0, "256-QAM", QAM(256)),
    )
    rows = []
    for eta, name, modulation in schemes:
        limit = shannon_ebn0_db(eta)
        uncoded = ebn0_at_ber(modulation, TARGET_BER)
        rows.append((eta, name, limit, uncoded, uncoded - limit))
    return rows


def lattice_table() -> list[tuple[str, int, float, float, float, float]]:
    """Return name, n, d^2, det, Hermite parameter, and nominal gain in dB."""
    return [
        (
            lattice.name,
            lattice.n,
            lattice.min_squared_norm,
            lattice.determinant,
            hermite_parameter(lattice),
            nominal_coding_gain_db(lattice),
        )
        for lattice in LATTICES
    ]


def shaping_table() -> list[tuple[str, float, float]]:
    """Return region name, normalized second moment, and shaping gain in dB."""
    rows = [
        (name, moment, shaping_gain_db(moment))
        for name, moment in NORMALIZED_SECOND_MOMENTS.items()
    ]
    rows.append(("sphere", 1.0 / (2.0 * np.pi * np.e), SPHERE_SHAPING_GAIN_DB))
    return rows


def main() -> None:
    """Print the three tables."""
    print(f"Gap to capacity at BER = {TARGET_BER:.0e} (uncoded, closed form)")
    print(f"{'eta':>4} {'scheme':>9} {'Shannon':>8} {'uncoded':>8} {'gap':>7}")
    for eta, name, limit, uncoded, gap in gap_table():
        print(f"{eta:4.0f} {name:>9} {limit:8.3f} {uncoded:8.3f} {gap:7.3f}")

    print()
    print("Lattice invariants (asymptotic, in the lattice's own dimension)")
    print(f"{'lattice':>8} {'n':>3} {'d^2':>5} {'det':>5} {'gamma':>7} {'gain':>7}")
    for name, n, d2, det, gamma, gain in lattice_table():
        print(f"{name:>8} {n:3d} {d2:5.1f} {det:5.1f} {gamma:7.4f} {gain:7.3f}")

    print()
    print("Centre densities")
    for lattice in LATTICES:
        print(f"{lattice.name:>8} {center_density(lattice):10.5f}")

    print()
    print("Shaping gain of a region over the cube of equal volume")
    for name, moment, gain in shaping_table():
        print(f"{name:>8} G={moment:.7f} gain={gain:7.3f} dB")

    print()
    print(f"shaping ceiling pi*e/6 = {SPHERE_SHAPING_GAIN_DB:.4f} dB")


if __name__ == "__main__":
    main()
