# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Bit-error-rate simulations over Eb/N0."""

from collections.abc import Mapping, Sequence

import numpy as np
from numpy.typing import NDArray

from modcode.channel.awgn import AWGNChannel
from modcode.fec.base import FEC
from modcode.modulation.base import Modulation
from modcode.pipeline.chain import Link, run_link


def ber_sweep(
    modulations: Mapping[str, Modulation] | Sequence[Link],
    ebn0_db: NDArray[np.floating] | Sequence[float],
    n_bits: int = 100_000,
    fec: FEC | None = None,
    seed: int | None = None,
) -> dict[str, NDArray[np.float64]]:
    """Estimate BER for each named modulation at each Eb/N0 point.

    When a sequence of :class:`Link` objects is provided, each link's own FEC
    and channel are used. A mapping of modulation names uses AWGN by default
    and applies the optional shared FEC.
    """
    points = np.asarray(ebn0_db, dtype=np.float64)
    if points.ndim != 1 or points.size == 0:
        raise ValueError("ebn0_db must be a non-empty one-dimensional sequence")
    if not np.all(np.isfinite(points)):
        raise ValueError("ebn0_db values must be finite")
    if n_bits <= 0:
        raise ValueError("n_bits must be positive")

    if isinstance(modulations, Mapping):
        if not modulations:
            raise ValueError("at least one modulation is required")
        links = {
            name: Link(modulation, AWGNChannel(), fec)
            for name, modulation in modulations.items()
        }
    else:
        if not modulations:
            raise ValueError("at least one link is required")
        links = {
            f"{type(link.modulation).__name__}"
            + (f"-{link.modulation.order}" if hasattr(link.modulation, "order") else "")
            + (f"-{index}" if index else ""): link
            for index, link in enumerate(modulations)
        }

    generator = np.random.default_rng(seed)
    results: dict[str, NDArray[np.float64]] = {}
    for name, link in links.items():
        ber = np.empty(points.size, dtype=np.float64)
        for index, point in enumerate(points):
            ber[index] = run_link(
                link, n_bits=n_bits, ebn0_db=float(point), rng=generator
            ).ber
        results[name] = ber
    return results
