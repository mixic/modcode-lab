# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Bit-error-rate simulations over Eb/N0."""

from __future__ import annotations

import hashlib
import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from statistics import NormalDist

import numpy as np
from numpy.typing import NDArray

from modcode.channel.awgn import AWGNChannel
from modcode.fec.base import FEC
from modcode.modulation.base import Modulation
from modcode.pipeline.chain import Link, run_link


@dataclass(frozen=True)
class BERCurve:
    """A BER estimate with counts and a binomial confidence interval."""

    ber: NDArray[np.float64]
    errors: NDArray[np.int64]
    bits: NDArray[np.int64]
    ci_low: NDArray[np.float64]
    ci_high: NDArray[np.float64]


def _wilson_interval(
    errors: int, bits: int, z: float
) -> tuple[float, float]:
    """Return the Wilson score interval for a binomial proportion."""
    if bits <= 0:
        return 0.0, 1.0
    p = errors / bits
    denominator = 1.0 + z * z / bits
    centre = (p + z * z / (2.0 * bits)) / denominator
    half = (
        z
        * math.sqrt(p * (1.0 - p) / bits + z * z / (4.0 * bits * bits))
        / denominator
    )
    return max(0.0, centre - half), min(1.0, centre + half)


def _link_seed(seed: int | None, name: str) -> np.random.SeedSequence:
    """Make a link's stream depend on its name, not on dict insertion order."""
    if seed is None:
        return np.random.SeedSequence()
    digest = hashlib.blake2b(name.encode("utf-8"), digest_size=8).digest()
    return np.random.SeedSequence(
        [int(seed), int.from_bytes(digest, "little")]
    )


def _links(
    modulations: Mapping[str, Modulation] | Sequence[Link],
    fec: FEC | None,
) -> dict[str, Link]:
    if isinstance(modulations, Mapping):
        if not modulations:
            raise ValueError("at least one modulation is required")
        return {
            name: Link(modulation, AWGNChannel(), fec)
            for name, modulation in modulations.items()
        }
    if not modulations:
        raise ValueError("at least one link is required")
    return {
        f"{type(link.modulation).__name__}-{index}": link
        for index, link in enumerate(modulations)
    }


def _simulate(
    link: Link,
    ebn0_db: float,
    *,
    min_errors: int,
    max_bits: int,
    batch_bits: int,
    generator: np.random.Generator,
) -> tuple[int, int]:
    """Run one link at one Eb/N0 until an error or bit target is reached."""
    errors = 0
    bits = 0
    while bits < max_bits and (bits == 0 or errors < min_errors):
        request = min(batch_bits, max_bits - bits)
        result = run_link(
            link,
            n_bits=request,
            ebn0_db=ebn0_db,
            rng=generator,
        )
        errors += result.bit_errors
        bits += request
    return errors, bits


def ber_sweep_detailed(
    modulations: Mapping[str, Modulation] | Sequence[Link],
    ebn0_db: NDArray[np.floating] | Sequence[float],
    fec: FEC | None = None,
    seed: int | None = None,
    *,
    min_errors: int = 100,
    max_bits: int = 1_000_000,
    batch_bits: int = 50_000,
    confidence: float = 0.95,
) -> dict[str, BERCurve]:
    """Estimate BER with a stopping rule and binomial confidence intervals.

    A point stops when it has accumulated ``min_errors`` errors or
    ``max_bits`` bits, whichever comes first.  Bits are simulated in chunks
    of ``batch_bits``.  The confidence interval is the Wilson score interval,
    which remains informative at zero errors: a zero-error point still has a
    nonzero upper bound and can be plotted on a log axis using ``ci_high``.
    """
    points = np.asarray(ebn0_db, dtype=np.float64)
    if points.ndim != 1 or points.size == 0:
        raise ValueError("ebn0_db must be a non-empty one-dimensional sequence")
    if not np.all(np.isfinite(points)):
        raise ValueError("ebn0_db values must be finite")
    if min_errors < 0:
        raise ValueError("min_errors must be non-negative")
    if max_bits <= 0 or batch_bits <= 0:
        raise ValueError("max_bits and batch_bits must be positive")
    if not 0.0 < confidence < 1.0:
        raise ValueError("confidence must be in the interval (0, 1)")

    z = NormalDist().inv_cdf(0.5 + confidence / 2.0)
    links = _links(modulations, fec)
    results: dict[str, BERCurve] = {}
    for name, link in links.items():
        ber = np.empty(points.size, dtype=np.float64)
        errors_out = np.empty(points.size, dtype=np.int64)
        bits_out = np.empty(points.size, dtype=np.int64)
        low_out = np.empty(points.size, dtype=np.float64)
        high_out = np.empty(points.size, dtype=np.float64)
        point_seeds = _link_seed(seed, name).spawn(points.size)
        for index, point in enumerate(points):
            errors, bits = _simulate(
                link,
                float(point),
                min_errors=min_errors,
                max_bits=max_bits,
                batch_bits=batch_bits,
                generator=np.random.default_rng(point_seeds[index]),
            )
            low, high = _wilson_interval(errors, bits, z)
            errors_out[index] = errors
            bits_out[index] = bits
            ber[index] = errors / bits if bits else 0.0
            low_out[index] = low
            high_out[index] = high
        results[name] = BERCurve(
            ber=ber,
            errors=errors_out,
            bits=bits_out,
            ci_low=low_out,
            ci_high=high_out,
        )
    return results


def ber_sweep(
    modulations: Mapping[str, Modulation] | Sequence[Link],
    ebn0_db: NDArray[np.floating] | Sequence[float],
    n_bits: int = 100_000,
    fec: FEC | None = None,
    seed: int | None = None,
    *,
    min_errors: int = 0,
    max_bits: int | None = None,
    batch_bits: int | None = None,
) -> dict[str, NDArray[np.float64]]:
    """Estimate BER and return only the curves.

    With the defaults this reproduces the historical fixed-``n_bits``
    behaviour: one batch, no early stopping.  Set ``min_errors`` and
    ``max_bits`` for an adaptive run, or call :func:`ber_sweep_detailed` for
    the counts and confidence intervals as well.
    """
    if n_bits <= 0:
        raise ValueError("n_bits must be positive")
    target_max = n_bits if max_bits is None else max_bits
    target_batch = target_max if batch_bits is None else batch_bits
    detailed = ber_sweep_detailed(
        modulations,
        ebn0_db,
        fec=fec,
        seed=seed,
        min_errors=min_errors,
        max_bits=target_max,
        batch_bits=target_batch,
    )
    return {name: curve.ber for name, curve in detailed.items()}
