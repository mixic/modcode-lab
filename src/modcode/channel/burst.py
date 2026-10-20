# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Synthetic burst channel for interleaver demonstrations."""

from __future__ import annotations

from numbers import Integral

import numpy as np
from numpy.typing import NDArray

from modcode.channel.awgn import AWGNChannel
from modcode.channel.base import Channel


class BurstChannel(Channel):
    """AWGN with an occasional contiguous block of erased symbols.

    The burst is a modelling device, not a physical fading model: it makes
    the interleaver's effect visible by putting several symbol errors next to
    each other.  ``burst_probability`` is the long-run fraction of symbols
    erased, so the number of bursts scales with the block length.
    ``last_burst`` stores the most recent ``(start, stop)`` interval and
    ``last_bursts`` stores all of them for diagnostics.
    """

    def __init__(
        self,
        burst_length: int = 8,
        burst_probability: float = 0.05,
    ) -> None:
        if (
            isinstance(burst_length, bool)
            or not isinstance(burst_length, Integral)
            or burst_length < 1
        ):
            raise ValueError("burst_length must be a positive integer")
        probability = float(burst_probability)
        if not 0.0 <= probability <= 1.0:
            raise ValueError("burst_probability must be in [0, 1]")
        self.burst_length = int(burst_length)
        self.burst_probability = probability
        self.last_burst: tuple[int, int] | None = None
        self.last_bursts: list[tuple[int, int]] = []

    @staticmethod
    def noise_variance(
        ebn0_db: float, bits_per_symbol: int, code_rate: float = 1.0
    ) -> float:
        return AWGNChannel.noise_variance(ebn0_db, bits_per_symbol, code_rate)

    def transmit(
        self,
        symbols: NDArray[np.complexfloating],
        noise_variance: float,
        rng: np.random.Generator | None = None,
    ) -> NDArray[np.complex128]:
        values = np.asarray(symbols, dtype=np.complex128)
        if values.ndim != 1:
            raise ValueError("symbols must be a one-dimensional array")
        if not np.all(np.isfinite(values)):
            raise ValueError("symbols must contain only finite values")
        variance = float(noise_variance)
        if not np.isfinite(variance) or variance < 0.0:
            raise ValueError("noise_variance must be finite and non-negative")
        generator = rng if rng is not None else np.random.default_rng()
        noise = np.sqrt(variance) * (
            generator.standard_normal(values.size)
            + 1j * generator.standard_normal(values.size)
        )
        received = values + noise
        self.last_burst = None
        self.last_bursts = []
        if self.burst_probability > 0.0:
            expected = values.size * self.burst_probability / self.burst_length
            bursts = int(generator.poisson(expected))
            for _ in range(bursts):
                start = int(generator.integers(0, values.size))
                stop = min(values.size, start + self.burst_length)
                received[start:stop] = 0.0
                self.last_burst = (start, stop)
                self.last_bursts.append((start, stop))
        return received


__all__ = ["BurstChannel"]
