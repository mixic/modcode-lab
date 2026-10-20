# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Block Rayleigh flat fading with known channel state information."""

from __future__ import annotations

from numbers import Integral

import numpy as np
from numpy.typing import NDArray

from modcode.channel.awgn import AWGNChannel
from modcode.channel.base import Channel


class BlockRayleighChannel(Channel):
    """Flat Rayleigh fading, constant over blocks of symbols.

    A new complex Gaussian coefficient is drawn for each block.  The receiver
    is assumed to know the coefficient and applies zero-forcing equalization
    when ``equalize`` is true.  The equalized noise is coloured by
    ``1/|h|^2``; for soft decoding the caller should use the per-symbol
    variance, not the scalar input variance.

    ``last_coefficients`` stores the coefficients from the most recent call.
    It is intended for diagnostics and plotting, not for concurrent use.
    """

    def __init__(self, block_length: int = 1, *, equalize: bool = True) -> None:
        if (
            isinstance(block_length, bool)
            or not isinstance(block_length, Integral)
            or block_length < 1
        ):
            raise ValueError("block_length must be a positive integer")
        self.block_length = int(block_length)
        self.equalize = bool(equalize)
        self.last_coefficients: NDArray[np.complex128] | None = None

    @staticmethod
    def noise_variance(
        ebn0_db: float, bits_per_symbol: int, code_rate: float = 1.0
    ) -> float:
        return AWGNChannel.noise_variance(ebn0_db, bits_per_symbol, code_rate)

    def coefficients(
        self, size: int, rng: np.random.Generator
    ) -> NDArray[np.complex128]:
        """Return one Rayleigh coefficient per block, repeated per symbol."""
        if size <= 0:
            raise ValueError("size must be positive")
        blocks = (size + self.block_length - 1) // self.block_length
        values = (
            rng.standard_normal(blocks) + 1j * rng.standard_normal(blocks)
        ) / np.sqrt(2.0)
        return np.repeat(values, self.block_length)[:size]

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
        coefficients = self.coefficients(values.size, generator)
        noise = np.sqrt(variance) * (
            generator.standard_normal(values.size)
            + 1j * generator.standard_normal(values.size)
        )
        received = coefficients * values + noise
        self.last_coefficients = coefficients
        if self.equalize:
            return received / coefficients
        return received


__all__ = ["BlockRayleighChannel"]
