# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Shared interface and input checks for digital modulation schemes."""

from abc import ABC, abstractmethod

import numpy as np
from numpy.typing import NDArray

BitArray = NDArray[np.uint8]
SymbolArray = NDArray[np.complex128]


class Modulation(ABC):
    """Map binary data to constellation symbols and back."""

    @property
    @abstractmethod
    def bits_per_symbol(self) -> int:
        """Number of binary digits represented by one symbol."""

    @property
    @abstractmethod
    def constellation(self) -> dict[tuple[int, ...], complex]:
        """Return the mapping from bit labels to complex symbols."""

    @abstractmethod
    def modulate(self, bits: NDArray[np.integer]) -> SymbolArray:
        """Map a one-dimensional bit sequence to complex symbols."""

    @abstractmethod
    def demodulate(self, symbols: NDArray[np.complexfloating]) -> BitArray:
        """Make hard decisions and return the corresponding bits."""

    def soft_demodulate(
        self,
        symbols: NDArray[np.complexfloating],
        noise_variance: float,
    ) -> NDArray[np.float64]:
        """Return max-log bit LLRs for passed symbols.

        The convention is ``L = log P(b=0|y) / P(b=1|y)``: positive values
        favour bit 0.  ``noise_variance`` is the variance per real component
        of the circular complex AWGN, matching
        :meth:`modcode.channel.awgn.AWGNChannel.noise_variance`.

        The implementation uses the constellation and its bit labels, so it
        applies to every modulation with a finite constellation.  It is
        max-log rather than exact log-sum-exp; for BPSK and QPSK the two
        agree because each bit has one constellation point per value.
        """
        values = self._validated_symbols(symbols)
        variance = float(noise_variance)
        if not np.isfinite(variance) or variance <= 0.0:
            raise ValueError("noise_variance must be finite and positive")
        if values.size == 0:
            return np.empty(0, dtype=np.float64)

        points = np.asarray(
            list(self.constellation.values()), dtype=np.complex128
        )
        labels = np.asarray(
            list(self.constellation.keys()), dtype=np.uint8
        )
        if points.ndim != 1 or labels.ndim != 2:
            raise ValueError("constellation must be a finite bit-label mapping")
        if labels.shape != (points.size, self.bits_per_symbol):
            raise ValueError("constellation labels do not match bits_per_symbol")

        distances = np.abs(values[:, None] - points[None, :]) ** 2
        llrs = np.empty(
            (values.size, self.bits_per_symbol), dtype=np.float64
        )
        for bit in range(self.bits_per_symbol):
            zero = labels[:, bit] == 0
            if not np.any(zero) or np.all(zero):
                raise ValueError("each bit position must take both values")
            nearest_zero = np.min(distances[:, zero], axis=1)
            nearest_one = np.min(distances[:, ~zero], axis=1)
            llrs[:, bit] = (nearest_one - nearest_zero) / (2.0 * variance)
        return llrs.reshape(-1)

    def _validated_bits(self, bits: NDArray[np.integer]) -> BitArray:
        values = np.asarray(bits)
        if values.ndim != 1:
            raise ValueError("bits must be a one-dimensional array")
        if not np.all((values == 0) | (values == 1)):
            raise ValueError("bits must contain only 0 and 1")
        return values.astype(np.uint8, copy=False)

    def _validated_symbols(
        self, symbols: NDArray[np.complexfloating]
    ) -> NDArray[np.complex128]:
        values = np.asarray(symbols, dtype=np.complex128)
        if values.ndim != 1:
            raise ValueError("symbols must be a one-dimensional array")
        if not np.all(np.isfinite(values)):
            raise ValueError("symbols must contain only finite values")
        return values
