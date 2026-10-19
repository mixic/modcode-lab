# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Nonnegative intensity modulation formats for optical IM/DD experiments."""

import numpy as np
from numpy.typing import NDArray

from modcode.modulation.base import BitArray


class IntensityPAM:
    """Natural-binary M-PAM with nonnegative, unit-mean intensity levels."""

    def __init__(self, order: int = 4) -> None:
        if isinstance(order, bool) or not isinstance(order, int):
            raise ValueError("PAM order must be an integer")
        if order < 2 or order & (order - 1):
            raise ValueError("PAM order must be a power of two >= 2")
        self.order = order
        self.bits_per_symbol = order.bit_length() - 1
        self.levels = np.linspace(0.0, 2.0, order, dtype=np.float64)
        self._validate_levels(self.levels)

    @property
    def constellation(self) -> NDArray[np.float64]:
        """Return a copy of the nonnegative optical intensity levels."""
        return self.levels.copy()

    def modulate(self, bits: NDArray[np.integer]) -> NDArray[np.float64]:
        values = _validate_bits(bits)
        if values.size % self.bits_per_symbol:
            raise ValueError(
                f"{self.order}-PAM requires a multiple of "
                f"{self.bits_per_symbol} input bits"
            )
        labels = values.reshape(-1, self.bits_per_symbol)
        gray_indexes = np.zeros(labels.shape[0], dtype=np.intp)
        for bit_index in range(self.bits_per_symbol):
            gray_indexes = (gray_indexes << 1) | labels[:, bit_index]
        indexes = gray_indexes.copy()
        shift = gray_indexes >> 1
        while np.any(shift):
            indexes ^= shift
            shift >>= 1
        return self.levels[indexes]

    def demodulate(self, symbols: NDArray[np.floating]) -> BitArray:
        values = _validate_intensities(symbols)
        indexes = np.argmin(np.abs(values[:, None] - self.levels[None, :]), axis=1)
        gray_indexes = indexes ^ (indexes >> 1)
        bits = np.empty((values.size, self.bits_per_symbol), dtype=np.uint8)
        for bit_index in range(self.bits_per_symbol):
            shift = self.bits_per_symbol - bit_index - 1
            bits[:, bit_index] = (gray_indexes >> shift) & 1
        return bits.reshape(-1)

    @staticmethod
    def _validate_levels(levels: NDArray[np.float64]) -> None:
        if levels.ndim != 1 or levels.size < 2:
            raise ValueError("intensity levels must be a one-dimensional array")
        if not np.all(np.isfinite(levels)) or np.any(levels < 0):
            raise ValueError("intensity levels must be finite and nonnegative")
        if np.any(np.diff(levels) <= 0):
            raise ValueError("intensity levels must be strictly increasing")
        if not np.isclose(levels.mean(), 1.0):
            raise ValueError("intensity levels must have unit mean power")


class OOK(IntensityPAM):
    """On-off keying: the two unit-mean intensity levels are 0 and 2."""

    def __init__(self) -> None:
        super().__init__(order=2)


class NonUniformPAM(IntensityPAM):
    """M-PAM with explicitly chosen increasing, unit-mean intensity levels."""

    def __init__(self, levels: NDArray[np.floating]) -> None:
        values = np.asarray(levels, dtype=np.float64)
        if values.size < 2 or values.size & (values.size - 1):
            raise ValueError("the number of PAM levels must be a power of two >= 2")
        self.order = int(values.size)
        self.bits_per_symbol = self.order.bit_length() - 1
        self.levels = values.copy()
        self._validate_levels(self.levels)


def _validate_bits(bits: NDArray[np.integer]) -> BitArray:
    values = np.asarray(bits)
    if values.ndim != 1:
        raise ValueError("bits must be a one-dimensional array")
    if not np.all((values == 0) | (values == 1)):
        raise ValueError("bits must contain only 0 and 1")
    return values.astype(np.uint8, copy=False)


def _validate_intensities(
    symbols: NDArray[np.floating],
) -> NDArray[np.float64]:
    values = np.asarray(symbols, dtype=np.float64)
    if values.ndim != 1:
        raise ValueError("intensity symbols must be one-dimensional")
    if not np.all(np.isfinite(values)):
        raise ValueError("intensity symbols must be finite")
    return values
