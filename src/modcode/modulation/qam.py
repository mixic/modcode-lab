# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Gray-mapped square QAM constellations."""

import math
from numbers import Integral

import numpy as np
from numpy.typing import NDArray

from modcode.modulation.base import BitArray, Modulation, SymbolArray


class QAM(Modulation):
    """Square Gray-coded QAM with unit average symbol energy.

    The I-axis label comes first, followed by the Q-axis label. Each axis
    independently uses natural-order Gray labels.

    Hard slicing is per axis: the received coordinate is scaled to the odd
    integer level grid and rounded, so the cost is O(N) in the number of
    symbols and no N-by-M distance matrix is allocated.  The soft LLR method
    uses the same per-axis decomposition.
    """

    def __init__(self, order: int = 16) -> None:
        if isinstance(order, bool) or not isinstance(order, Integral):
            raise ValueError("QAM order must be an integer")
        order = int(order)
        side = math.isqrt(order)
        if order < 4 or side * side != order or side & (side - 1):
            raise ValueError("QAM order must be a square power of two >= 4")
        self.order = order
        self._side = side
        self._axis_bits = int(math.log2(side))
        self._bits_per_symbol = 2 * self._axis_bits
        self._scale = math.sqrt(2.0 * (order - 1) / 3.0)
        self._labels = self._build_labels()
        self._axis_levels = (
            2 * np.arange(side, dtype=np.float64) - side + 1
        )
        self._axis_gray = self._binary_to_gray(
            np.arange(side, dtype=np.intp)
        )
        self._axis_bit_labels = self._index_to_bits(self._axis_gray)

    @property
    def bits_per_symbol(self) -> int:
        return self._bits_per_symbol

    @property
    def constellation(self) -> dict[tuple[int, ...], complex]:
        return self._labels.copy()

    def modulate(self, bits: NDArray[np.integer]) -> SymbolArray:
        values = self._validated_bits(bits)
        if values.size % self.bits_per_symbol:
            raise ValueError(
                f"{self.order}-QAM requires a multiple of "
                f"{self.bits_per_symbol} input bits"
            )
        labels = values.reshape(-1, self.bits_per_symbol)
        i_index = self._axis_index(labels[:, : self._axis_bits])
        q_index = self._axis_index(labels[:, self._axis_bits :])
        i_level = self._axis_levels[i_index] / self._scale
        q_level = self._axis_levels[q_index] / self._scale
        return (i_level + 1j * q_level).astype(np.complex128)

    def demodulate(self, symbols: NDArray[np.complexfloating]) -> BitArray:
        values = self._validated_symbols(symbols)
        if values.size == 0:
            return np.empty(0, dtype=np.uint8)
        i_index = self._slice_axis(values.real)
        q_index = self._slice_axis(values.imag)
        i_bits = self._index_to_bits(
            self._binary_to_gray(i_index)
        )
        q_bits = self._index_to_bits(
            self._binary_to_gray(q_index)
        )
        return np.concatenate([i_bits, q_bits], axis=1).reshape(-1)

    def soft_demodulate(
        self,
        symbols: NDArray[np.complexfloating],
        noise_variance: float,
    ) -> NDArray[np.float64]:
        """Return per-axis max-log LLRs without an N-by-M distance matrix."""
        values = self._validated_symbols(symbols)
        variance = float(noise_variance)
        if not np.isfinite(variance) or variance <= 0.0:
            raise ValueError("noise_variance must be finite and positive")
        if values.size == 0:
            return np.empty(0, dtype=np.float64)

        level_variance = variance * self._scale**2
        llrs = np.empty(
            (values.size, self.bits_per_symbol), dtype=np.float64
        )
        for axis, coordinate in enumerate((values.real, values.imag)):
            levels = coordinate * self._scale
            distances = (
                np.abs(levels[:, None] - self._axis_levels[None, :]) ** 2
            )
            for bit in range(self._axis_bits):
                zero = self._axis_bit_labels[:, bit] == 0
                nearest_zero = np.min(distances[:, zero], axis=1)
                nearest_one = np.min(distances[:, ~zero], axis=1)
                llrs[:, axis * self._axis_bits + bit] = (
                    nearest_one - nearest_zero
                ) / (2.0 * level_variance)
        return llrs.reshape(-1)

    def _slice_axis(self, coordinate: NDArray[np.float64]) -> NDArray[np.intp]:
        levels = coordinate * self._scale
        index = np.rint((levels + self._side - 1) / 2.0).astype(np.intp)
        return np.clip(index, 0, self._side - 1)

    def _axis_index(self, labels: NDArray[np.uint8]) -> NDArray[np.intp]:
        weights = (
            1 << np.arange(self._axis_bits - 1, -1, -1)
        ).astype(np.intp)
        gray = labels.astype(np.intp) @ weights
        return self._gray_to_binary_array(gray)

    def _gray_to_binary_array(
        self, gray: NDArray[np.intp]
    ) -> NDArray[np.intp]:
        value = gray.copy()
        shift = 1
        while shift < (1 << self._axis_bits):
            value ^= value >> shift
            shift <<= 1
        return value

    def _index_to_bits(self, indices: NDArray[np.intp]) -> NDArray[np.uint8]:
        shifts = np.arange(self._axis_bits - 1, -1, -1)
        return ((indices[:, None] >> shifts) & 1).astype(np.uint8)

    @staticmethod
    def _binary_to_gray(indices: NDArray[np.intp]) -> NDArray[np.intp]:
        return indices ^ (indices >> 1)

    def _build_labels(self) -> dict[tuple[int, ...], complex]:
        labels: dict[tuple[int, ...], complex] = {}
        side = self._side
        scale = self._scale
        for integer in range(self.order):
            bits = tuple(
                (integer >> shift) & 1
                for shift in range(self.bits_per_symbol - 1, -1, -1)
            )
            i_index = self._gray_to_binary(
                self._bits_to_integer(bits[: self._axis_bits])
            )
            q_index = self._gray_to_binary(
                self._bits_to_integer(bits[self._axis_bits :])
            )
            i_level = 2 * i_index - side + 1
            q_level = 2 * q_index - side + 1
            labels[bits] = complex(i_level / scale, q_level / scale)
        return labels

    @staticmethod
    def _bits_to_integer(bits: tuple[int, ...]) -> int:
        value = 0
        for bit in bits:
            value = (value << 1) | bit
        return value

    @staticmethod
    def _gray_to_binary(gray: int) -> int:
        value = gray
        while gray:
            gray >>= 1
            value ^= gray
        return value
