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
    """

    def __init__(self, order: int = 16) -> None:
        if isinstance(order, bool) or not isinstance(order, Integral):
            raise ValueError("QAM order must be an integer")
        order = int(order)
        side = math.isqrt(order)
        if order < 4 or side * side != order or side & (side - 1):
            raise ValueError("QAM order must be a square power of two >= 4")
        self.order = order
        self._axis_bits = int(math.log2(side))
        self._bits_per_symbol = 2 * self._axis_bits
        self._scale = math.sqrt(2.0 * (order - 1) / 3.0)
        self._labels = self._build_labels()
        self._points = np.asarray(list(self._labels.values()), dtype=np.complex128)
        self._label_bits = np.asarray(list(self._labels.keys()), dtype=np.uint8)

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
        indexes = self._bits_to_indexes(labels)
        return self._points[indexes]

    def demodulate(self, symbols: NDArray[np.complexfloating]) -> BitArray:
        values = self._validated_symbols(symbols)
        if values.size == 0:
            return np.empty(0, dtype=np.uint8)
        distances = np.abs(values[:, None] - self._points[None, :]) ** 2
        nearest = np.argmin(distances, axis=1)
        return self._label_bits[nearest].reshape(-1)

    def _build_labels(self) -> dict[tuple[int, ...], complex]:
        labels: dict[tuple[int, ...], complex] = {}
        side = math.isqrt(self.order)
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

    def _bits_to_indexes(self, labels: NDArray[np.uint8]) -> NDArray[np.intp]:
        label_indexes = np.zeros(labels.shape[0], dtype=np.intp)
        for index, point_bits in enumerate(self._label_bits):
            matches = np.all(labels == point_bits, axis=1)
            label_indexes[matches] = index
        return label_indexes

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
