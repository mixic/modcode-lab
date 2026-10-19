# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Concentric-ring APSK constellations with explicit ring phase actions."""

from collections.abc import Sequence
from numbers import Integral

import numpy as np
from numpy.typing import NDArray

from modcode.modulation.base import BitArray, Modulation, SymbolArray


class APSK(Modulation):
    """Amplitude phase-shift keying with configurable ring populations.

    Each ring's equally spaced phases form an orbit under a cyclic rotation
    action. The complete multi-ring constellation is not asserted to be a
    group under complex multiplication.
    """

    def __init__(
        self,
        ring_sizes: Sequence[int] = (4, 12),
        radii: Sequence[float] = (1.0, 2.0),
        phase_offsets: Sequence[float] | None = None,
    ) -> None:
        sizes = tuple(ring_sizes)
        radial_levels = np.asarray(radii, dtype=np.float64)
        if not sizes or len(sizes) != radial_levels.size:
            raise ValueError("ring_sizes and radii must have the same nonzero length")
        if any(
            isinstance(size, bool) or not isinstance(size, Integral) or size < 1
            for size in sizes
        ):
            raise ValueError("each ring size must be a positive integer")
        sizes = tuple(int(size) for size in sizes)
        order = sum(sizes)
        if order & (order - 1):
            raise ValueError("the total constellation size must be a power of two")
        if not np.all(np.isfinite(radial_levels)) or np.any(radial_levels <= 0):
            raise ValueError("ring radii must be finite and positive")
        if np.any(np.diff(radial_levels) <= 0):
            raise ValueError("ring radii must be strictly increasing")
        if phase_offsets is None:
            offsets = np.zeros(len(sizes), dtype=np.float64)
        else:
            offsets = np.asarray(phase_offsets, dtype=np.float64)
            if offsets.shape != (len(sizes),) or not np.all(np.isfinite(offsets)):
                raise ValueError("phase_offsets must contain one finite angle per ring")

        self.ring_sizes = sizes
        self.radii = radial_levels
        self.phase_offsets = offsets
        self.order = order
        self._bits_per_symbol = order.bit_length() - 1
        points = np.concatenate(
            [
                radius
                * np.exp(
                    1j
                    * (offset + 2.0 * np.pi * np.arange(size, dtype=np.float64) / size)
                )
                for size, radius, offset in zip(
                    sizes, radial_levels, offsets, strict=True
                )
            ]
        )
        points /= np.sqrt(np.mean(np.abs(points) ** 2))
        self._points = points.astype(np.complex128)
        self._labels = self._build_labels()
        self._label_bits = np.asarray(list(self._labels), dtype=np.uint8)

    @property
    def bits_per_symbol(self) -> int:
        return self._bits_per_symbol

    @property
    def constellation(self) -> dict[tuple[int, ...], complex]:
        """Return Gray labels mapped to normalized complex constellation points."""
        return self._labels.copy()

    def modulate(self, bits: NDArray[np.integer]) -> SymbolArray:
        values = self._validated_bits(bits)
        if values.size % self.bits_per_symbol:
            raise ValueError(
                f"{self.order}-APSK requires a multiple of "
                f"{self.bits_per_symbol} input bits"
            )
        labels = values.reshape(-1, self.bits_per_symbol)
        gray_indexes = np.zeros(labels.shape[0], dtype=np.intp)
        for bit_index in range(self.bits_per_symbol):
            gray_indexes = (gray_indexes << 1) | labels[:, bit_index]
        point_indexes = self._gray_to_binary(gray_indexes)
        return self._points[point_indexes]

    def demodulate(self, symbols: NDArray[np.complexfloating]) -> BitArray:
        values = self._validated_symbols(symbols)
        if values.size == 0:
            return np.empty(0, dtype=np.uint8)
        distances = np.abs(values[:, None] - self._points[None, :]) ** 2
        nearest = np.argmin(distances, axis=1)
        return self._label_bits[nearest].reshape(-1)

    def ring_points(self, ring_index: int) -> NDArray[np.complex128]:
        """Return the points on one ring in increasing phase order."""
        if isinstance(ring_index, bool) or not isinstance(ring_index, Integral):
            raise ValueError("ring_index must be an integer")
        if not 0 <= ring_index < len(self.ring_sizes):
            raise IndexError("ring_index is outside the constellation")
        start = sum(self.ring_sizes[:ring_index])
        stop = start + self.ring_sizes[ring_index]
        return self._points[start:stop].copy()

    def rotate_ring(self, ring_index: int, steps: int = 1) -> NDArray[np.complex128]:
        """Apply the ring's cyclic action by an integer number of phase steps."""
        if isinstance(steps, bool) or not isinstance(steps, Integral):
            raise ValueError("steps must be an integer")
        return np.roll(self.ring_points(ring_index), -int(steps))

    def _build_labels(self) -> dict[tuple[int, ...], complex]:
        labels: dict[tuple[int, ...], complex] = {}
        for point_index, point in enumerate(self._points):
            gray_index = point_index ^ (point_index >> 1)
            label = tuple(
                (gray_index >> shift) & 1
                for shift in range(self.bits_per_symbol - 1, -1, -1)
            )
            labels[label] = complex(point)
        return labels

    @staticmethod
    def _gray_to_binary(gray: NDArray[np.intp]) -> NDArray[np.intp]:
        binary = gray.copy()
        shift = gray >> 1
        while np.any(shift):
            binary ^= shift
            shift >>= 1
        return binary
