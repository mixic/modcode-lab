# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Multi-dimensional lattice modulation with spherical cuts.

This is a standalone appendix, not a drop-in ``Modulation`` implementation:
a D4 symbol is two complex components and an E8 symbol is four, while the
core ``Link`` pipeline carries one complex value per symbol.  The class
therefore exposes a vector ``modulate``/``demodulate`` API and is simulated by
:mod:`modcode.metrics.lattice`, not by ``run_link``.

D4 is the set of integer points in R^4 with even coordinate sum.  E8 is the
union of the integer points in R^8 with even sum and the half-integer points
with even sum.  Points are enumerated in order of squared norm and the M
lowest-norm points form the spherical cut.  The class normalizes to unit
average *symbol* energy; :meth:`per_dimension_points` rescales to unit energy
per real dimension for a fair comparison with 2D QAM.
"""

from __future__ import annotations

import itertools
from functools import lru_cache
from numbers import Integral

import numpy as np
from numpy.typing import ArrayLike, NDArray


@lru_cache(maxsize=16)
def _d4_candidates(max_shell: int) -> NDArray[np.float64]:
    points = [
        point
        for point in itertools.product(
            range(-max_shell, max_shell + 1), repeat=4
        )
        if sum(point) % 2 == 0 and any(point)
    ]
    array = np.asarray(points, dtype=np.float64)
    return array[np.argsort(np.sum(array**2, axis=1), kind="stable")]


@lru_cache(maxsize=8)
def _e8_candidates(max_shell: int) -> NDArray[np.float64]:
    integer = [
        point
        for point in itertools.product(
            range(-max_shell, max_shell + 1), repeat=8
        )
        if sum(point) % 2 == 0 and any(point)
    ]
    half = [
        point
        for point in itertools.product(
            np.arange(-max_shell, max_shell) + 0.5, repeat=8
        )
        if sum(point) % 2 == 0
    ]
    array = np.asarray(integer + half, dtype=np.float64)
    return array[np.argsort(np.sum(array**2, axis=1), kind="stable")]


def _candidates(kind: str, count: int) -> NDArray[np.float64]:
    for max_shell in range(1, 7):
        points = (
            _d4_candidates(max_shell)
            if kind == "D4"
            else _e8_candidates(max_shell)
        )
        if points.shape[0] >= count:
            return points
    raise ValueError(f"{kind} candidate enumeration is too small for M={count}")


class LatticeModulation:
    """Spherical-cut D4 or E8 constellation.

    Parameters
    ----------
    kind
        ``"D4"`` (four real dimensions) or ``"E8"`` (eight).
    m_bits
        Information bits per multi-dimensional symbol.
    labeler
        ``"lex"`` assigns lexicographic bit patterns; ``"gray"`` uses a
        nearest-neighbour walk and the standard binary-reflected Gray code.
        The Gray option is a heuristic, not a proven optimal labeling.
    """

    def __init__(
        self,
        kind: str,
        m_bits: int,
        labeler: str = "lex",
    ) -> None:
        if kind not in {"D4", "E8"}:
            raise ValueError("kind must be 'D4' or 'E8'")
        if (
            isinstance(m_bits, bool)
            or not isinstance(m_bits, Integral)
            or m_bits < 1
        ):
            raise ValueError("m_bits must be a positive integer")
        self.kind = kind
        self.K = 4 if kind == "D4" else 8
        self.M = 1 << int(m_bits)
        self.bits_per_symbol = int(m_bits)
        self.dimensions = self.K
        self.complex_dimensions = self.K // 2

        points = _candidates(kind, self.M)[: self.M].copy()
        self.raw_points = points
        energy = float(np.mean(np.sum(points**2, axis=1)))
        self._unit_symbol_scale = 1.0 / np.sqrt(energy)
        self.points = points * self._unit_symbol_scale

        self.constellation = (
            self.points[:, 0::2] + 1j * self.points[:, 1::2]
        )
        self.bit_patterns = self._labels(labeler)
        self.labeler = labeler
        self._pattern_to_index = {
            tuple(pattern): index
            for index, pattern in enumerate(self.bit_patterns)
        }

    def _labels(self, labeler: str) -> NDArray[np.uint8]:
        if labeler == "lex":
            return np.asarray(
                [
                    [
                        (index >> (self.bits_per_symbol - 1 - bit)) & 1
                        for bit in range(self.bits_per_symbol)
                    ]
                    for index in range(self.M)
                ],
                dtype=np.uint8,
            )
        if labeler != "gray":
            raise ValueError("labeler must be 'lex' or 'gray'")

        visited = np.zeros(self.M, dtype=bool)
        order = np.empty(self.M, dtype=np.intp)
        current = 0
        visited[current] = True
        order[0] = current
        for step in range(1, self.M):
            distances = np.sum(
                (self.points - self.points[current]) ** 2, axis=1
            )
            distances[visited] = np.inf
            current = int(np.argmin(distances))
            visited[current] = True
            order[step] = current

        patterns = np.empty(
            (self.M, self.bits_per_symbol), dtype=np.uint8
        )
        for step, point_index in enumerate(order):
            gray = step ^ (step >> 1)
            patterns[point_index] = [
                (gray >> (self.bits_per_symbol - 1 - bit)) & 1
                for bit in range(self.bits_per_symbol)
            ]
        return patterns

    def modulate(self, bits: NDArray[np.integer]) -> NDArray[np.complex128]:
        values = np.asarray(bits)
        if values.ndim != 1:
            raise ValueError("bits must be a one-dimensional array")
        if not np.all((values == 0) | (values == 1)):
            raise ValueError("bits must contain only 0 and 1")
        if values.size == 0 or values.size % self.bits_per_symbol:
            raise ValueError(
                "bit count must be a nonempty multiple of bits_per_symbol"
            )
        patterns = values.reshape(-1, self.bits_per_symbol).astype(np.uint8)
        indices = np.asarray(
            [
                self._pattern_to_index[tuple(pattern)]
                for pattern in patterns
            ],
            dtype=np.intp,
        )
        return self.constellation[indices]

    def demodulate(self, symbols: ArrayLike) -> NDArray[np.uint8]:
        values = np.asarray(symbols, dtype=np.complex128)
        if values.ndim != 2 or values.shape[1] != self.complex_dimensions:
            raise ValueError(
                "symbols must have shape (N, complex_dimensions)"
            )
        if values.shape[0] == 0 or not np.all(np.isfinite(values)):
            raise ValueError("symbols must be nonempty and finite")
        real = np.empty((values.shape[0], self.K), dtype=np.float64)
        real[:, 0::2] = values.real
        real[:, 1::2] = values.imag
        distances = np.sum(
            (real[:, None, :] - self.points[None, :, :]) ** 2, axis=2
        )
        indices = np.argmin(distances, axis=1)
        return self.bit_patterns[indices].reshape(-1)

    def average_energy(self) -> float:
        """Return the average symbol energy (normalized to 1)."""
        return float(np.mean(np.sum(self.points**2, axis=1)))

    def per_dimension_points(
        self, energy_per_dimension: float = 1.0
    ) -> NDArray[np.float64]:
        """Return points rescaled to a chosen energy per real dimension."""
        if not np.isfinite(energy_per_dimension) or energy_per_dimension <= 0:
            raise ValueError("energy_per_dimension must be positive")
        return self.points * np.sqrt(self.K * energy_per_dimension)

    def min_distance(
        self, energy_per_dimension: float | None = None
    ) -> float:
        """Return the minimum Euclidean distance.

        With ``energy_per_dimension=None`` the points are the unit-average-
        symbol-energy points used by :meth:`modulate`.  Supply
        ``1.0`` for a fair comparison with 2D QAM at equal energy per real
        dimension.
        """
        points = self.points
        if energy_per_dimension is not None:
            points = self.per_dimension_points(energy_per_dimension)
        distances = np.sum(
            (points[:, None, :] - points[None, :, :]) ** 2, axis=2
        )
        np.fill_diagonal(distances, np.inf)
        return float(np.sqrt(np.min(distances)))

    @property
    def bits_per_real_dimension(self) -> float:
        return self.bits_per_symbol / self.K


__all__ = ["LatticeModulation"]
