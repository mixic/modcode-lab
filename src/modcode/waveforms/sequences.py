# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""PSK constellations and Zadoff–Chu sequences.

The two constructions are often described as “cyclic,” but they use different
cyclic structures. With zero phase offset, the order-``M`` PSK points are the
``M``-th roots of unity and form the cyclic group ``C_M`` under complex
multiplication. A Zadoff–Chu sequence is a quadratic-phase sequence whose
periodic autocorrelation vanishes at every nonzero shift; its “root” is a
phase index, not a primitive root.

Autocorrelation functions here are normalized by the sequence length, so a
perfect sequence has peak one and zero sidelobes.
"""

from __future__ import annotations

from math import gcd, pi

import numpy as np
from numpy.typing import ArrayLike, NDArray


def _positive_integer(value: int, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, (int, np.integer)):
        raise ValueError(f"{name} must be an integer")
    number = int(value)
    if number < 1:
        raise ValueError(f"{name} must be positive")
    return number


def _finite_vector(values: ArrayLike, name: str) -> NDArray[np.complex128]:
    array = np.asarray(values, dtype=np.complex128)
    if array.ndim != 1:
        raise ValueError(f"{name} must be a one-dimensional array")
    if array.size == 0:
        raise ValueError(f"{name} must be nonempty")
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must contain only finite values")
    return array


def psk_constellation(
    order: int, phase_offset: float = 0.0
) -> NDArray[np.complex128]:
    """Return the ``order``-PSK points in increasing phase order.

    The points are ``exp(j*(phase_offset + 2*pi*k/order))``. With
    ``phase_offset=0`` they are the ``order``-th roots of unity and form the
    cyclic group ``C_order`` under multiplication. A nonzero global offset
    makes them a coset of that group, not a subgroup.
    """
    order = _positive_integer(order, "order")
    if order < 2:
        raise ValueError("order must be at least 2")
    offset = float(phase_offset)
    if not np.isfinite(offset):
        raise ValueError("phase_offset must be finite")
    indexes = np.arange(order, dtype=np.float64)
    return np.exp(1j * (offset + 2.0 * pi * indexes / order))


def zadoff_chu(length: int, root: int = 1) -> NDArray[np.complex128]:
    """Return the odd-length Zadoff–Chu sequence with the given root index.

    The sequence is ``s_root[n] = exp(-j*pi*root*n*(n+1)/length)`` for
    ``n = 0, ..., length - 1``. ``length`` must be odd and ``root`` must be
    coprime to it; under those conditions the periodic autocorrelation is a
    delta. The root index is a quadratic phase coefficient, not a primitive
    root modulo ``length``.
    """
    length = _positive_integer(length, "length")
    if length < 3 or length % 2 == 0:
        raise ValueError("length must be an odd integer at least 3")
    root = _positive_integer(root, "root")
    if root >= length or gcd(root, length) != 1:
        raise ValueError("root must lie in [1, length) and be coprime to length")
    indexes = np.arange(length, dtype=np.float64)
    return np.exp(-1j * pi * root * indexes * (indexes + 1.0) / length)


def periodic_autocorrelation(
    sequence: ArrayLike, *, normalize: bool = True
) -> NDArray[np.complex128]:
    """Return the periodic autocorrelation for every cyclic lag.

    ``result[lag] = sum_n sequence[n]* sequence[(n + lag) mod N]``. With
    ``normalize=True`` the result is divided by the sequence length, so a
    perfect sequence has ``result[0] == 1``.
    """
    values = _finite_vector(sequence, "sequence")
    result = np.asarray(
        [np.vdot(values, np.roll(values, lag)) for lag in range(values.size)],
        dtype=np.complex128,
    )
    return result / values.size if normalize else result


def cyclic_cross_correlation(
    first: ArrayLike, second: ArrayLike, *, normalize: bool = True
) -> NDArray[np.complex128]:
    """Return the periodic cross-correlation for every cyclic lag.

    ``result[lag] = sum_n first[n]* second[(n + lag) mod N]``. The two
    sequences must have the same length. With ``normalize=True`` the result
    is divided by that length.
    """
    left = _finite_vector(first, "first")
    right = _finite_vector(second, "second")
    if left.size != right.size:
        raise ValueError("sequences must have the same length")
    result = np.asarray(
        [np.vdot(left, np.roll(right, lag)) for lag in range(left.size)],
        dtype=np.complex128,
    )
    return result / left.size if normalize else result
