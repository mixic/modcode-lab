# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Build inspectable binary quasi-cyclic LDPC parity-check matrices."""

from numbers import Integral

import numpy as np
from numpy.typing import NDArray


def circulant_identity(size: int, shift: int) -> NDArray[np.uint8]:
    """Return I cyclically shifted right by ``shift`` columns modulo ``size``."""
    if isinstance(size, bool) or not isinstance(size, Integral) or size < 1:
        raise ValueError("size must be a positive integer")
    if isinstance(shift, bool) or not isinstance(shift, Integral):
        raise ValueError("shift must be an integer")
    return np.roll(np.eye(int(size), dtype=np.uint8), int(shift) % int(size), axis=1)


def expand_qc_base_matrix(
    base_matrix: NDArray[np.integer],
    lifting_factor: int,
) -> NDArray[np.uint8]:
    """Expand a shift-value base matrix to its binary QC parity-check matrix.

    A base-matrix entry of ``-1`` denotes a zero block; a nonnegative entry
    ``s`` denotes the circulant permutation matrix ``I`` shifted right by
    ``s`` columns modulo the lifting factor.
    """
    if (
        isinstance(lifting_factor, bool)
        or not isinstance(lifting_factor, Integral)
        or lifting_factor < 1
    ):
        raise ValueError("lifting_factor must be a positive integer")
    lifting_factor = int(lifting_factor)
    base = np.asarray(base_matrix)
    if base.ndim != 2 or base.size == 0 or not np.issubdtype(base.dtype, np.integer):
        raise ValueError("base_matrix must be a nonempty two-dimensional integer array")
    if np.any(base < -1):
        raise ValueError("base-matrix entries must be -1 or nonnegative shifts")
    if np.any(base >= lifting_factor):
        raise ValueError("shift entries must be smaller than lifting_factor")

    rows, columns = base.shape
    expanded = np.zeros(
        (rows * lifting_factor, columns * lifting_factor), dtype=np.uint8
    )
    identity = np.eye(lifting_factor, dtype=np.uint8)
    for row in range(rows):
        row_slice = slice(row * lifting_factor, (row + 1) * lifting_factor)
        for column in range(columns):
            shift = int(base[row, column])
            if shift >= 0:
                column_slice = slice(
                    column * lifting_factor, (column + 1) * lifting_factor
                )
                expanded[row_slice, column_slice] = np.roll(identity, shift, axis=1)
    return expanded
