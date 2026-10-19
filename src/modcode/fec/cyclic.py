# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Small utilities for exploring cyclic binary codes."""

import numpy as np
from numpy.typing import NDArray

from modcode.fec.base import FEC, BitArray


def cyclic_shift(bits: NDArray[np.integer], shifts: int = 1) -> BitArray:
    """Rotate a one-dimensional binary vector to the left."""
    values = FEC.validated_bits(bits)
    if isinstance(shifts, bool) or not isinstance(shifts, int):
        raise ValueError("shifts must be an integer")
    if values.size == 0:
        return values.copy()
    return np.roll(values, -(shifts % values.size))
