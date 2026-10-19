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
