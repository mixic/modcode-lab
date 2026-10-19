# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Gray-mapped, unit-energy quadrature phase-shift keying."""

import numpy as np
from numpy.typing import NDArray

from modcode.modulation.base import BitArray, Modulation, SymbolArray


class QPSK(Modulation):
    """QPSK mapping each bit to the sign of the I or Q component."""

    @property
    def bits_per_symbol(self) -> int:
        return 2

    @property
    def constellation(self) -> dict[tuple[int, ...], complex]:
        scale = 1.0 / np.sqrt(2.0)
        return {
            (i, q): complex((1.0 - 2.0 * i) * scale, (1.0 - 2.0 * q) * scale)
            for i in (0, 1)
            for q in (0, 1)
        }

    def modulate(self, bits: NDArray[np.integer]) -> SymbolArray:
        values = self._validated_bits(bits)
        if values.size % self.bits_per_symbol:
            raise ValueError("QPSK requires a multiple of 2 input bits")
        pairs = values.reshape(-1, 2)
        scale = 1.0 / np.sqrt(2.0)
        return (
            (1.0 - 2.0 * pairs[:, 0]) * scale + 1j * (1.0 - 2.0 * pairs[:, 1]) * scale
        ).astype(np.complex128)

    def demodulate(self, symbols: NDArray[np.complexfloating]) -> BitArray:
        values = self._validated_symbols(symbols)
        bits = np.empty((values.size, 2), dtype=np.uint8)
        bits[:, 0] = values.real < 0
        bits[:, 1] = values.imag < 0
        return bits.reshape(-1)
