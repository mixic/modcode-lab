# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Binary phase-shift keying."""

import numpy as np
from numpy.typing import NDArray

from modcode.modulation.base import BitArray, Modulation, SymbolArray


class BPSK(Modulation):
    """Unit-energy BPSK with bit 0 mapped to +1 and bit 1 to -1."""

    @property
    def bits_per_symbol(self) -> int:
        return 1

    @property
    def constellation(self) -> dict[tuple[int, ...], complex]:
        return {(0,): 1.0 + 0.0j, (1,): -1.0 + 0.0j}

    def modulate(self, bits: NDArray[np.integer]) -> SymbolArray:
        values = self._validated_bits(bits)
        return (1.0 - 2.0 * values).astype(np.complex128)

    def demodulate(self, symbols: NDArray[np.complexfloating]) -> BitArray:
        values = self._validated_symbols(symbols)
        return (values.real < 0).astype(np.uint8)
