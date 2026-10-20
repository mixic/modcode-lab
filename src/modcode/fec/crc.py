# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Systematic cyclic-redundancy-check encoding and error detection.

A CRC appends the remainder of ``x**width * m(x)`` divided by a fixed
generator polynomial ``g(x)``. It is a cyclic code in systematic form: the
message bits are unchanged and the parity bits are the remainder. Unlike
Hamming or Reed–Solomon codes, a CRC is normally used for *detection*, not
correction; there is no general single-error correction rule.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray


def _validated_bits(bits: ArrayLike, name: str = "bits") -> NDArray[np.uint8]:
    values = np.asarray(bits)
    if values.ndim != 1:
        raise ValueError(f"{name} must be a one-dimensional array")
    if not np.all((values == 0) | (values == 1)):
        raise ValueError(f"{name} must contain only 0 and 1")
    return values.astype(np.uint8, copy=False)


@dataclass(frozen=True)
class CRC:
    """A bit-oriented CRC with an explicit generator polynomial.

    Parameters
    ----------
    width
        Number of parity bits, i.e. the degree of the generator polynomial.
    polynomial
        The lower ``width`` coefficients of the generator, as an integer. The
        implicit leading coefficient ``x**width`` is not stored. For example,
        CRC-8 uses ``width=8`` and ``polynomial=0x07`` for
        ``x**8 + x**2 + x + 1``.
    """

    width: int
    polynomial: int

    def __post_init__(self) -> None:
        if isinstance(self.width, bool) or not isinstance(self.width, int):
            raise ValueError("width must be an integer")
        if not 1 <= self.width <= 64:
            raise ValueError("width must be between 1 and 64")
        if isinstance(self.polynomial, bool) or not isinstance(self.polynomial, int):
            raise ValueError("polynomial must be an integer")
        if not 0 < self.polynomial < (1 << self.width):
            raise ValueError("polynomial must fit in width bits and be positive")
        if self.polynomial & 1 == 0:
            raise ValueError("polynomial must have a nonzero constant term")

    @property
    def _poly_bits(self) -> NDArray[np.uint8]:
        return np.asarray(
            [
                (self.polynomial >> shift) & 1
                for shift in range(self.width, -1, -1)
            ],
            dtype=np.uint8,
        )

    def remainder(self, bits: ArrayLike) -> int:
        """Return the CRC remainder of a bit sequence as an integer."""
        values = _validated_bits(bits)
        if values.size == 0:
            raise ValueError("bits must be nonempty")
        work = np.concatenate(
            [values, np.zeros(self.width, dtype=np.uint8)]
        )
        polynomial = self._poly_bits
        for index in range(values.size):
            if work[index]:
                work[index : index + self.width + 1] ^= polynomial
        remainder = 0
        for bit in work[-self.width :]:
            remainder = (remainder << 1) | int(bit)
        return remainder

    def encode(self, bits: ArrayLike) -> NDArray[np.uint8]:
        """Return the systematic codeword: message bits followed by the CRC."""
        values = _validated_bits(bits)
        if values.size == 0:
            raise ValueError("bits must be nonempty")
        remainder = self.remainder(values)
        parity = np.asarray(
            [(remainder >> shift) & 1 for shift in range(self.width - 1, -1, -1)],
            dtype=np.uint8,
        )
        return np.concatenate([values, parity])

    def check(self, codeword: ArrayLike) -> bool:
        """Return whether a received word has zero CRC remainder."""
        values = _validated_bits(codeword, "codeword")
        if values.size < self.width:
            raise ValueError("codeword must contain at least width bits")
        return self.remainder(values) == 0

    def syndrome(self, codeword: ArrayLike) -> int:
        """Return the CRC remainder of a received word."""
        return self.remainder(codeword)

    def division_trace(self, bits: ArrayLike) -> NDArray[np.uint8]:
        """Return the register contents after each message-bit step.

        The result has shape ``(len(bits), width)``. It is intended for
        visualising the shift-register division that computes the remainder.
        """
        values = _validated_bits(bits)
        if values.size == 0:
            raise ValueError("bits must be nonempty")
        work = np.concatenate(
            [values, np.zeros(self.width, dtype=np.uint8)]
        )
        polynomial = self._poly_bits
        trace = np.empty((values.size, self.width), dtype=np.uint8)
        for index in range(values.size):
            if work[index]:
                work[index : index + self.width + 1] ^= polynomial
            trace[index] = work[-self.width :]
        return trace


CRC8 = CRC(width=8, polynomial=0x07)
"""CRC-8 with ``g(x) = x**8 + x**2 + x + 1``."""

CRC16 = CRC(width=16, polynomial=0x1021)
"""CRC-16/CCITT with ``g(x) = x**16 + x**12 + x**5 + 1``."""

CRC32 = CRC(width=32, polynomial=0x04C11DB7)
"""CRC-32 with the IEEE 802.3 generator polynomial."""
