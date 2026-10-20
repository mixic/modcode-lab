# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Small binary-extension Galois fields and polynomial helpers.

Polynomials are stored low-to-high: ``[c0, c1, ...]`` represents
``c0 + c1 x + ...``.  This module is intentionally small and scalar-oriented;
it is used by the pedagogical BCH and Reed-Solomon codes, not by a
performance-critical inner loop.
"""

from __future__ import annotations

from numbers import Integral

import numpy as np
from numpy.typing import NDArray


class GF2m:
    """The field ``GF(2^m)`` with a primitive polynomial."""

    def __init__(self, m: int, primitive: int) -> None:
        if isinstance(m, bool) or not isinstance(m, Integral) or m < 2:
            raise ValueError("m must be an integer >= 2")
        if (
            isinstance(primitive, bool)
            or not isinstance(primitive, Integral)
            or not (1 << m) <= primitive < (1 << (m + 1))
        ):
            raise ValueError("primitive must have degree m")
        self.m = int(m)
        self.primitive = int(primitive)
        self.size = 1 << self.m
        self.order = self.size - 1
        self.exp = np.zeros(2 * self.size, dtype=np.intp)
        self.log = np.full(self.size, -1, dtype=np.intp)
        value = 1
        for power in range(self.order):
            self.exp[power] = value
            self.log[value] = power
            value <<= 1
            if value & self.size:
                value ^= self.primitive
            value &= self.size - 1
        if value != 1:
            raise ValueError("primitive polynomial does not generate the field")
        for power in range(self.order, 2 * self.size):
            self.exp[power] = self.exp[power - self.order]

    def add(self, left: int, right: int) -> int:
        return int(left) ^ int(right)

    def mul(self, left: int, right: int) -> int:
        if left == 0 or right == 0:
            return 0
        return int(self.exp[self.log[left] + self.log[right]])

    def inv(self, value: int) -> int:
        if value == 0:
            raise ZeroDivisionError("zero has no multiplicative inverse")
        return int(self.exp[self.order - self.log[value]])

    def div(self, left: int, right: int) -> int:
        if right == 0:
            raise ZeroDivisionError("division by zero")
        if left == 0:
            return 0
        return self.mul(left, self.inv(right))

    def pow(self, value: int, exponent: int) -> int:
        if exponent < 0:
            return self.pow(self.inv(value), -exponent)
        if value == 0:
            return 0 if exponent else 1
        return int(self.exp[(self.log[value] * exponent) % self.order])

    def poly_eval(self, coefficients: NDArray[np.integer], x: int) -> int:
        result = 0
        for coefficient in np.asarray(coefficients, dtype=np.intp)[::-1]:
            result = self.mul(result, x) ^ int(coefficient)
        return result

    def poly_mul(
        self,
        left: NDArray[np.integer],
        right: NDArray[np.integer],
    ) -> NDArray[np.intp]:
        a = np.asarray(left, dtype=np.intp)
        b = np.asarray(right, dtype=np.intp)
        result = np.zeros(a.size + b.size - 1, dtype=np.intp)
        for i, ai in enumerate(a):
            if ai == 0:
                continue
            for j, bj in enumerate(b):
                result[i + j] ^= self.mul(int(ai), int(bj))
        return self._trim(result)

    def poly_mod(
        self,
        dividend: NDArray[np.integer],
        divisor: NDArray[np.integer],
    ) -> NDArray[np.intp]:
        remainder = self._trim(np.asarray(dividend, dtype=np.intp).copy())
        divisor = self._trim(np.asarray(divisor, dtype=np.intp))
        if divisor.size == 0 or np.all(divisor == 0):
            raise ZeroDivisionError("polynomial division by zero")
        while remainder.size >= divisor.size and not np.all(remainder == 0):
            coefficient = self.div(
                int(remainder[-1]), int(divisor[-1])
            )
            shift = remainder.size - divisor.size
            for index, value in enumerate(divisor):
                remainder[index + shift] ^= self.mul(coefficient, int(value))
            remainder = self._trim(remainder)
        return remainder

    @staticmethod
    def _trim(poly: NDArray[np.intp]) -> NDArray[np.intp]:
        end = poly.size
        while end > 1 and poly[end - 1] == 0:
            end -= 1
        return poly[:end]


__all__ = ["GF2m"]
