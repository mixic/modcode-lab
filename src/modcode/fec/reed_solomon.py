# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Pedagogical Reed-Solomon code over GF(2^m).

The default is RS(15,11) over GF(16), correcting two symbol errors.  Encoding
is systematic; hard decoding uses Berlekamp-Massey, Chien search, and a small
GF linear solve for the error values.  This is a lab code, not a
standards-compliant RS implementation.
"""

from __future__ import annotations

from numbers import Integral

import numpy as np
from numpy.typing import NDArray

from modcode.fec.base import FEC, BitArray
from modcode.fec.galois import GF2m


class ReedSolomon(FEC):
    """Systematic RS(2^m - 1, 2^m - 1 - 2t) code over GF(2^m)."""

    def __init__(
        self,
        m: int = 4,
        t: int = 2,
        primitive: int = 0b10011,
    ) -> None:
        if isinstance(m, bool) or not isinstance(m, Integral) or m < 2:
            raise ValueError("m must be an integer >= 2")
        if isinstance(t, bool) or not isinstance(t, Integral) or t < 1:
            raise ValueError("t must be a positive integer")
        self.m = int(m)
        self.t = int(t)
        self.n = (1 << self.m) - 1
        self.k = self.n - 2 * self.t
        if self.k < 1:
            raise ValueError("t is too large for the field size")
        self.gf = GF2m(self.m, primitive)
        self.generator = self._build_generator()

    def _build_generator(self) -> NDArray[np.intp]:
        generator = np.asarray([1], dtype=np.intp)
        for power in range(1, 2 * self.t + 1):
            factor = np.asarray(
                [self.gf.exp[power], 1], dtype=np.intp
            )
            generator = self.gf.poly_mul(generator, factor)
        return generator

    @property
    def data_bits_per_block(self) -> int:
        return self.k * self.m

    @property
    def encoded_bits_per_block(self) -> int:
        return self.n * self.m

    @property
    def rate(self) -> float:
        return self.k / self.n

    def encode(self, bits: NDArray[np.integer]) -> BitArray:
        values = self.validated_bits(bits)
        block_bits = self.k * self.m
        if values.size % block_bits:
            raise ValueError(
                f"RS({self.n},{self.k}) requires a multiple of "
                f"{block_bits} data bits"
            )
        symbols = self._bits_to_symbols(values).reshape(-1, self.k)
        encoded = np.empty((symbols.shape[0], self.n), dtype=np.uint8)
        parity = 2 * self.t
        for index, message in enumerate(symbols):
            polynomial = message[::-1].astype(np.intp)
            shifted = np.concatenate(
                [np.zeros(parity, dtype=np.intp), polynomial]
            )
            remainder = self.gf.poly_mod(shifted, self.generator)
            codeword = shifted.copy()
            codeword[: remainder.size] ^= remainder
            encoded[index] = codeword[::-1].astype(np.uint8)
        return self._symbols_to_bits(encoded.reshape(-1))

    def decode(self, bits: NDArray[np.integer]) -> BitArray:
        values = self.validated_bits(bits)
        block_bits = self.n * self.m
        if values.size % block_bits:
            raise ValueError(
                f"RS({self.n},{self.k}) requires a multiple of "
                f"{block_bits} code bits"
            )
        words = self._bits_to_symbols(values).reshape(-1, self.n)
        decoded = np.empty((words.shape[0], self.k), dtype=np.uint8)
        for index, word in enumerate(words):
            decoded[index] = self._decode_symbols(word)[: self.k]
        return self._symbols_to_bits(decoded.reshape(-1))

    def _decode_symbols(self, word: NDArray[np.uint8]) -> NDArray[np.uint8]:
        received = word.astype(np.intp)
        polynomial = received[::-1]
        syndromes = np.asarray(
            [
                self.gf.poly_eval(polynomial, self.gf.exp[power])
                for power in range(1, 2 * self.t + 1)
            ],
            dtype=np.intp,
        )
        if np.all(syndromes == 0):
            return word.copy()

        locator = self._berlekamp_massey(syndromes)
        degree = locator.size - 1
        if degree == 0 or degree > self.t:
            return word.copy()

        positions = []
        for position in range(self.n):
            root = int(self.gf.exp[(self.n - position) % self.n])
            if self.gf.poly_eval(locator, root) == 0:
                positions.append(position)
        if len(positions) != degree:
            return word.copy()

        locations = np.asarray(
            [self.gf.exp[position] for position in positions], dtype=np.intp
        )
        matrix = np.empty((2 * self.t, degree), dtype=np.intp)
        for row in range(2 * self.t):
            for column, location in enumerate(locations):
                matrix[row, column] = self.gf.pow(
                    int(location), row + 1
                )
        values = self._solve(matrix, syndromes)
        if values is None:
            return word.copy()

        corrected = received.copy()
        for position, value in zip(positions, values, strict=True):
            corrected[self.n - 1 - position] ^= int(value)
        return corrected.astype(np.uint8)

    def _berlekamp_massey(
        self, syndromes: NDArray[np.intp]
    ) -> NDArray[np.intp]:
        connection = np.asarray([1], dtype=np.intp)
        previous = np.asarray([1], dtype=np.intp)
        degree = 0
        shift = 1
        previous_discrepancy = 1

        for step in range(syndromes.size):
            discrepancy = int(syndromes[step])
            for index in range(1, degree + 1):
                if index < connection.size:
                    discrepancy ^= self.gf.mul(
                        int(connection[index]), int(syndromes[step - index])
                    )
            if discrepancy == 0:
                shift += 1
                continue

            saved = connection.copy()
            coefficient = self.gf.div(discrepancy, previous_discrepancy)
            needed = previous.size + shift
            if connection.size < needed:
                connection = np.pad(
                    connection, (0, needed - connection.size)
                )
            for index, value in enumerate(previous):
                connection[index + shift] ^= self.gf.mul(
                    coefficient, int(value)
                )

            if 2 * degree <= step:
                degree = step + 1 - degree
                previous = saved
                previous_discrepancy = discrepancy
                shift = 1
            else:
                shift += 1

        return connection[: degree + 1]

    def _solve(
        self,
        matrix: NDArray[np.intp],
        target: NDArray[np.intp],
    ) -> NDArray[np.intp] | None:
        rows, columns = matrix.shape
        augmented = np.concatenate(
            [matrix.copy(), target.reshape(-1, 1)], axis=1
        )
        for column in range(columns):
            pivot = next(
                (
                    row
                    for row in range(column, rows)
                    if augmented[row, column] != 0
                ),
                None,
            )
            if pivot is None:
                return None
            if pivot != column:
                augmented[[column, pivot]] = augmented[[pivot, column]]
            inverse = self.gf.inv(int(augmented[column, column]))
            for index in range(column, columns + 1):
                augmented[column, index] = self.gf.mul(
                    int(augmented[column, index]), inverse
                )
            for row in range(rows):
                if row == column:
                    continue
                factor = int(augmented[row, column])
                if factor == 0:
                    continue
                for index in range(column, columns + 1):
                    augmented[row, index] ^= self.gf.mul(
                        factor, int(augmented[column, index])
                    )
        return augmented[:columns, -1].astype(np.intp)

    def _bits_to_symbols(self, bits: NDArray[np.uint8]) -> NDArray[np.uint8]:
        if bits.size % self.m:
            raise ValueError("bit count must be a multiple of m")
        groups = bits.reshape(-1, self.m).astype(np.intp)
        weights = (1 << np.arange(self.m - 1, -1, -1)).astype(np.intp)
        return (groups @ weights).astype(np.uint8)

    def _symbols_to_bits(self, symbols: NDArray[np.uint8]) -> BitArray:
        shifts = np.arange(self.m - 1, -1, -1)
        return ((symbols[:, None].astype(np.intp) >> shifts) & 1).astype(
            np.uint8
        ).reshape(-1)


__all__ = ["ReedSolomon"]
