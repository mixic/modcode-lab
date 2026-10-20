# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Pedagogical binary BCH(15,7) code with hard and soft decoding.

The code is the primitive binary BCH code of length 15 and designed distance
5 (``t = 2``), built from roots ``alpha^1 ... alpha^4`` over GF(16).  Hard
decoding uses Berlekamp-Massey and a Chien search.  Soft decoding enumerates
the 128 codewords, which is practical at this length and makes the LLR gain
visible.
"""

from __future__ import annotations

import itertools

import numpy as np
from numpy.typing import NDArray

from modcode.fec.base import FEC, BitArray
from modcode.fec.galois import GF2m


class BCH15(FEC):
    """Binary BCH(15,7,5) code."""

    n = 15
    k = 7
    t = 2

    def __init__(self) -> None:
        self.gf = GF2m(4, 0b10011)
        self.generator = self._build_generator()
        parity_bits = self.n - self.k
        if self.generator.size != parity_bits + 1:
            raise RuntimeError("unexpected BCH generator degree")
        self._data_words = np.asarray(
            list(itertools.product((0, 1), repeat=self.k)), dtype=np.uint8
        )
        self._codewords = self.encode(
            self._data_words.reshape(-1)
        ).reshape(-1, self.n)

    def _build_generator(self) -> NDArray[np.uint8]:
        roots: set[int] = set()
        for root in range(1, 2 * self.t + 1):
            current = root
            while current not in roots:
                roots.add(current)
                current = (2 * current) % self.n
        generator = np.asarray([1], dtype=np.intp)
        for root in sorted(roots):
            factor = np.asarray(
                [self.gf.exp[root], 1], dtype=np.intp
            )
            generator = self.gf.poly_mul(generator, factor)
        if np.any((generator != 0) & (generator != 1)):
            raise RuntimeError("BCH generator is not binary")
        return generator.astype(np.uint8)

    @property
    def data_bits_per_block(self) -> int:
        return self.k

    @property
    def encoded_bits_per_block(self) -> int:
        return self.n

    @property
    def rate(self) -> float:
        return self.k / self.n

    def encode(self, bits: NDArray[np.integer]) -> BitArray:
        values = self.validated_bits(bits)
        if values.size % self.k:
            raise ValueError("BCH(15,7) requires a multiple of 7 data bits")
        words = values.reshape(-1, self.k)
        encoded = np.empty((words.shape[0], self.n), dtype=np.uint8)
        parity_bits = self.n - self.k
        generator = self.generator.astype(np.intp)
        for index, word in enumerate(words):
            message = word[::-1].astype(np.intp)
            shifted = np.concatenate(
                [np.zeros(parity_bits, dtype=np.intp), message]
            )
            remainder = self.gf.poly_mod(shifted, generator)
            codeword = shifted.copy()
            codeword[: remainder.size] ^= remainder
            encoded[index] = codeword[::-1].astype(np.uint8)
        return encoded.reshape(-1)

    def decode(self, bits: NDArray[np.integer]) -> BitArray:
        values = self.validated_bits(bits)
        if values.size % self.n:
            raise ValueError("BCH(15,7) requires a multiple of 15 code bits")
        words = values.reshape(-1, self.n)
        decoded = np.empty((words.shape[0], self.k), dtype=np.uint8)
        for index, word in enumerate(words):
            decoded[index] = self._decode_word(word)[: self.k]
        return decoded.reshape(-1)

    def decode_soft(self, llrs: NDArray[np.floating]) -> BitArray:
        values = np.asarray(llrs, dtype=np.float64)
        if values.ndim != 1:
            raise ValueError("llrs must be a one-dimensional array")
        if not np.all(np.isfinite(values)):
            raise ValueError("llrs must contain only finite values")
        if values.size % self.n:
            raise ValueError("BCH(15,7) requires a multiple of 15 LLRs")
        if values.size == 0:
            return np.empty(0, dtype=np.uint8)
        words = values.reshape(-1, self.n)
        signs = 1.0 - 2.0 * self._codewords.astype(np.float64)
        best = np.argmax(words @ signs.T, axis=1)
        return self._data_words[best].reshape(-1)

    def _decode_word(self, word: NDArray[np.uint8]) -> NDArray[np.uint8]:
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

        corrected = word.copy()
        for position in positions:
            corrected[self.n - 1 - position] ^= 1
        return corrected

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

        return self.gf._trim(connection)


__all__ = ["BCH15"]
