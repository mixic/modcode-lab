# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Systematic cyclic Hamming(7,4) code."""

import itertools

import numpy as np
from numpy.typing import NDArray

from modcode.fec.base import FEC, BitArray


class Hamming74(FEC):
    """Cyclic Hamming(7,4), using generator polynomial x^3 + x + 1.

    Data bits occupy the first four positions. Decoding corrects any single
    bit error in each codeword.  Soft decoding uses maximum likelihood over
    the 16 codewords, which also makes the best use of unreliable bits.
    """

    _generator = np.array([1, 0, 1, 1], dtype=np.uint8)
    _parity_checks = np.array(
        [
            [1, 1, 1, 0, 1, 0, 0],
            [0, 1, 1, 1, 0, 1, 0],
            [1, 1, 0, 1, 0, 0, 1],
        ],
        dtype=np.uint8,
    )

    def __init__(self) -> None:
        self._data_words = np.asarray(
            list(itertools.product((0, 1), repeat=4)), dtype=np.uint8
        )
        self._codewords = self.encode(self._data_words.reshape(-1)).reshape(
            -1, 7
        )
        self._syndrome_to_position = self._build_syndrome_lut()

    @classmethod
    def _build_syndrome_lut(cls) -> NDArray[np.intp]:
        """Map each nonzero syndrome to a 1-based error position."""
        lut = np.zeros(8, dtype=np.intp)
        for position in range(7):
            column = cls._parity_checks[:, position]
            value = int(column[0]) * 4 + int(column[1]) * 2 + int(column[2])
            lut[value] = position + 1
        if np.any(lut[1:] == 0):
            raise RuntimeError("Hamming parity-check columns are not distinct")
        return lut

    @property
    def data_bits_per_block(self) -> int:
        return 4

    @property
    def encoded_bits_per_block(self) -> int:
        return 7

    @property
    def rate(self) -> float:
        return 4.0 / 7.0

    def encode(self, bits: NDArray[np.integer]) -> BitArray:
        values = self.validated_bits(bits)
        if values.size % 4:
            raise ValueError("Hamming(7,4) requires a multiple of 4 data bits")
        words = values.reshape(-1, 4)
        encoded = np.empty((words.shape[0], 7), dtype=np.uint8)
        encoded[:, :4] = words
        encoded[:, 4:] = self._remainders(words)
        return encoded.reshape(-1)

    def decode(self, bits: NDArray[np.integer]) -> BitArray:
        values = self.validated_bits(bits)
        if values.size % 7:
            raise ValueError("Hamming(7,4) requires a multiple of 7 code bits")
        words = values.reshape(-1, 7).copy()
        syndromes = (words @ self._parity_checks.T) % 2
        indices = syndromes @ np.asarray([4, 2, 1], dtype=np.uint8)
        positions = self._syndrome_to_position[indices]
        error_rows = np.flatnonzero(positions)
        words[error_rows, positions[error_rows] - 1] ^= 1
        return words[:, :4].reshape(-1)

    def decode_soft(
        self, llrs: NDArray[np.floating]
    ) -> BitArray:
        """Return the maximum-likelihood data word for each 7-bit LLR block."""
        values = np.asarray(llrs, dtype=np.float64)
        if values.ndim != 1:
            raise ValueError("llrs must be a one-dimensional array")
        if not np.all(np.isfinite(values)):
            raise ValueError("llrs must contain only finite values")
        if values.size % 7:
            raise ValueError("Hamming(7,4) requires a multiple of 7 LLRs")
        if values.size == 0:
            return np.empty(0, dtype=np.uint8)
        words = values.reshape(-1, 7)
        signs = 1.0 - 2.0 * self._codewords.astype(np.float64)
        metrics = words @ signs.T
        best = np.argmax(metrics, axis=1)
        return self._data_words[best].reshape(-1)

    @classmethod
    def _remainders(cls, words: NDArray[np.uint8]) -> NDArray[np.uint8]:
        work = np.pad(words.copy(), ((0, 0), (0, 3)))
        for offset in range(4):
            active = work[:, offset] == 1
            work[active, offset : offset + 4] ^= cls._generator
        return work[:, -3:]
