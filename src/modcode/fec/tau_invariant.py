# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Involution-invariant block codes with RM(1,3) as the canonical example.

A code C is invariant under a coordinate involution tau if tau(C) = C.  This
is a **permutation-invariance** condition, not a convolution: cyclic codes are
invariant under a cyclic permutation, while RM(1,3) is invariant under the
non-cyclic coordinate flip i -> i XOR 1.

Cyclic codes are not a subclass of tau-invariant codes for a fixed tau; the
correct statement is that both are families of permutation-invariant codes,
indexed by different permutations.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from modcode.fec.base import FEC, BitArray


def _rm13_generator() -> NDArray[np.uint8]:
    """Return the 4-by-8 generator of RM(1,3)."""
    return np.asarray(
        [
            [1, 1, 1, 1, 1, 1, 1, 1],
            [0, 0, 0, 0, 1, 1, 1, 1],
            [0, 0, 1, 1, 0, 0, 1, 1],
            [0, 1, 0, 1, 0, 1, 0, 1],
        ],
        dtype=np.uint8,
    )


def _tau_permutation(n: int = 8) -> NDArray[np.intp]:
    """Return the fixed-point-free involution i -> i XOR 1."""
    return np.asarray([index ^ 1 for index in range(n)], dtype=np.intp)


class TauInvariantCode(FEC):
    """A binary code invariant under a coordinate involution.

    Parameters
    ----------
    generator
        ``(k, n)`` generator matrix over GF(2), with independent rows.
    tau_perm
        ``(n,)`` permutation array satisfying ``tau_perm[tau_perm] = arange(n)``.
    """

    def __init__(
        self,
        generator: NDArray[np.integer],
        tau_perm: NDArray[np.integer],
    ) -> None:
        matrix = np.asarray(generator, dtype=np.uint8)
        permutation = np.asarray(tau_perm, dtype=np.intp)
        if matrix.ndim != 2 or matrix.size == 0:
            raise ValueError("generator must be a nonempty 2-D matrix")
        if not np.all((matrix == 0) | (matrix == 1)):
            raise ValueError("generator must be binary")
        if permutation.ndim != 1 or permutation.size != matrix.shape[1]:
            raise ValueError("tau_perm length must equal the code length")
        if not np.array_equal(np.sort(permutation), np.arange(permutation.size)):
            raise ValueError("tau_perm must be a permutation")
        if not np.array_equal(permutation[permutation], np.arange(permutation.size)):
            raise ValueError("tau_perm is not an involution")

        self.generator = matrix
        self.tau_perm = permutation
        self.k, self.n = matrix.shape
        self.codebook = self._enumerate_codebook()
        self._codebook_set = {tuple(row) for row in self.codebook}
        if len(self._codebook_set) != (1 << self.k):
            raise ValueError("generator rows are not independent")
        self._verify_invariance()

    def _enumerate_codebook(self) -> NDArray[np.uint8]:
        rows = 1 << self.k
        codebook = np.empty((rows, self.n), dtype=np.uint8)
        for index in range(rows):
            coefficients = np.asarray(
                [
                    (index >> (self.k - 1 - bit)) & 1
                    for bit in range(self.k)
                ],
                dtype=np.uint8,
            )
            codebook[index] = (coefficients @ self.generator) % 2
        return codebook

    def _verify_invariance(self) -> None:
        for codeword in self.codebook:
            if tuple(codeword[self.tau_perm]) not in self._codebook_set:
                raise ValueError(
                    "code is not tau-invariant: the image of a codeword "
                    "is not in the codebook"
                )

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
            raise ValueError(f"bit length must be a multiple of k={self.k}")
        blocks = values.reshape(-1, self.k)
        return ((blocks.astype(np.intp) @ self.generator) % 2).reshape(-1)

    def decode(self, bits: NDArray[np.integer]) -> BitArray:
        values = self.validated_bits(bits)
        if values.size % self.n:
            raise ValueError(f"bit length must be a multiple of n={self.n}")
        blocks = values.reshape(-1, self.n)
        decoded = np.empty((blocks.shape[0], self.k), dtype=np.uint8)
        for index, block in enumerate(blocks):
            distances = np.sum(self.codebook != block, axis=1)
            decoded[index] = self._index_to_bits(int(np.argmin(distances)))
        return decoded.reshape(-1)

    def decode_soft(self, llrs: NDArray[np.floating]) -> BitArray:
        values = np.asarray(llrs, dtype=np.float64)
        if values.ndim != 1:
            raise ValueError("llrs must be a one-dimensional array")
        if not np.all(np.isfinite(values)):
            raise ValueError("llrs must contain only finite values")
        if values.size % self.n:
            raise ValueError(f"LLR length must be a multiple of n={self.n}")
        blocks = values.reshape(-1, self.n)
        signs = 1.0 - 2.0 * self.codebook.astype(np.float64)
        best = np.argmax(blocks @ signs.T, axis=1)
        return np.asarray(
            [self._index_to_bits(int(index)) for index in best],
            dtype=np.uint8,
        ).reshape(-1)

    def min_distance(self) -> int:
        distance = self.n
        for left in range(self.codebook.shape[0]):
            distances = np.sum(
                self.codebook[left] != self.codebook[left + 1 :], axis=1
            )
            if distances.size:
                distance = min(distance, int(np.min(distances)))
        return distance

    def tau_action(self, codeword: NDArray[np.integer]) -> NDArray[np.uint8]:
        values = np.asarray(codeword, dtype=np.uint8)
        if values.ndim != 1 or values.size != self.n:
            raise ValueError("codeword must have length n")
        return values[self.tau_perm]

    def _index_to_bits(self, index: int) -> NDArray[np.uint8]:
        return np.asarray(
            [
                (index >> (self.k - 1 - bit)) & 1
                for bit in range(self.k)
            ],
            dtype=np.uint8,
        )


def rm13() -> TauInvariantCode:
    """Return RM(1,3), invariant under ``i -> i XOR 1`` and not cyclic."""
    return TauInvariantCode(_rm13_generator(), _tau_permutation(8))


__all__ = ["TauInvariantCode", "rm13"]
