# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""A toy QC-MDPC McEliece-style demonstration.

QC-MDPC McEliece uses a parity-check matrix ``H = [H0 | H1]`` whose blocks are
circulant, so ``H`` is determined by two sparse first rows. The public
generator is ``G = [I | Q]`` with ``Q = (H1^{-1} H0)^T``, and ``Q`` is again
circulant. The public key therefore has ``n`` bits instead of the ``n**2``
bits of a dense generator.

This module is deliberately tiny and is **not cryptography**. The parameters
are far below any secure size, the decoder is a toy bit-flipping/search
decoder, and no IND-CCA, side-channel, or fault-attack protections are
modelled. It exists to make the circulant-block structure and the key-size
argument visible and testable.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from math import ceil, log2

import numpy as np
from numpy.typing import ArrayLike, NDArray


def circulant(first_row: ArrayLike) -> NDArray[np.uint8]:
    """Return the circulant matrix whose first row is ``first_row``."""
    row = np.asarray(first_row, dtype=np.uint8)
    if row.ndim != 1 or row.size < 1:
        raise ValueError("first_row must be a nonempty one-dimensional array")
    if not np.all((row == 0) | (row == 1)):
        raise ValueError("first_row must contain only 0 and 1")
    return np.asarray(
        [np.roll(row, shift) for shift in range(row.size)],
        dtype=np.uint8,
    )


def is_circulant(matrix: ArrayLike) -> bool:
    """Return whether a square binary matrix is circulant."""
    values = np.asarray(matrix, dtype=np.uint8)
    if values.ndim != 2 or values.shape[0] != values.shape[1]:
        return False
    if not np.all((values == 0) | (values == 1)):
        return False
    first_row = values[0]
    return all(
        np.array_equal(values[index], np.roll(first_row, index))
        for index in range(values.shape[0])
    )


def _validated_bits(bits: ArrayLike, name: str = "bits") -> NDArray[np.uint8]:
    values = np.asarray(bits)
    if values.ndim != 1:
        raise ValueError(f"{name} must be a one-dimensional array")
    if not np.all((values == 0) | (values == 1)):
        raise ValueError(f"{name} must contain only 0 and 1")
    return values.astype(np.uint8, copy=False)


def _invert_binary(matrix: ArrayLike) -> NDArray[np.uint8] | None:
    values = np.asarray(matrix, dtype=np.uint8) % 2
    if values.ndim != 2 or values.shape[0] != values.shape[1]:
        raise ValueError("matrix must be square")
    size = values.shape[0]
    augmented = np.hstack(
        [values, np.eye(size, dtype=np.uint8)]
    )
    row = 0
    for column in range(size):
        pivot = next(
            (
                candidate
                for candidate in range(row, size)
                if augmented[candidate, column]
            ),
            None,
        )
        if pivot is None:
            return None
        if pivot != row:
            augmented[[row, pivot]] = augmented[[pivot, row]]
        for other in range(size):
            if other != row and augmented[other, column]:
                augmented[other] ^= augmented[row]
        row += 1
    return augmented[:, size:]


def _random_circulant(
    rng: np.random.Generator, block_length: int, weight: int
) -> NDArray[np.uint8]:
    first_row = np.zeros(block_length, dtype=np.uint8)
    positions = rng.choice(block_length, size=weight, replace=False)
    first_row[positions] = 1
    return circulant(first_row)


@dataclass(frozen=True)
class DecodeResult:
    """Result of the toy QC-MDPC decoder."""

    message: NDArray[np.uint8]
    corrected: NDArray[np.uint8]
    iterations: int
    used_search: bool
    syndrome_history: tuple[int, ...]


class QCMDPC:
    """A small quasi-cyclic moderate-density parity-check code.

    Parameters
    ----------
    block_length
        Size of each circulant block. The code length is ``2 * block_length``
        and the dimension is ``block_length``.
    row_weight
        Number of ones in the first row of each circulant block. Each block
        has column weight ``row_weight``; the full parity check has column
        weight ``2 * row_weight``.
    seed
        Seed for deterministic key generation.
    max_attempts
        Number of circulant pairs to try before giving up on an invertible
        ``H1``.
    """

    def __init__(
        self,
        block_length: int = 47,
        row_weight: int = 3,
        seed: int = 0,
        max_attempts: int = 1000,
    ) -> None:
        if isinstance(block_length, bool) or not isinstance(block_length, int):
            raise ValueError("block_length must be an integer")
        if block_length < 3:
            raise ValueError("block_length must be at least 3")
        if isinstance(row_weight, bool) or not isinstance(row_weight, int):
            raise ValueError("row_weight must be an integer")
        if not 1 <= row_weight < block_length:
            raise ValueError("row_weight must be in [1, block_length)")
        if row_weight % 2 == 0:
            raise ValueError(
                "row_weight must be odd for an invertible circulant H1"
            )
        if isinstance(seed, bool) or not isinstance(seed, int):
            raise ValueError("seed must be an integer")
        if max_attempts < 1:
            raise ValueError("max_attempts must be positive")

        self.block_length = block_length
        self.row_weight = row_weight
        self.seed = seed
        self.code_length = 2 * block_length
        self.dimension = block_length
        self.rate = 0.5

        rng = np.random.default_rng(seed)
        inverse: NDArray[np.uint8] | None = None
        h0 = h1 = np.empty((0, 0), dtype=np.uint8)
        for _ in range(max_attempts):
            h0 = _random_circulant(rng, block_length, row_weight)
            h1 = _random_circulant(rng, block_length, row_weight)
            inverse = _invert_binary(h1)
            if inverse is not None:
                break
        if inverse is None:
            raise RuntimeError("could not generate an invertible H1")

        self.h0 = h0
        self.h1 = h1
        self.h1_inverse = inverse
        self.parity_check = np.hstack([h0, h1])
        # G = [I | Q], with Q = (H1^{-1} H0)^T, so H G^T = 0 over GF(2).
        self.q_block = (inverse @ h0).T % 2
        self.generator = np.hstack(
            [np.eye(self.dimension, dtype=np.uint8), self.q_block]
        )

    @property
    def key_size_bits(self) -> dict[str, int]:
        """Return dense and quasi-cyclic key sizes in bits."""
        return {
            "dense_generator": self.dimension**2,
            "qc_generator": self.block_length,
            "sparse_parity": 2
            * self.row_weight
            * ceil(log2(self.block_length)),
        }

    @property
    def compression_ratio(self) -> float:
        """Dense generator bits divided by quasi-cyclic generator bits."""
        return self.dimension**2 / self.block_length

    def syndrome(self, bits: ArrayLike) -> NDArray[np.uint8]:
        """Return the binary syndrome ``H * bits``."""
        values = _validated_bits(bits)
        if values.size != self.code_length:
            raise ValueError("bits must have length 2 * block_length")
        return (
            self.parity_check.astype(np.int64) @ values.astype(np.int64)
        ).astype(np.uint8) % 2

    def encode(self, message: ArrayLike) -> NDArray[np.uint8]:
        """Encode ``block_length`` message bits to a codeword."""
        values = _validated_bits(message, "message")
        if values.size != self.dimension:
            raise ValueError("message must have block_length bits")
        return (
            values.astype(np.int64) @ self.generator.astype(np.int64)
        ).astype(np.uint8) % 2

    def encrypt(
        self,
        message: ArrayLike,
        error_weight: int = 2,
        *,
        seed: int | None = None,
    ) -> tuple[NDArray[np.uint8], NDArray[np.uint8]]:
        """Return ``(codeword + error, error)`` for a toy McEliece ciphertext."""
        if isinstance(error_weight, bool) or not isinstance(error_weight, int):
            raise ValueError("error_weight must be an integer")
        if not 0 <= error_weight <= self.code_length:
            raise ValueError("error_weight must be in [0, code_length]")
        codeword = self.encode(message)
        error = np.zeros(self.code_length, dtype=np.uint8)
        if error_weight:
            rng = np.random.default_rng(seed)
            positions = rng.choice(
                self.code_length, size=error_weight, replace=False
            )
            error[positions] = 1
        return (codeword ^ error).astype(np.uint8), error

    def bit_flip_decode(
        self, ciphertext: ArrayLike, max_iterations: int = 100
    ) -> tuple[NDArray[np.uint8], int, tuple[int, ...], bool]:
        """Run a simple bit-flipping decoder.

        Returns ``(corrected, iterations, syndrome_history, success)``.
        """
        values = _validated_bits(ciphertext, "ciphertext")
        if values.size != self.code_length:
            raise ValueError("ciphertext must have 2 * block_length bits")
        corrected = values.copy()
        syndrome = self.syndrome(corrected)
        history = [int(syndrome.sum())]
        iterations = 0
        for _ in range(max_iterations):
            if not syndrome.any():
                return corrected, iterations, tuple(history), True
            counts = (
                self.parity_check.T.astype(np.int64)
                @ syndrome.astype(np.int64)
            )
            largest = int(counts.max())
            if largest == 0:
                break
            position = int(np.flatnonzero(counts == largest)[0])
            corrected[position] ^= 1
            syndrome ^= self.parity_check[:, position]
            iterations += 1
            history.append(int(syndrome.sum()))
        return corrected, iterations, tuple(history), not syndrome.any()

    def decode(
        self,
        ciphertext: ArrayLike,
        error_weight: int = 2,
        max_iterations: int = 100,
    ) -> DecodeResult:
        """Decode a ciphertext, with a small exhaustive-search fallback."""
        corrected, iterations, history, success = self.bit_flip_decode(
            ciphertext, max_iterations
        )
        if success:
            return DecodeResult(
                message=corrected[: self.dimension].copy(),
                corrected=corrected,
                iterations=iterations,
                used_search=False,
                syndrome_history=history,
            )

        values = _validated_bits(ciphertext, "ciphertext")
        target = self.syndrome(values)
        if error_weight > 3:
            raise ValueError(
                "toy fallback supports at most three error positions"
            )
        for positions in combinations(range(self.code_length), error_weight):
            candidate = np.bitwise_xor.reduce(
                self.parity_check[:, positions], axis=1
            )
            if np.array_equal(candidate, target):
                repaired = values.copy()
                repaired[list(positions)] ^= 1
                return DecodeResult(
                    message=repaired[: self.dimension].copy(),
                    corrected=repaired,
                    iterations=iterations,
                    used_search=True,
                    syndrome_history=history,
                )
        raise ValueError("toy decoder failed to find an error pattern")

    def decrypt(
        self,
        ciphertext: ArrayLike,
        error_weight: int = 2,
        max_iterations: int = 100,
    ) -> NDArray[np.uint8]:
        """Return the recovered message bits."""
        return self.decode(ciphertext, error_weight, max_iterations).message
