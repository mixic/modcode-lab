# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Short regular (3,6) LDPC code with a normalized min-sum decoder.

The construction is quasi-cyclic: six block rows and twelve block columns of
circulant permutations with block size 84, giving N = 1008 and M = 504 checks.
Two checks are linearly dependent, so the rank is 502 and the code carries
506 information bits; the rate is 506/1008 = 0.502, approximately one half.

The decoder is normalized min-sum, a practical approximation to sum-product
belief propagation.  It is deliberately small and pedagogical, not a 5G LDPC
implementation: there is no lifting-family structure, rate matching, CRC-aided
decoding, or layered schedule.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from modcode.fec.base import FEC, BitArray


class LDPC36(FEC):
    """Quasi-cyclic (3,6) LDPC code, N=1008, approximately rate 1/2."""

    block_size = 84
    block_rows = 6
    block_columns = 12
    n = block_size * block_columns
    m = block_size * block_rows

    def __init__(
        self,
        *,
        seed: int = 0,
        alpha: float = 0.8,
        max_iterations: int = 50,
    ) -> None:
        if not 0.0 < alpha <= 1.0:
            raise ValueError("alpha must be in (0, 1]")
        if max_iterations < 1:
            raise ValueError("max_iterations must be positive")
        self.alpha = float(alpha)
        self.max_iterations = int(max_iterations)
        self.H = self._build_parity_check(seed)
        self._rref, self._pivots = self._row_reduce(self.H)
        self.rank = len(self._pivots)
        self._free = np.asarray(
            [
                column
                for column in range(self.n)
                if column not in set(self._pivots)
            ],
            dtype=np.intp,
        )
        self.k = self._free.size
        self.G = self._build_generator()
        self._check_vars = [
            np.flatnonzero(self.H[row]).astype(np.intp)
            for row in range(self.m)
        ]
        self._build_edges()

    def _build_parity_check(self, seed: int) -> NDArray[np.uint8]:
        rng = np.random.default_rng(seed)
        shifts = rng.integers(
            0, self.block_size, size=(self.block_rows, self.block_columns)
        )
        H = np.zeros((self.m, self.n), dtype=np.uint8)
        for row_block in range(self.block_rows):
            columns = [
                (2 * row_block + offset) % self.block_columns
                for offset in range(6)
            ]
            for column_block in columns:
                shift = int(shifts[row_block, column_block])
                for index in range(self.block_size):
                    H[
                        row_block * self.block_size + index,
                        column_block * self.block_size
                        + ((index + shift) % self.block_size),
                    ] = 1
        return H

    @staticmethod
    def _row_reduce(
        matrix: NDArray[np.uint8],
    ) -> tuple[NDArray[np.uint8], list[int]]:
        reduced = matrix.copy()
        rows, columns = reduced.shape
        pivots: list[int] = []
        pivot_row = 0
        for column in range(columns):
            if pivot_row >= rows:
                break
            candidates = np.flatnonzero(reduced[pivot_row:, column])
            if candidates.size == 0:
                continue
            selected = pivot_row + int(candidates[0])
            if selected != pivot_row:
                reduced[[pivot_row, selected]] = reduced[[selected, pivot_row]]
            other_rows = np.flatnonzero(reduced[:, column])
            other_rows = other_rows[other_rows != pivot_row]
            reduced[other_rows] ^= reduced[pivot_row]
            pivots.append(column)
            pivot_row += 1
        return reduced, pivots

    def _build_generator(self) -> NDArray[np.uint8]:
        generator = np.zeros((self.k, self.n), dtype=np.uint8)
        for row, free_column in enumerate(self._free):
            generator[row, free_column] = 1
            generator[row, self._pivots] = self._rref[
                : self.rank, free_column
            ]
        return generator

    def _build_edges(self) -> None:
        self._check_edges: list[NDArray[np.intp]] = []
        self._var_edges: list[list[int]] = [[] for _ in range(self.n)]
        self._edge_var: list[int] = []
        edge = 0
        for variables in self._check_vars:
            ids = []
            for variable in variables:
                self._edge_var.append(int(variable))
                self._var_edges[int(variable)].append(edge)
                ids.append(edge)
                edge += 1
            self._check_edges.append(np.asarray(ids, dtype=np.intp))
        self._edge_count = edge
        self._var_edge_arrays = [
            np.asarray(ids, dtype=np.intp) for ids in self._var_edges
        ]
        self._edge_var_array = np.asarray(self._edge_var, dtype=np.intp)

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
            raise ValueError(f"LDPC requires a multiple of {self.k} data bits")
        blocks = values.reshape(-1, self.k)
        encoded = np.empty((blocks.shape[0], self.n), dtype=np.uint8)
        for index, block in enumerate(blocks):
            encoded[index] = (block.astype(np.intp) @ self.G) % 2
        return encoded.reshape(-1)

    def decode(self, bits: NDArray[np.integer]) -> BitArray:
        values = self.validated_bits(bits)
        if values.size % self.n:
            raise ValueError(f"LDPC requires a multiple of {self.n} code bits")
        blocks = values.reshape(-1, self.n)
        decoded = np.empty((blocks.shape[0], self.k), dtype=np.uint8)
        for index, block in enumerate(blocks):
            llr = (1.0 - 2.0 * block.astype(np.float64)) * 4.0
            decoded[index] = self._decode_llr(llr)[self._free]
        return decoded.reshape(-1)

    def decode_soft(self, llrs: NDArray[np.floating]) -> BitArray:
        values = np.asarray(llrs, dtype=np.float64)
        if values.ndim != 1:
            raise ValueError("llrs must be a one-dimensional array")
        if not np.all(np.isfinite(values)):
            raise ValueError("llrs must contain only finite values")
        if values.size % self.n:
            raise ValueError(f"LDPC requires a multiple of {self.n} LLRs")
        blocks = values.reshape(-1, self.n)
        decoded = np.empty((blocks.shape[0], self.k), dtype=np.uint8)
        for index, block in enumerate(blocks):
            decoded[index] = self._decode_llr(block)[self._free]
        return decoded.reshape(-1)

    def _decode_llr(self, channel_llr: NDArray[np.float64]) -> BitArray:
        v2c = channel_llr[self._edge_var_array].copy()
        c2v = np.zeros(self._edge_count, dtype=np.float64)
        for _ in range(self.max_iterations):
            for edge_ids in self._check_edges:
                values = v2c[edge_ids]
                signs = np.where(values < 0.0, -1.0, 1.0)
                magnitudes = np.abs(values)
                order = np.argsort(magnitudes)
                minimum = magnitudes[order[0]]
                second = magnitudes[order[1]]
                product = float(np.prod(signs))
                for position, edge in enumerate(edge_ids):
                    c2v[edge] = (
                        self.alpha
                        * product
                        * signs[position]
                        * (second if position == order[0] else minimum)
                    )
            total = channel_llr + np.asarray(
                [c2v[ids].sum() for ids in self._var_edge_arrays]
            )
            hard = (total < 0.0).astype(np.uint8)
            if self._syndrome_ok(hard):
                return hard
            for variable, edge_ids in enumerate(self._var_edge_arrays):
                v2c[edge_ids] = total[variable] - c2v[edge_ids]
        return (channel_llr < 0.0).astype(np.uint8)

    def _syndrome_ok(self, bits: NDArray[np.uint8]) -> bool:
        return all(
            int(bits[variables].sum()) % 2 == 0
            for variables in self._check_vars
        )


__all__ = ["LDPC36"]
