# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""A small tail-biting, cyclically equivariant 8-PSK trellis-coded example."""

import numpy as np
from numpy.typing import NDArray

from modcode.modulation.base import BitArray, SymbolArray


class TailBitingTCM8PSK:
    """Rate-2/3 four-state trellis code mapped to Gray-labeled 8-PSK.

    Each trellis step consumes two bits and emits a three-bit 8-PSK label.
    State is the previous two-bit input pair. The encoder uses tail-biting
    initialization, so rotating the sequence of input pairs rotates its
    encoded symbol sequence. This is a pedagogical construction; it is not a
    claim that every scheme called cyclic TCM has this property.
    """

    bits_per_symbol = 2
    output_label_bits = 3
    state_count = 4
    constellation = np.exp(2j * np.pi * np.arange(8) / 8)

    def encode(self, bits: NDArray[np.integer]) -> SymbolArray:
        """Encode an even-length bit sequence to unit-energy 8-PSK symbols."""
        values = self._validate_bits(bits)
        if values.size == 0 or values.size % 2:
            raise ValueError("tail-biting TCM requires a nonempty even bit count")
        pairs = values.reshape(-1, 2)
        previous_pairs = np.roll(pairs, 1, axis=0)
        parity = pairs[:, 0] ^ previous_pairs[:, 0] ^ previous_pairs[:, 1]
        labels = np.column_stack((pairs, parity)).astype(np.uint8)
        gray_indexes = self._bits_to_integers(labels)
        phase_indexes = self._gray_to_binary(gray_indexes)
        return self.constellation[phase_indexes]

    def decode(self, symbols: NDArray[np.complexfloating]) -> BitArray:
        """Decode by tail-biting Viterbi search over the four start states."""
        received = np.asarray(symbols, dtype=np.complex128)
        if received.ndim != 1 or received.size == 0:
            raise ValueError("symbols must be a nonempty one-dimensional array")
        if not np.all(np.isfinite(received)):
            raise ValueError("symbols must contain only finite values")

        best_metric = float("inf")
        best_pairs: NDArray[np.uint8] | None = None
        for initial_state in range(self.state_count):
            metric = np.full(self.state_count, np.inf, dtype=np.float64)
            metric[initial_state] = 0.0
            predecessor = np.full((received.size, self.state_count), -1, dtype=np.int8)
            inputs = np.full((received.size, self.state_count), -1, dtype=np.int8)
            for time_index, sample in enumerate(received):
                next_metric = np.full(self.state_count, np.inf, dtype=np.float64)
                for state in range(self.state_count):
                    if not np.isfinite(metric[state]):
                        continue
                    for input_pair in range(self.state_count):
                        symbol_index = self._transition_symbol(state, input_pair)
                        branch_metric = (
                            abs(sample - self.constellation[symbol_index]) ** 2
                        )
                        candidate = metric[state] + branch_metric
                        if candidate < next_metric[input_pair]:
                            next_metric[input_pair] = candidate
                            predecessor[time_index, input_pair] = state
                            inputs[time_index, input_pair] = input_pair
                metric = next_metric

            if metric[initial_state] < best_metric:
                best_metric = float(metric[initial_state])
                path = np.empty(received.size, dtype=np.uint8)
                state = initial_state
                for time_index in range(received.size - 1, -1, -1):
                    path[time_index] = inputs[time_index, state]
                    state = int(predecessor[time_index, state])
                best_pairs = np.column_stack(((path >> 1) & 1, path & 1)).astype(
                    np.uint8
                )

        if best_pairs is None:
            raise RuntimeError("no valid tail-biting trellis path found")
        return best_pairs.reshape(-1)

    @staticmethod
    def _validate_bits(bits: NDArray[np.integer]) -> BitArray:
        values = np.asarray(bits)
        if values.ndim != 1:
            raise ValueError("bits must be a one-dimensional array")
        if not np.all((values == 0) | (values == 1)):
            raise ValueError("bits must contain only 0 and 1")
        return values.astype(np.uint8, copy=False)

    @staticmethod
    def _bits_to_integers(labels: NDArray[np.uint8]) -> NDArray[np.intp]:
        indexes = np.zeros(labels.shape[0], dtype=np.intp)
        for bit_index in range(labels.shape[1]):
            indexes = (indexes << 1) | labels[:, bit_index]
        return indexes

    @staticmethod
    def _gray_to_binary(gray: NDArray[np.intp]) -> NDArray[np.intp]:
        binary = gray.copy()
        shift = gray >> 1
        while np.any(shift):
            binary ^= shift
            shift >>= 1
        return binary

    @staticmethod
    def _transition_symbol(state: int, input_pair: int) -> int:
        first = (input_pair >> 1) & 1
        second = input_pair & 1
        previous_first = (state >> 1) & 1
        previous_second = state & 1
        parity = first ^ previous_first ^ previous_second
        gray_label = (first << 2) | (second << 1) | parity
        return TailBitingTCM8PSK._gray_to_binary(
            np.asarray([gray_label], dtype=np.intp)
        )[0]
