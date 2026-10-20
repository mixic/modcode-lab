# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""A four-state rate-2/3 8-PSK trellis-coded modulation example.

Two information bits enter each trellis step.  The encoder appends one parity
bit and maps the three-bit label to Gray-labeled 8-PSK.  The code is a small
pedagogical TCM, not a standardized Ungerboeck code; its purpose is to show a
real trellis-coded modulation gain against uncoded QPSK at the same
2 bit/symbol spectral efficiency.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

from modcode.modulation.base import BitArray, SymbolArray


class TCM8PSK:
    """Four-state, rate-2/3, Gray-labeled 8-PSK TCM."""

    bits_per_symbol = 2
    output_label_bits = 3
    state_count = 4
    constellation = np.exp(2j * np.pi * np.arange(8) / 8)

    def encode(self, bits: NDArray[np.integer]) -> SymbolArray:
        values = self._validate_bits(bits)
        if values.size == 0 or values.size % 2:
            raise ValueError("TCM requires a nonempty even bit count")
        pairs = values.reshape(-1, 2)
        symbols = np.empty(pairs.shape[0], dtype=np.complex128)
        state = 0
        for index, (first, second) in enumerate(pairs):
            symbols[index] = self.constellation[
                self._transition_symbol(state, int(first), int(second))
            ]
            state = ((int(first) << 1) | int(second)) & 0b11
        return symbols

    def decode(self, symbols: ArrayLike) -> BitArray:
        """Soft Viterbi decode by Euclidean distance."""
        received = np.asarray(symbols, dtype=np.complex128)
        if received.ndim != 1 or received.size == 0:
            raise ValueError("symbols must be a nonempty one-dimensional array")
        if not np.all(np.isfinite(received)):
            raise ValueError("symbols must contain only finite values")

        transitions = [
            (
                state,
                first,
                second,
                ((first << 1) | second) & 0b11,
                self.constellation[
                    self._transition_symbol(state, first, second)
                ],
            )
            for state in range(self.state_count)
            for first in (0, 1)
            for second in (0, 1)
        ]

        metric = np.full(self.state_count, np.inf)
        metric[0] = 0.0
        predecessor = np.full(
            (received.size, self.state_count), -1, dtype=np.int8
        )
        inputs = np.full(
            (received.size, self.state_count), -1, dtype=np.int8
        )
        for time_index, sample in enumerate(received):
            next_metric = np.full(self.state_count, np.inf)
            for state, first, second, next_state, expected in transitions:
                if not np.isfinite(metric[state]):
                    continue
                candidate = metric[state] + abs(sample - expected) ** 2
                if candidate < next_metric[next_state]:
                    next_metric[next_state] = candidate
                    predecessor[time_index, next_state] = state
                    inputs[time_index, next_state] = (first << 1) | second
            metric = next_metric

        state = int(np.argmin(metric))
        path = np.empty(received.size, dtype=np.uint8)
        for time_index in range(received.size - 1, -1, -1):
            path[time_index] = inputs[time_index, state]
            state = int(predecessor[time_index, state])
        return np.column_stack(((path >> 1) & 1, path & 1)).reshape(-1)

    @staticmethod
    def _validate_bits(bits: NDArray[np.integer]) -> BitArray:
        values = np.asarray(bits)
        if values.ndim != 1:
            raise ValueError("bits must be a one-dimensional array")
        if not np.all((values == 0) | (values == 1)):
            raise ValueError("bits must contain only 0 and 1")
        return values.astype(np.uint8, copy=False)

    @staticmethod
    def _transition_symbol(state: int, first: int, second: int) -> int:
        previous_first = (state >> 1) & 1
        previous_second = state & 1
        parity = first ^ second ^ previous_first ^ previous_second
        gray_label = (first << 2) | (second << 1) | parity
        return TCM8PSK._gray_to_binary(gray_label)

    @staticmethod
    def _gray_to_binary(gray: int) -> int:
        value = gray
        shift = gray >> 1
        while shift:
            value ^= shift
            shift >>= 1
        return value


__all__ = ["TCM8PSK"]
