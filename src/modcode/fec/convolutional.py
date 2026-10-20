# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Rate-1/2 convolutional code with hard and soft Viterbi decoding.

The code has constraint length 3, memory 2, and generator polynomials
``g1 = 1 + D + D^2`` (octal 7) and ``g2 = 1 + D^2`` (octal 5).  Each block
of ``block_bits`` data bits is terminated by two zero bits, so the encoded
block is ``2 * (block_bits + 2)`` bits and the rate is slightly below 1/2.

Soft decoding expects LLRs with the same convention as
:meth:`modcode.modulation.base.Modulation.soft_demodulate`: positive means
bit 0 is more likely.  The Viterbi metric is the correlation between the
expected bipolar codeword and the LLRs.
"""

from __future__ import annotations

from numbers import Integral

import numpy as np
from numpy.typing import NDArray

from modcode.fec.base import FEC, BitArray


class ConvolutionalCode(FEC):
    """Terminated rate-1/2 convolutional code with Viterbi decoding."""

    _memory = 2
    _states = 4

    def __init__(self, block_bits: int = 64) -> None:
        if isinstance(block_bits, bool) or not isinstance(block_bits, Integral):
            raise ValueError("block_bits must be an integer")
        value = int(block_bits)
        if value < 1:
            raise ValueError("block_bits must be positive")
        self.block_bits = value
        self._next_state = np.zeros((self._states, 2), dtype=np.intp)
        self._output = np.zeros((self._states, 2, 2), dtype=np.uint8)
        for state in range(self._states):
            previous_one = (state >> 1) & 1
            previous_two = state & 1
            for bit in (0, 1):
                self._next_state[state, bit] = (
                    (bit << 1) | previous_one
                ) & 0b11
                self._output[state, bit, 0] = (
                    bit ^ previous_one ^ previous_two
                )
                self._output[state, bit, 1] = bit ^ previous_two

    @property
    def data_bits_per_block(self) -> int:
        return self.block_bits

    @property
    def encoded_bits_per_block(self) -> int:
        return 2 * (self.block_bits + self._memory)

    @property
    def rate(self) -> float:
        return self.block_bits / self.encoded_bits_per_block

    def encode(self, bits: NDArray[np.integer]) -> BitArray:
        values = self.validated_bits(bits)
        if values.size % self.block_bits:
            raise ValueError(
                "convolutional encoder requires a multiple of block_bits"
            )
        blocks = values.reshape(-1, self.block_bits)
        encoded = np.empty(
            (blocks.shape[0], self.encoded_bits_per_block), dtype=np.uint8
        )
        for index, block in enumerate(blocks):
            encoded[index] = self._encode_block(block)
        return encoded.reshape(-1)

    def decode(self, bits: NDArray[np.integer]) -> BitArray:
        values = self.validated_bits(bits)
        return self._decode_blocks(values.astype(np.float64), soft=False)

    def decode_soft(self, llrs: NDArray[np.floating]) -> BitArray:
        values = np.asarray(llrs, dtype=np.float64)
        if values.ndim != 1:
            raise ValueError("llrs must be a one-dimensional array")
        if not np.all(np.isfinite(values)):
            raise ValueError("llrs must contain only finite values")
        return self._decode_blocks(values, soft=True)

    def _encode_block(self, block: NDArray[np.uint8]) -> NDArray[np.uint8]:
        padded = np.concatenate(
            [block, np.zeros(self._memory, dtype=np.uint8)]
        )
        encoded = np.empty(2 * padded.size, dtype=np.uint8)
        state = 0
        for step, bit in enumerate(padded):
            encoded[2 * step : 2 * step + 2] = self._output[state, bit]
            state = self._next_state[state, bit]
        return encoded

    def _decode_blocks(
        self, values: NDArray[np.float64], *, soft: bool
    ) -> BitArray:
        if values.size % self.encoded_bits_per_block:
            raise ValueError(
                "convolutional decoder requires whole encoded blocks"
            )
        blocks = values.reshape(-1, self.encoded_bits_per_block)
        decoded = np.empty((blocks.shape[0], self.block_bits), dtype=np.uint8)
        for index, block in enumerate(blocks):
            decoded[index] = self._viterbi(block, soft=soft)
        return decoded.reshape(-1)

    def _viterbi(
        self, received: NDArray[np.float64], *, soft: bool
    ) -> NDArray[np.uint8]:
        steps = self.block_bits + self._memory
        metrics = np.full((steps + 1, self._states), np.inf)
        previous_state = np.zeros((steps + 1, self._states), dtype=np.intp)
        previous_input = np.zeros((steps + 1, self._states), dtype=np.uint8)
        metrics[0, 0] = 0.0

        for step in range(steps):
            pair = received[2 * step : 2 * step + 2]
            for state in range(self._states):
                base = metrics[step, state]
                if not np.isfinite(base):
                    continue
                for bit in (0, 1):
                    expected = self._output[state, bit]
                    if soft:
                        cost = -float(
                            np.sum((1.0 - 2.0 * expected) * pair)
                        )
                    else:
                        cost = float(
                            np.count_nonzero(pair != expected)
                        )
                    candidate = base + cost
                    next_state = self._next_state[state, bit]
                    if candidate < metrics[step + 1, next_state]:
                        metrics[step + 1, next_state] = candidate
                        previous_state[step + 1, next_state] = state
                        previous_input[step + 1, next_state] = bit

        state = 0
        if not np.isfinite(metrics[steps, state]):
            state = int(np.argmin(metrics[steps]))
        decoded = np.empty(steps, dtype=np.uint8)
        for step in range(steps, 0, -1):
            decoded[step - 1] = previous_input[step, state]
            state = previous_state[step, state]
        return decoded[: self.block_bits]


__all__ = ["ConvolutionalCode"]
