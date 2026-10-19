# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Educational Polar coding for the binary erasure channel (BEC)."""

from numbers import Integral

import numpy as np
from numpy.typing import NDArray

from modcode.fec.base import BitArray


def bec_synthetic_erasure_probabilities(
    block_length: int,
    erasure_probability: float,
) -> NDArray[np.float64]:
    """Return the BEC erasure probability of each synthesized bit-channel.

    The recursion is exact for the Arikan transform under successive
    cancellation: ``epsilon_minus = 2*epsilon - epsilon**2`` and
    ``epsilon_plus = epsilon**2``.
    """
    _validate_block_length(block_length)
    _validate_erasure_probability(erasure_probability)
    probabilities = np.asarray([erasure_probability], dtype=np.float64)
    while probabilities.size < block_length:
        epsilon = probabilities
        probabilities = np.column_stack(
            (2.0 * epsilon - epsilon**2, epsilon**2)
        ).reshape(-1)
    return probabilities


def bec_capacity(erasure_probability: float) -> float:
    """Return the Shannon capacity in bits/use of a binary erasure channel."""
    _validate_erasure_probability(erasure_probability)
    return 1.0 - erasure_probability


class PolarBEC:
    """A short Polar code with BEC-designed information positions and SC decode.

    Information positions are selected as the ``k`` synthetic channels with
    the smallest exact BEC erasure probabilities. Frozen bits are fixed to
    zero. This pedagogical implementation does not include CRC-aided list
    decoding, rate matching, or construction for general channels.
    """

    def __init__(
        self,
        block_length: int,
        information_bits: int,
        erasure_probability: float,
    ) -> None:
        _validate_block_length(block_length)
        if (
            isinstance(information_bits, bool)
            or not isinstance(information_bits, Integral)
            or not 0 <= information_bits <= block_length
        ):
            raise ValueError("information_bits must be an integer in [0, block_length]")
        _validate_erasure_probability(erasure_probability)
        self.block_length = int(block_length)
        self.information_bits = int(information_bits)
        self.erasure_probability = float(erasure_probability)
        self.synthetic_erasure_probabilities = bec_synthetic_erasure_probabilities(
            self.block_length, self.erasure_probability
        )
        ranked = np.argsort(self.synthetic_erasure_probabilities, kind="stable")
        self.information_indices = np.sort(ranked[: self.information_bits])
        self.frozen_mask = np.ones(self.block_length, dtype=bool)
        self.frozen_mask[self.information_indices] = False

    @property
    def rate(self) -> float:
        """Return the code rate."""
        return self.information_bits / self.block_length

    def encode(self, information: NDArray[np.integer]) -> BitArray:
        """Place information bits, set frozen bits to zero, and polar transform."""
        values = _validated_bits(information)
        if values.size != self.information_bits:
            raise ValueError(
                f"expected {self.information_bits} information bits, got {values.size}"
            )
        source = np.zeros(self.block_length, dtype=np.uint8)
        source[self.information_indices] = values
        return _polar_transform(source)

    def decode_bec(self, received: NDArray[np.integer]) -> BitArray:
        """Successively-cancel decode a BEC word using -1 for erased symbols."""
        values = np.asarray(received)
        if values.ndim != 1 or values.size != self.block_length:
            raise ValueError("received must be one-dimensional with block_length items")
        if not np.all((values == -1) | (values == 0) | (values == 1)):
            raise ValueError("received values must be 0, 1, or -1 for erasure")
        decoded_source = _sc_decode(
            values.astype(np.int8, copy=False),
            self.frozen_mask,
        )
        return decoded_source[self.information_indices].astype(np.uint8)


def _sc_decode(
    received: NDArray[np.int8],
    frozen_mask: NDArray[np.bool_],
) -> NDArray[np.uint8]:
    if received.size == 1:
        if frozen_mask[0]:
            return np.zeros(1, dtype=np.uint8)
        return np.asarray([0 if received[0] == -1 else received[0]], dtype=np.uint8)

    half = received.size // 2
    first = received[:half]
    second = received[half:]
    left_observations = _bec_xor(first, second)
    left_source = _sc_decode(left_observations, frozen_mask[:half])
    left_codeword = _polar_transform(left_source)

    adjusted_first = np.where(
        (first == -1) | (left_codeword == 0),
        first,
        first ^ 1,
    ).astype(np.int8)
    right_observations = _bec_repetition(adjusted_first, second)
    right_source = _sc_decode(right_observations, frozen_mask[half:])
    return np.concatenate((left_source, right_source))


def _bec_xor(
    first: NDArray[np.int8],
    second: NDArray[np.int8],
) -> NDArray[np.int8]:
    known = (first != -1) & (second != -1)
    output = np.full(first.size, -1, dtype=np.int8)
    output[known] = first[known] ^ second[known]
    return output


def _bec_repetition(
    first: NDArray[np.int8],
    second: NDArray[np.int8],
) -> NDArray[np.int8]:
    output = np.where(first == -1, second, first).astype(np.int8)
    conflicting = (first != -1) & (second != -1) & (first != second)
    output[conflicting] = -1
    return output


def _polar_transform(bits: NDArray[np.integer]) -> BitArray:
    values = np.asarray(bits, dtype=np.uint8)
    if values.size == 1:
        return values.copy()
    half = values.size // 2
    left = _polar_transform(values[:half])
    right = _polar_transform(values[half:])
    return np.concatenate((left ^ right, right))


def _validated_bits(bits: NDArray[np.integer]) -> BitArray:
    values = np.asarray(bits)
    if values.ndim != 1:
        raise ValueError("information must be a one-dimensional array")
    if not np.all((values == 0) | (values == 1)):
        raise ValueError("information bits must contain only 0 and 1")
    return values.astype(np.uint8, copy=False)


def _validate_block_length(block_length: int) -> None:
    if (
        isinstance(block_length, bool)
        or not isinstance(block_length, Integral)
        or block_length < 1
        or int(block_length) & (int(block_length) - 1)
    ):
        raise ValueError("block_length must be a positive power of two")


def _validate_erasure_probability(erasure_probability: float) -> None:
    if not np.isfinite(erasure_probability) or not 0.0 <= erasure_probability <= 1.0:
        raise ValueError("erasure_probability must be finite and in [0, 1]")
