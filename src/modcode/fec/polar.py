# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Polar coding: exact BEC polarization and a short BI-AWGN SC/SCL decoder.

``PolarBEC`` is the appendix-level erasure-channel construction: it uses exact
BEC synthesis and successive cancellation on erased symbols.  ``PolarAWGN`` is
the core, soft-decision variant: Gaussian-approximation construction, SC and
successive-cancellation list (SCL) decoding on LLRs, and an optional CRC for
CRC-aided SCL.
"""

from numbers import Integral

import numpy as np
from numpy.typing import NDArray

from modcode.fec.base import FEC, BitArray
from modcode.fec.crc import CRC, CRC8


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


# ---------------------------------------------------------------------------
# Gaussian-approximation construction for the BI-AWGN channel
# ---------------------------------------------------------------------------

_QUADRATURE_NODES = np.linspace(-8.0, 8.0, 401)
_QUADRATURE_WEIGHTS = np.exp(-0.5 * _QUADRATURE_NODES**2)
_QUADRATURE_WEIGHTS /= _QUADRATURE_WEIGHTS.sum()


def gaussian_phi(x: NDArray[np.floating] | float) -> NDArray[np.float64]:
    """Return the Gaussian-approximation function ``phi(x)``.

    ``phi(x) = 1 - E[tanh(U / 2)]`` with ``U ~ N(x, 2x)`` is the mean update
    for a check-node output whose variable-node LLRs have mean ``x``.  The
    expectation is evaluated on a 401-point uniform grid of the normalised
    variable with Gaussian weights; the result is clipped to ``[1e-300, 1]``.
    """
    values = np.asarray(x, dtype=np.float64)
    positive = np.maximum(values, 1e-12)
    spread = np.sqrt(2.0 * positive)
    samples = positive[..., None] + spread[..., None] * _QUADRATURE_NODES
    means = np.tanh(0.5 * samples) @ _QUADRATURE_WEIGHTS
    return np.clip(1.0 - means, 1e-300, 1.0)


_PHI_ARGUMENTS = np.logspace(-6.0, 4.0, 2000)
_PHI_VALUES = gaussian_phi(_PHI_ARGUMENTS)
_PHI_LOG_VALUES = np.log(_PHI_VALUES)[::-1]
_PHI_LOG_ARGUMENTS = np.log(_PHI_ARGUMENTS)[::-1]


def gaussian_phi_inverse(
    probability: NDArray[np.floating] | float,
) -> NDArray[np.float64]:
    """Invert :func:`gaussian_phi` by log-log interpolation of a fixed table."""
    values = np.asarray(probability, dtype=np.float64)
    logs = np.log(np.clip(values, 1e-300, 1.0))
    return np.exp(np.interp(logs, _PHI_LOG_VALUES, _PHI_LOG_ARGUMENTS))


def awgn_synthetic_llr_means(
    block_length: int,
    design_ebn0_db: float,
    rate: float,
) -> NDArray[np.float64]:
    """Return the Gaussian-approximation LLR mean of each synthesized channel.

    The recursion is the standard Gaussian approximation: the left (check-like)
    child has mean ``phi**-1(1 - (1 - phi(m))**2)`` and the right
    (variable-like) child has mean ``2m``.  The initial mean is
    ``2 / sigma**2`` for BPSK with unit energy, where
    ``sigma**2 = 1 / (2 * rate * Eb/N0)``.
    """
    _validate_block_length(block_length)
    if not np.isfinite(design_ebn0_db):
        raise ValueError("design_ebn0_db must be finite")
    if not 0.0 < rate <= 1.0:
        raise ValueError("rate must be in (0, 1]")
    sigma_squared = 1.0 / (2.0 * rate * 10.0 ** (design_ebn0_db / 10.0))
    means = np.asarray([2.0 / sigma_squared], dtype=np.float64)
    while means.size < block_length:
        checked = gaussian_phi_inverse(1.0 - (1.0 - gaussian_phi(means)) ** 2)
        means = np.column_stack((checked, 2.0 * means)).reshape(-1)
    return means


# ---------------------------------------------------------------------------
# Successive-cancellation list decoder for the BI-AWGN channel
# ---------------------------------------------------------------------------


class _Node:
    """One node on the depth-first decoding stack of a polar factor tree."""

    __slots__ = ("alpha", "size", "next_child", "left_beta")

    def __init__(self, alpha: NDArray[np.float64], size: int) -> None:
        self.alpha = alpha
        self.size = int(size)
        self.next_child = 0
        self.left_beta: NDArray[np.uint8] | None = None


def _leaf_cost(llr: NDArray[np.float64], bit: int) -> NDArray[np.float64]:
    """Return the path-metric increment ``ln(1 + exp(-(1-2b)L))``."""
    return np.logaddexp(0.0, -(1.0 - 2.0 * bit) * llr)


class PolarAWGN(FEC):
    """A short Polar code for the BI-AWGN channel with SC/SCL decoding.

    Information positions are the ``data_bits + crc.width`` positions with the
    largest Gaussian-approximation LLR mean at ``design_ebn0_db``.  Frozen
    positions are fixed to zero.  ``list_size=1`` is successive cancellation;
    ``list_size > 1`` is successive-cancellation list decoding, and with a CRC
    the decoder keeps the minimum-metric path whose CRC checks.

    This is a pedagogical block code: no rate matching, no 5G NR sequence, no
    distributed CRC, and no claim of 3GPP conformance.
    """

    def __init__(
        self,
        block_length: int = 1024,
        data_bits: int = 512,
        *,
        crc: CRC | None = CRC8,
        list_size: int = 4,
        design_ebn0_db: float = 2.0,
    ) -> None:
        _validate_block_length(block_length)
        if isinstance(data_bits, bool) or not isinstance(data_bits, Integral):
            raise ValueError("data_bits must be an integer")
        if isinstance(list_size, bool) or not isinstance(list_size, Integral):
            raise ValueError("list_size must be an integer")
        if not np.isfinite(design_ebn0_db):
            raise ValueError("design_ebn0_db must be finite")
        self.block_length = int(block_length)
        self._data_bits = int(data_bits)
        self.crc = crc
        self.crc_width = 0 if crc is None else crc.width
        self.information_bits = self._data_bits + self.crc_width
        if self._data_bits < 1:
            raise ValueError("data_bits must be positive")
        if self.information_bits > self.block_length:
            raise ValueError("data_bits plus the CRC must not exceed block_length")
        if list_size < 1:
            raise ValueError("list_size must be positive")
        self.list_size = int(list_size)
        self.design_ebn0_db = float(design_ebn0_db)
        self.channel_means = awgn_synthetic_llr_means(
            self.block_length, self.design_ebn0_db, self.rate
        )
        ranked = np.argsort(self.channel_means, kind="stable")
        self.information_indices = np.sort(ranked[-self.information_bits :])
        self.frozen_mask = np.ones(self.block_length, dtype=bool)
        self.frozen_mask[self.information_indices] = False
        self.last_work = 0

    @property
    def data_bits_per_block(self) -> int:
        return self._data_bits

    @property
    def encoded_bits_per_block(self) -> int:
        return self.block_length

    @property
    def rate(self) -> float:
        return self._data_bits / self.block_length

    def encode(self, bits: NDArray[np.integer]) -> BitArray:
        values = self.validated_bits(bits)
        if values.size % self._data_bits:
            raise ValueError(
                f"Polar requires a multiple of {self._data_bits} data bits"
            )
        blocks = values.reshape(-1, self._data_bits)
        encoded = np.empty((blocks.shape[0], self.block_length), dtype=np.uint8)
        for index, block in enumerate(blocks):
            information = block if self.crc is None else self.crc.encode(block)
            source = np.zeros(self.block_length, dtype=np.uint8)
            source[self.information_indices] = information
            encoded[index] = _polar_transform(source)
        return encoded.reshape(-1)

    def decode(self, bits: NDArray[np.integer]) -> BitArray:
        values = self.validated_bits(bits)
        if values.size % self.block_length:
            raise ValueError(
                f"Polar requires a multiple of {self.block_length} code bits"
            )
        blocks = values.reshape(-1, self.block_length)
        decoded = np.empty(
            (blocks.shape[0], self._data_bits), dtype=np.uint8
        )
        for index, block in enumerate(blocks):
            llr = (1.0 - 2.0 * block.astype(np.float64)) * 4.0
            decoded[index] = self._decode_block(llr)
        return decoded.reshape(-1)

    def decode_soft(self, llrs: NDArray[np.floating]) -> BitArray:
        values = np.asarray(llrs)
        if np.iscomplexobj(values):
            raise ValueError("llrs must be real, not complex")
        values = values.astype(np.float64, copy=False)
        if values.ndim != 1:
            raise ValueError("llrs must be a one-dimensional array")
        if not np.all(np.isfinite(values)):
            raise ValueError("llrs must contain only finite values")
        if values.size % self.block_length:
            raise ValueError(
                f"Polar requires a multiple of {self.block_length} LLRs"
            )
        blocks = values.reshape(-1, self.block_length)
        decoded = np.empty(
            (blocks.shape[0], self._data_bits), dtype=np.uint8
        )
        for index, block in enumerate(blocks):
            decoded[index] = self._decode_block(block)
        return decoded.reshape(-1)

    def _decode_block(self, llr: NDArray[np.float64]) -> BitArray:
        self.last_work = 0
        decisions = np.zeros((1, self.block_length), dtype=np.uint8)
        metric = np.zeros(1, dtype=np.float64)
        stack = [_Node(llr.reshape(1, -1).copy(), self.block_length)]
        paths = 1
        for index in range(self.block_length):
            leaf = self._descend(stack)
            if self.frozen_mask[index]:
                metric = metric + _leaf_cost(leaf, 0)
                self._decide(stack, np.zeros(paths, dtype=np.uint8))
                continue
            candidates = np.concatenate(
                (metric + _leaf_cost(leaf, 0), metric + _leaf_cost(leaf, 1))
            )
            keep = np.argsort(candidates, kind="stable")[: self.list_size]
            parents = keep % paths
            bits = (keep // paths).astype(np.uint8)
            metric = candidates[keep]
            paths = keep.size
            decisions = decisions[parents]
            decisions[np.arange(paths), index] = bits
            for node in stack:
                node.alpha = node.alpha[parents]
                if node.left_beta is not None:
                    node.left_beta = node.left_beta[parents]
            self._decide(stack, bits)
        return self._select(decisions, metric)

    def _descend(self, stack: list[_Node]) -> NDArray[np.float64]:
        while stack[-1].size > 1:
            node = stack[-1]
            half = node.size // 2
            first = node.alpha[:, :half]
            second = node.alpha[:, half:]
            if node.next_child == 0:
                alpha = self._check_llr(first, second)
                node.next_child = 1
            else:
                alpha = self._variable_llr(first, second, node.left_beta)
                node.next_child = 2
            stack.append(_Node(alpha, half))
        return stack[-1].alpha[:, 0]

    def _decide(self, stack: list[_Node], bits: NDArray[np.uint8]) -> None:
        beta = bits.reshape(-1, 1).astype(np.uint8)
        stack.pop()
        while stack:
            node = stack[-1]
            half = node.size // 2
            if node.next_child == 1:
                node.left_beta = beta
                alpha = self._variable_llr(
                    node.alpha[:, :half], node.alpha[:, half:], beta
                )
                node.next_child = 2
                stack.append(_Node(alpha, half))
                return
            beta = np.concatenate((node.left_beta ^ beta, beta), axis=1)
            stack.pop()

    def _check_llr(
        self,
        first: NDArray[np.float64],
        second: NDArray[np.float64],
    ) -> NDArray[np.float64]:
        """Return the exact sum-product check-node LLR.

        ``2 atanh(tanh(a/2) tanh(b/2))`` is evaluated through the equivalent
        Jacobian-logarithm form, which is numerically stable and invariant to a
        positive rescaling of the input LLRs.
        """
        self.last_work += int(first.size)
        magnitude = np.minimum(np.abs(first), np.abs(second))
        correction = np.log1p(np.exp(-np.abs(first + second))) - np.log1p(
            np.exp(-np.abs(first - second))
        )
        return np.sign(first) * np.sign(second) * magnitude + correction

    def _variable_llr(
        self,
        first: NDArray[np.float64],
        second: NDArray[np.float64],
        beta: NDArray[np.uint8],
    ) -> NDArray[np.float64]:
        self.last_work += int(first.size)
        return np.where(beta == 0, first + second, second - first)

    def _select(
        self,
        decisions: NDArray[np.uint8],
        metric: NDArray[np.float64],
    ) -> BitArray:
        order = np.argsort(metric, kind="stable")
        if self.crc is None:
            return decisions[order[0]][self.information_indices][
                : self._data_bits
            ]
        for index in order:
            information = decisions[index][self.information_indices]
            if self.crc.check(information):
                return information[: self._data_bits]
        return decisions[order[0]][self.information_indices][: self._data_bits]


__all__ = [
    "PolarAWGN",
    "PolarBEC",
    "awgn_synthetic_llr_means",
    "bec_capacity",
    "bec_synthetic_erasure_probabilities",
    "gaussian_phi",
    "gaussian_phi_inverse",
]
