# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Cyclic Trellis-Constrained Index Modulation (CT-IM).

The index pattern of a CT-IM group is a codeword of a constant-weight cyclic
code. For length 7, the weight-3 codewords of Hamming(7,4) form exactly one
cyclic orbit of size 7; they have minimum Hamming distance 4 and therefore
correct any single index-bit error.

This module implements the codebook, the OFDM-IM mapper, and a threshold
detector. The name "trellis-constrained" is inherited from the research
proposal; the implemented constraint is the cyclic constant-weight codebook.
A minimal trellis realization is left as future work.
"""

from __future__ import annotations

import itertools

import numpy as np
from numpy.typing import ArrayLike, NDArray

from modcode.fec.hamming import Hamming74
from modcode.modulation.qam import QAM


def _validated_bits(bits: ArrayLike, name: str = "bits") -> NDArray[np.uint8]:
    values = np.asarray(bits)
    if values.ndim != 1:
        raise ValueError(f"{name} must be a one-dimensional array")
    if not np.all((values == 0) | (values == 1)):
        raise ValueError(f"{name} must contain only 0 and 1")
    return values.astype(np.uint8, copy=False)


class ConstantWeightCyclicCode:
    """The weight-3 cyclic orbit of Hamming(7,4).

    The seven weight-3 Hamming codewords are closed under cyclic shifts and
    have pairwise Hamming distance 4. The set is not a linear subspace (it has
    no zero word), but it is a constant-weight cyclic code in the coding-theory
    sense of a cyclic-invariant constant-weight codebook.
    """

    def __init__(self) -> None:
        self.n = 7
        self.weight = 3
        self.minimum_distance = 4
        hamming = Hamming74()
        data = np.asarray(
            list(itertools.product((0, 1), repeat=4)), dtype=np.uint8
        )
        words = {
            tuple(hamming.encode(bits).tolist()) for bits in data
        }
        patterns = sorted(word for word in words if sum(word) == self.weight)
        if len(patterns) != 7:
            raise RuntimeError("unexpected Hamming weight distribution")
        self.patterns = np.asarray(patterns, dtype=np.uint8)
        self._labels = {
            pattern: label for label, pattern in enumerate(patterns)
        }

    def __len__(self) -> int:
        return self.patterns.shape[0]

    def encode(self, label: int) -> NDArray[np.uint8]:
        """Return the index pattern for a 0-based label in ``[0, 7)``."""
        if isinstance(label, bool) or not isinstance(label, int):
            raise ValueError("label must be an integer")
        if not 0 <= label < len(self):
            raise ValueError("label must be in [0, 7)")
        return self.patterns[label].copy()

    def decode(self, bits: ArrayLike) -> tuple[int, int]:
        """Return ``(nearest_label, hamming_distance)`` for a hard pattern."""
        values = _validated_bits(bits, "index pattern")
        if values.size != self.n:
            raise ValueError("index pattern must have length 7")
        distances = np.sum(values != self.patterns, axis=1)
        label = int(np.argmin(distances))
        return label, int(distances[label])

    def is_codeword(self, bits: ArrayLike) -> bool:
        """Return whether ``bits`` is exactly one of the seven patterns."""
        values = _validated_bits(bits, "index pattern")
        return tuple(values.tolist()) in self._labels

    def shift(self, label: int, steps: int) -> NDArray[np.uint8]:
        """Return a cyclically shifted codeword."""
        if isinstance(steps, bool) or not isinstance(steps, int):
            raise ValueError("steps must be an integer")
        return np.roll(self.encode(label), steps)

    def distance_matrix(self) -> NDArray[np.int64]:
        """Return the pairwise Hamming distance matrix of the codebook."""
        return np.asarray(
            [
                [np.sum(left != right) for right in self.patterns]
                for left in self.patterns
            ],
            dtype=np.int64,
        )


class CTIM:
    """Cyclic Trellis-Constrained Index Modulation for one 7-subcarrier group.

    Parameters
    ----------
    M
        QAM order. ``M=4`` is QPSK.
    cyclic_shift
        If true, ``encode`` selects the cyclic shift with the lowest PAPR. The
        receiver must be told the shift; the syndrome cannot identify it,
        because every shift of a codeword is again a codeword.
    threshold
        Activity-detection threshold on ``|y_j|``. The default is suitable
        for unit-energy QPSK; other constellations need calibration.
    """

    n = 7
    active = 3
    index_bits = 3

    def __init__(
        self,
        M: int = 4,
        *,
        cyclic_shift: bool = False,
        threshold: float = 0.5,
    ) -> None:
        self.index_code = ConstantWeightCyclicCode()
        self.symbol_mod = QAM(M)
        self.M = int(M)
        self.cyclic_shift = cyclic_shift
        self.threshold = float(threshold)
        if not np.isfinite(self.threshold) or self.threshold <= 0.0:
            raise ValueError("threshold must be finite and positive")
        self.symbol_bits = self.active * self.symbol_mod.bits_per_symbol
        self.bits_per_block = self.index_bits + self.symbol_bits

    @property
    def spectral_efficiency(self) -> float:
        """Bits per subcarrier, counting the 3-bit index label."""
        return self.bits_per_block / self.n

    @property
    def index_rate(self) -> float:
        """Index bits per subcarrier (3-bit label, 8 labels, 7 valid)."""
        return self.index_bits / self.n

    def encode(
        self, bits: ArrayLike, *, return_shift: bool = False
    ) -> NDArray[np.complex128] | tuple[NDArray[np.complex128], int]:
        """Map one block of bits to a frequency-domain CT-IM vector.

        The first three bits select one of the seven index labels. The value
        ``7`` is reserved; it has no codeword and is rejected.
        """
        values = _validated_bits(bits)
        if values.size != self.bits_per_block:
            raise ValueError(
                f"expected {self.bits_per_block} bits, got {values.size}"
            )
        label = int(
            (values[0] << 2) | (values[1] << 1) | values[2]
        )
        if label >= len(self.index_code):
            raise ValueError("index label 7 is reserved")

        pattern = self.index_code.encode(label)
        symbols = self.symbol_mod.modulate(values[self.index_bits :])
        vector = np.zeros(self.n, dtype=np.complex128)
        vector[pattern.astype(bool)] = symbols

        shift = self.best_shift(vector) if self.cyclic_shift else 0
        if shift:
            vector = np.roll(vector, shift)
        if return_shift:
            return vector, shift
        return vector

    def decode(self, received: ArrayLike, *, shift: int = 0) -> NDArray[np.uint8]:
        """Recover the bit block from a frequency-domain observation."""
        values = np.asarray(received, dtype=np.complex128)
        if values.ndim != 1 or values.size != self.n:
            raise ValueError("received must be a one-dimensional length-7 vector")
        if not np.all(np.isfinite(values)):
            raise ValueError("received must contain only finite values")
        if isinstance(shift, bool) or not isinstance(shift, int):
            raise ValueError("shift must be an integer")
        if shift:
            values = np.roll(values, -shift)

        activity = (np.abs(values) > self.threshold).astype(np.uint8)
        label, _ = self.index_code.decode(activity)
        pattern = self.index_code.encode(label)

        symbols = values[pattern.astype(bool)]
        symbol_bits = self.symbol_mod.demodulate(symbols)
        index_bits = np.asarray(
            [(label >> 2) & 1, (label >> 1) & 1, label & 1],
            dtype=np.uint8,
        )
        return np.concatenate([index_bits, symbol_bits])

    def decode_index(
        self, received: ArrayLike, *, shift: int = 0
    ) -> tuple[int, int]:
        """Return the nearest index label and its distance."""
        values = np.asarray(received, dtype=np.complex128)
        if shift:
            values = np.roll(values, -shift)
        activity = (np.abs(values) > self.threshold).astype(np.uint8)
        return self.index_code.decode(activity)

    @staticmethod
    def papr_db(vector: ArrayLike) -> float:
        """Return the PAPR of the OFDM symbol in dB."""
        values = np.asarray(vector, dtype=np.complex128)
        waveform = np.fft.ifft(values)
        power = np.abs(waveform) ** 2
        mean_power = float(power.mean())
        if mean_power == 0.0:
            raise ValueError("vector must not be all zero")
        return float(10.0 * np.log10(power.max() / mean_power))

    def best_shift(self, vector: ArrayLike) -> int:
        """Return the cyclic shift with the lowest PAPR."""
        values = np.asarray(vector, dtype=np.complex128)
        scores = [
            self.papr_db(np.roll(values, shift)) for shift in range(self.n)
        ]
        return int(np.argmin(scores))


__all__ = ["CTIM", "ConstantWeightCyclicCode"]
