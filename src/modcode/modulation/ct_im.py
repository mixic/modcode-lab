# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Cyclic Trellis-Constrained Index Modulation (CT-IM) and its baselines.

The index pattern of a CT-IM group is a codeword of a constant-weight cyclic
code. For length 7, the weight-3 codewords of Hamming(7,4) form exactly one
cyclic orbit of size 7; they have minimum Hamming distance 4 and therefore
correct any single index-bit error.  The module also contains the shortened
constant-weight variant, a BCH(15,7) weight-7 orbit, and the classical
OFDM-IM baseline used for the BER comparison.

The name "trellis-constrained" is inherited from the research proposal; the
implemented constraint is the cyclic constant-weight codebook.  A minimal
trellis realization is left as future work.
"""

from __future__ import annotations

import itertools
import math
from abc import ABC, abstractmethod

import numpy as np
from numpy.typing import ArrayLike, NDArray

from modcode.fec.bch import BCH15
from modcode.fec.hamming import Hamming74
from modcode.modulation.qam import QAM


def _validated_bits(bits: ArrayLike, name: str = "bits") -> NDArray[np.uint8]:
    values = np.asarray(bits)
    if values.ndim != 1:
        raise ValueError(f"{name} must be a one-dimensional array")
    if not np.all((values == 0) | (values == 1)):
        raise ValueError(f"{name} must contain only 0 and 1")
    return values.astype(np.uint8, copy=False)


class IndexCode:
    """A constant-weight codebook used as the index alphabet."""

    n: int
    weight: int
    minimum_distance: int

    def __init__(self, patterns: NDArray[np.uint8]) -> None:
        values = np.asarray(patterns, dtype=np.uint8)
        if values.ndim != 2 or values.shape[0] == 0:
            raise ValueError("patterns must be a non-empty two-dimensional array")
        if not np.all((values == 0) | (values == 1)):
            raise ValueError("patterns must contain only 0 and 1")
        if not np.all(values.sum(axis=1) == values.sum(axis=1)[0]):
            raise ValueError("patterns must have constant weight")
        self.patterns = values
        self.n = int(values.shape[1])
        self.weight = int(values.sum(axis=1)[0])

    def __len__(self) -> int:
        return self.patterns.shape[0]

    @property
    def index_bits(self) -> int:
        """Fixed-length label width: ceil(log2(number of patterns))."""
        return int(math.ceil(math.log2(len(self))))

    def encode(self, label: int) -> NDArray[np.uint8]:
        """Return the index pattern for a 0-based label."""
        if isinstance(label, bool) or not isinstance(label, int):
            raise ValueError("label must be an integer")
        if not 0 <= label < len(self):
            raise ValueError(f"label must be in [0, {len(self)})")
        return self.patterns[label].copy()

    def decode(self, bits: ArrayLike) -> tuple[int, int]:
        """Return ``(nearest_label, hamming_distance)`` for a hard pattern."""
        values = _validated_bits(bits, "index pattern")
        if values.size != self.n:
            raise ValueError(f"index pattern must have length {self.n}")
        distances = np.sum(values != self.patterns, axis=1)
        label = int(np.argmin(distances))
        return label, int(distances[label])

    def is_codeword(self, bits: ArrayLike) -> bool:
        """Return whether ``bits`` is exactly one of the patterns."""
        values = _validated_bits(bits, "index pattern")
        if values.size != self.n:
            return False
        return bool(np.any(np.all(self.patterns == values, axis=1)))

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


class ConstantWeightCyclicCode(IndexCode):
    """The weight-3 cyclic orbit of Hamming(7,4).

    The seven weight-3 Hamming codewords are closed under cyclic shifts and
    have pairwise Hamming distance 4. The set is not a linear subspace (it has
    no zero word), but it is a constant-weight cyclic code in the coding-theory
    sense of a cyclic-invariant constant-weight codebook.
    """

    def __init__(self) -> None:
        hamming = Hamming74()
        data = np.asarray(
            list(itertools.product((0, 1), repeat=4)), dtype=np.uint8
        )
        words = {tuple(hamming.encode(bits).tolist()) for bits in data}
        patterns = sorted(word for word in words if sum(word) == 3)
        if len(patterns) != 7:
            raise RuntimeError("unexpected Hamming weight distribution")
        super().__init__(np.asarray(patterns, dtype=np.uint8))
        self.minimum_distance = 4


class ShortenedConstantWeightCode(IndexCode):
    """CT-IM-1: the constant-weight shortening of the Hamming(7,4) orbit.

    Puncturing the weight-3 orbit does **not** keep a constant weight: deleting
    an active position gives weight 2 and deleting a silent position gives
    weight 3.  Shortening does: fix coordinate 0 to zero and delete it, which
    keeps the four words that are silent there.  The result is a length-6,
    weight-3 code with pairwise distance 4, i.e. a constant symbol budget of
    three active subcarriers per block.
    """

    def __init__(self) -> None:
        orbit = ConstantWeightCyclicCode()
        kept = [word for word in orbit.patterns if word[0] == 0]
        if len(kept) != 4:
            raise RuntimeError("unexpected shortened codebook size")
        patterns = sorted(tuple(word[1:].tolist()) for word in kept)
        super().__init__(np.asarray(patterns, dtype=np.uint8))
        self.minimum_distance = 4
        self.shortened_position = 0


class BCHIndexCode(IndexCode):
    """CT-IM-6: a constant-weight cyclic orbit of BCH(15,7).

    The binary BCH(15,7,d=5) code is linear, so its codeword weights vary; it
    is not constant-weight.  The weight-7 codewords form a single cyclic orbit
    of 15 words with pairwise distance 8, and that orbit is used here as a
    length-15 constant-weight index alphabet.
    """

    def __init__(self, weight: int = 7) -> None:
        bch = BCH15()
        words = np.unique(bch._codewords, axis=0)
        subset = words[words.sum(axis=1) == weight]
        if subset.shape[0] == 0:
            raise ValueError(f"BCH(15,7) has no weight-{weight} codewords")
        index = {tuple(row.tolist()): i for i, row in enumerate(subset)}
        best: list[tuple[int, ...]] = []
        seen: set[int] = set()
        for start, row in enumerate(subset):
            if start in seen:
                continue
            orbit = {
                np.roll(row, steps).tobytes(): tuple(np.roll(row, steps).tolist())
                for steps in range(bch.n)
            }
            members = sorted(set(orbit.values()))
            seen.update(index[member] for member in members)
            if len(members) > len(best):
                best = members
        if not best:
            raise RuntimeError("no cyclic orbit found")
        array = np.asarray(sorted(set(best)), dtype=np.uint8)
        super().__init__(array)
        distances = [
            int(np.sum(array[left] != array[right]))
            for left in range(array.shape[0])
            for right in range(left + 1, array.shape[0])
        ]
        self.minimum_distance = min(distances)


def _nearest_symbol_distances(
    received: NDArray[np.complex128],
    gains: NDArray[np.complex128] | None,
    constellation: NDArray[np.complex128],
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Return per-subcarrier squared distances to silent and active symbols."""
    if gains is not None:
        scaled = gains[:, None] * constellation[None, :]
        distances = np.abs(received[:, None] - scaled) ** 2
    else:
        distances = np.abs(received[:, None] - constellation[None, :]) ** 2
    return np.abs(received) ** 2, distances.min(axis=1)


def ml_label(
    received: ArrayLike,
    patterns: NDArray[np.uint8],
    constellation: NDArray[np.complex128],
    *,
    gains: NDArray[np.complexfloating] | None = None,
) -> int:
    """Return the maximum-likelihood index label for known CSI.

    The metric is the squared distance to the nearest active symbol on the
    pattern positions plus the squared distance to silence elsewhere.  For
    equal-energy symbols this is the exact ML index detector; the noise
    variance cancels because it is common to all patterns.
    """
    values = np.asarray(received, dtype=np.complex128)
    gain_values = (
        None if gains is None else np.asarray(gains, dtype=np.complex128)
    )
    if gain_values is not None and gain_values.shape != values.shape:
        raise ValueError("gains must match received")
    silent, active = _nearest_symbol_distances(
        values, gain_values, constellation
    )
    mask = patterns.astype(bool)
    metrics = mask.astype(np.float64) @ active + (~mask).astype(
        np.float64
    ) @ silent
    return int(np.argmin(metrics))


class IndexModulation(ABC):
    """Common mapper for a constant-weight index codebook with QAM symbols."""

    def __init__(self, code: IndexCode, M: int = 4) -> None:
        self.index_code = code
        self.symbol_mod = QAM(M)
        self.M = int(M)
        self.n = code.n
        self.active = code.weight
        self.index_bits = code.index_bits
        self.symbol_bits = self.active * self.symbol_mod.bits_per_symbol
        self.bits_per_block = self.index_bits + self.symbol_bits

    @property
    def spectral_efficiency(self) -> float:
        """Bits per subcarrier, counting the fixed-length index label."""
        return self.bits_per_block / self.n

    @property
    def energy_per_block(self) -> int:
        """Average symbol energy per block (one unit-energy symbol per active)."""
        return self.active

    def _label_from_bits(self, bits: NDArray[np.uint8]) -> int:
        label = 0
        for bit in bits[: self.index_bits]:
            label = (label << 1) | int(bit)
        if label >= len(self.index_code):
            raise ValueError(
                f"index label {label} is reserved: "
                f"the codebook has {len(self.index_code)} words"
            )
        return label

    def _encode_parts(
        self, bits: NDArray[np.uint8]
    ) -> tuple[int, NDArray[np.uint8]]:
        values = _validated_bits(bits)
        if values.size != self.bits_per_block:
            raise ValueError(
                f"expected {self.bits_per_block} bits, got {values.size}"
            )
        label = self._label_from_bits(values)
        return label, values[self.index_bits :]

    def _pattern_vector(
        self, pattern: NDArray[np.uint8], symbols: NDArray[np.complex128]
    ) -> NDArray[np.complex128]:
        vector = np.zeros(self.n, dtype=np.complex128)
        vector[pattern.astype(bool)] = symbols
        return vector

    def _bits_from_label(self, label: int) -> NDArray[np.uint8]:
        return np.asarray(
            [(label >> shift) & 1 for shift in range(self.index_bits - 1, -1, -1)],
            dtype=np.uint8,
        )

    def _resolve_gains(
        self, received: NDArray[np.complex128], gains: ArrayLike | None
    ) -> NDArray[np.complex128] | None:
        if gains is None:
            return None
        values = np.asarray(gains, dtype=np.complex128)
        if values.shape != received.shape:
            raise ValueError("gains must match received")
        if not np.all(np.isfinite(values)) or np.any(values == 0.0):
            raise ValueError("gains must be finite and non-zero")
        return values

    @staticmethod
    def _validate_received(received: ArrayLike, n: int) -> NDArray[np.complex128]:
        values = np.asarray(received, dtype=np.complex128)
        if values.ndim != 1 or values.size != n:
            raise ValueError(f"received must be a one-dimensional length-{n} vector")
        if not np.all(np.isfinite(values)):
            raise ValueError("received must contain only finite values")
        return values

    @abstractmethod
    def encode(self, bits: ArrayLike) -> NDArray[np.complex128]:
        """Map one block of bits to a frequency-domain vector."""

    @abstractmethod
    def decode(
        self,
        received: ArrayLike,
        *,
        gains: ArrayLike | None = None,
    ) -> NDArray[np.uint8]:
        """Recover the bit block from a frequency-domain observation."""


class CTIM(IndexModulation):
    """Cyclic Trellis-Constrained Index Modulation for one index group.

    Parameters
    ----------
    M
        QAM order. ``M=4`` is QPSK.
    code
        Index codebook; default is the seven-word weight-3 Hamming orbit.
    cyclic_shift
        If true, ``encode`` selects the cyclic shift with the lowest PAPR. The
        receiver must be told the shift; the syndrome cannot identify it,
        because every shift of a codeword is again a codeword.
    detector
        ``"threshold"`` uses the fixed magnitude threshold on ``|y_j|``;
        ``"ml"`` uses the exact maximum-likelihood index detector over the
        codebook.  The ML detector is the matched baseline for the experiments.
    threshold
        Activity-detection threshold on ``|y_j|`` for the threshold detector.
    """

    def __init__(
        self,
        M: int = 4,
        *,
        code: IndexCode | None = None,
        cyclic_shift: bool = False,
        detector: str = "threshold",
        threshold: float = 0.5,
    ) -> None:
        super().__init__(code if code is not None else ConstantWeightCyclicCode(), M)
        self.cyclic_shift = bool(cyclic_shift)
        if detector not in {"threshold", "ml"}:
            raise ValueError("detector must be 'threshold' or 'ml'")
        self.detector = detector
        self.threshold = float(threshold)
        if not np.isfinite(self.threshold) or self.threshold <= 0.0:
            raise ValueError("threshold must be finite and positive")

    @property
    def index_rate(self) -> float:
        """Index bits per subcarrier."""
        return self.index_bits / self.n

    def encode(
        self, bits: ArrayLike, *, return_shift: bool = False
    ) -> NDArray[np.complex128] | tuple[NDArray[np.complex128], int]:
        """Map one block of bits to a frequency-domain CT-IM vector."""
        label, symbol_bits = self._encode_parts(np.asarray(bits, dtype=np.uint8))
        pattern = self.index_code.encode(label)
        symbols = self.symbol_mod.modulate(symbol_bits)
        vector = self._pattern_vector(pattern, symbols)
        shift = self.best_shift(vector) if self.cyclic_shift else 0
        if shift:
            vector = np.roll(vector, shift)
        if return_shift:
            return vector, shift
        return vector

    def _detect_ml(
        self, received: NDArray[np.complex128], gains: NDArray[np.complex128] | None
    ) -> int:
        constellation = np.asarray(
            list(self.symbol_mod.constellation.values()), dtype=np.complex128
        )
        return ml_label(
            received, self.index_code.patterns, constellation, gains=gains
        )

    def _detect_threshold(
        self, received: NDArray[np.complex128], gains: NDArray[np.complex128] | None
    ) -> tuple[int, int]:
        scaled = np.abs(received) if gains is None else np.abs(received) / np.abs(
            gains
        )
        activity = (scaled > self.threshold).astype(np.uint8)
        return self.index_code.decode(activity)

    def decode_index(
        self,
        received: ArrayLike,
        *,
        shift: int = 0,
        gains: ArrayLike | None = None,
    ) -> tuple[int, int]:
        """Return the index label and its Hamming distance to the observation."""
        values = self._validate_received(received, self.n)
        if isinstance(shift, bool) or not isinstance(shift, int):
            raise ValueError("shift must be an integer")
        if shift:
            values = np.roll(values, -shift)
        gain_values = self._resolve_gains(values, gains)
        if self.detector == "ml":
            label = self._detect_ml(values, gain_values)
            pattern = self.index_code.encode(label)
            activity = self._activity(values, gain_values)
            return label, int(np.sum(pattern != activity))
        return self._detect_threshold(values, gain_values)

    def _activity(
        self, received: NDArray[np.complex128], gains: NDArray[np.complex128] | None
    ) -> NDArray[np.uint8]:
        scaled = np.abs(received) if gains is None else np.abs(received) / np.abs(
            gains
        )
        return (scaled > self.threshold).astype(np.uint8)

    def decode(
        self,
        received: ArrayLike,
        *,
        shift: int = 0,
        gains: ArrayLike | None = None,
    ) -> NDArray[np.uint8]:
        """Recover the bit block from a frequency-domain observation."""
        values = self._validate_received(received, self.n)
        if isinstance(shift, bool) or not isinstance(shift, int):
            raise ValueError("shift must be an integer")
        if shift:
            values = np.roll(values, -shift)
        gain_values = self._resolve_gains(values, gains)
        if self.detector == "ml":
            label = self._detect_ml(values, gain_values)
        else:
            label, _ = self._detect_threshold(values, gain_values)
        pattern = self.index_code.encode(label)
        observed = (
            values if gain_values is None else values / gain_values
        )
        symbols = observed[pattern.astype(bool)]
        symbol_bits = self.symbol_mod.demodulate(symbols)
        return np.concatenate(
            [self._bits_from_label(label), symbol_bits]
        ).astype(np.uint8)

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

    def best_shift(
        self, vector: ArrayLike, shifts: tuple[int, ...] | None = None
    ) -> int:
        """Return the shift with the lowest PAPR, within ``shifts`` if given."""
        values = np.asarray(vector, dtype=np.complex128)
        candidates = range(self.n) if shifts is None else shifts
        scores = [self.papr_db(np.roll(values, shift)) for shift in candidates]
        return int(list(candidates)[int(np.argmin(scores))])


class ClassicalOFDMIM(IndexModulation):
    """CT-IM-2 baseline: classical OFDM-IM with n=7, k=3, QPSK.

    There are 35 active-set patterns; the twelve that carry the 5-bit index
    alphabet are the first 32 in lexicographic order.  With 6 symbol bits the
    block carries 11 bits, and the index detector is the maximum-likelihood
    detector over the 32-pattern codebook (the matched baseline).
    """

    n = 7
    active = 3
    index_bits = 5
    bits_per_block = 11

    def __init__(self, M: int = 4) -> None:
        combinations = list(itertools.combinations(range(7), 3))[:32]
        patterns = np.zeros((len(combinations), 7), dtype=np.uint8)
        for row, active in enumerate(combinations):
            patterns[row, list(active)] = 1
        self.index_code = IndexCode(patterns)
        self.index_code.minimum_distance = int(
            min(
                np.sum(patterns[left] != patterns[right])
                for left in range(32)
                for right in range(left + 1, 32)
            )
        )
        self.symbol_mod = QAM(M)
        self.M = int(M)
        self.symbol_bits = 3 * self.symbol_mod.bits_per_symbol
        if self.symbol_bits + self.index_bits != self.bits_per_block:
            raise RuntimeError("classical OFDM-IM budget mismatch")

    def encode(self, bits: ArrayLike) -> NDArray[np.complex128]:
        label, symbol_bits = self._encode_parts(np.asarray(bits, dtype=np.uint8))
        pattern = self.index_code.patterns[label]
        return self._pattern_vector(pattern, self.symbol_mod.modulate(symbol_bits))

    def decode(
        self,
        received: ArrayLike,
        *,
        gains: ArrayLike | None = None,
    ) -> NDArray[np.uint8]:
        values = self._validate_received(received, self.n)
        gain_values = self._resolve_gains(values, gains)
        constellation = np.asarray(
            list(self.symbol_mod.constellation.values()), dtype=np.complex128
        )
        label = ml_label(
            values, self.index_code.patterns, constellation, gains=gain_values
        )
        pattern = self.index_code.patterns[label]
        observed = values if gain_values is None else values / gain_values
        symbol_bits = self.symbol_mod.demodulate(observed[pattern.astype(bool)])
        return np.concatenate([self._bits_from_label(label), symbol_bits]).astype(
            np.uint8
        )


__all__ = [
    "BCHIndexCode",
    "CTIM",
    "ClassicalOFDMIM",
    "ConstantWeightCyclicCode",
    "IndexCode",
    "IndexModulation",
    "ShortenedConstantWeightCode",
    "ml_label",
]
