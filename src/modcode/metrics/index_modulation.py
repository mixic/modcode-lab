# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Monte Carlo BER and PAPR for constant-weight index modulation.

The simulator is shared by the CT-IM experiments so that every scheme uses the
same channel, the same Eb/N0 definition, and the same detector budget.  Eb/N0
is per information bit: a block carries ``index_bits + symbol_bits`` bits at
block energy ``energy_per_block``, so the noise variance per real component is
``energy / (2 * bits_per_block * Eb/N0)``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray


@dataclass(frozen=True)
class IndexResult:
    """Error counts for one index-modulation scheme at one operating point."""

    index_errors: int
    symbol_errors: int
    block_errors: int
    blocks: int
    index_bits: int
    symbol_bits: int

    @property
    def bits_per_block(self) -> int:
        return self.index_bits + self.symbol_bits

    @property
    def index_ber(self) -> float:
        total = self.blocks * self.index_bits
        return self.index_errors / total if total else 0.0

    @property
    def symbol_ber(self) -> float:
        total = self.blocks * self.symbol_bits
        return self.symbol_errors / total if total else 0.0

    @property
    def ber(self) -> float:
        total = self.blocks * self.bits_per_block
        errors = self.index_errors + self.symbol_errors
        return errors / total if total else 0.0

    @property
    def block_error_rate(self) -> float:
        return self.block_errors / self.blocks if self.blocks else 0.0


def _bits_of(label: int, width: int) -> NDArray[np.uint8]:
    return np.asarray(
        [(label >> shift) & 1 for shift in range(width - 1, -1, -1)],
        dtype=np.uint8,
    )


def simulate_index_modulation(
    scheme: object,
    ebn0_db: float,
    *,
    blocks: int = 2000,
    rayleigh: bool = False,
    rng: np.random.Generator | None = None,
) -> IndexResult:
    """Simulate one index-modulation scheme over AWGN or known-CSI Rayleigh.

    ``scheme`` must expose ``n``, ``bits_per_block``, ``index_bits``,
    ``energy_per_block``, ``index_code``, ``encode(bits)``, and
    ``decode(received, gains=...)``, which is the interface of
    :class:`modcode.modulation.ct_im.CTIM` and
    :class:`modcode.modulation.ct_im.ClassicalOFDMIM`.
    """
    if not np.isfinite(ebn0_db):
        raise ValueError("ebn0_db must be finite")
    if blocks <= 0:
        raise ValueError("blocks must be positive")
    generator = rng if rng is not None else np.random.default_rng()
    n = int(scheme.n)  # type: ignore[attr-defined]
    bits_per_block = int(scheme.bits_per_block)  # type: ignore[attr-defined]
    index_bits = int(scheme.index_bits)  # type: ignore[attr-defined]
    symbol_bits = bits_per_block - index_bits
    energy = float(scheme.energy_per_block)  # type: ignore[attr-defined]
    codebook = scheme.index_code  # type: ignore[attr-defined]
    ebn0 = 10.0 ** (ebn0_db / 10.0)
    variance = energy / (2.0 * bits_per_block * ebn0)

    index_errors = 0
    symbol_errors = 0
    block_errors = 0
    for _ in range(int(blocks)):
        label = int(generator.integers(0, len(codebook)))
        bits = np.concatenate(
            (
                _bits_of(label, index_bits),
                generator.integers(0, 2, symbol_bits, dtype=np.uint8),
            )
        ).astype(np.uint8)
        vector = np.asarray(scheme.encode(bits), dtype=np.complex128)  # type: ignore[attr-defined]
        if rayleigh:
            gains = (
                generator.standard_normal(n)
                + 1j * generator.standard_normal(n)
            ) / np.sqrt(2.0)
            received = gains * vector
        else:
            gains = None
            received = vector.copy()
        received = received + np.sqrt(variance) * (
            generator.standard_normal(n) + 1j * generator.standard_normal(n)
        )
        decoded = np.asarray(
            scheme.decode(received, gains=gains),  # type: ignore[attr-defined]
            dtype=np.uint8,
        )
        errors = decoded != bits
        index_errors += int(errors[:index_bits].sum())
        symbol_errors += int(errors[index_bits:].sum())
        block_errors += int(errors.any())
    return IndexResult(
        index_errors=index_errors,
        symbol_errors=symbol_errors,
        block_errors=block_errors,
        blocks=int(blocks),
        index_bits=index_bits,
        symbol_bits=symbol_bits,
    )


def papr_values(vectors: ArrayLike) -> NDArray[np.float64]:
    """Return the PAPR in dB for each row of a stack of frequency vectors."""
    values = np.asarray(vectors, dtype=np.complex128)
    if values.ndim != 2:
        raise ValueError("vectors must be a two-dimensional array")
    waveform = np.fft.ifft(values, axis=1)
    power = np.abs(waveform) ** 2
    return 10.0 * np.log10(power.max(axis=1) / power.mean(axis=1))


def block_papr_values(
    vectors: ArrayLike,
    *,
    n_fft: int = 64,
    offset: int = 0,
) -> NDArray[np.float64]:
    """Return the PAPR of each block placed at ``offset`` in an ``n_fft`` band.

    A cyclic shift of the seven-subcarrier block inside its own length does not
    change the PAPR, because it multiplies the time waveform by a phase.  A
    shift of the block *within a larger band* does, and that is what this
    function measures.  ``offset`` is the first occupied bin.
    """
    values = np.asarray(vectors, dtype=np.complex128)
    if values.ndim != 2:
        raise ValueError("vectors must be a two-dimensional array")
    if not isinstance(n_fft, int) or n_fft < values.shape[1]:
        raise ValueError("n_fft must be at least the block length")
    if not 0 <= offset or offset + values.shape[1] > n_fft:
        raise ValueError("offset must place the whole block inside the band")
    grid = np.zeros((values.shape[0], n_fft), dtype=np.complex128)
    grid[:, offset : offset + values.shape[1]] = values
    waveform = np.fft.ifft(grid, axis=1)
    power = np.abs(waveform) ** 2
    return 10.0 * np.log10(power.max(axis=1) / power.mean(axis=1))


def papr_ccdf(
    vectors: ArrayLike, thresholds: NDArray[np.floating]
) -> NDArray[np.float64]:
    """Return the empirical CCDF ``P(PAPR > x)`` for a stack of vectors."""
    papr = papr_values(vectors)
    points = np.asarray(thresholds, dtype=np.float64)
    if points.ndim != 1:
        raise ValueError("thresholds must be one-dimensional")
    return np.asarray([np.mean(papr > point) for point in points])


__all__ = [
    "IndexResult",
    "block_papr_values",
    "papr_ccdf",
    "papr_values",
    "simulate_index_modulation",
]
