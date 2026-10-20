# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""OFDM over multipath with pilot-based channel estimation.

The model is the frequency-domain equivalent of a cyclic-prefix OFDM link:

$$
y_m = H_m x_m + n_m,
\\qquad
H_m = \\sum_{\\ell} h_{\\ell} e^{-j 2\\pi m d_{\\ell}/N},
$$

with one flat gain per subcarrier and circular complex AWGN.  The cyclic
prefix is assumed to be at least as long as the delay spread, so there is no
inter-symbol interference, and the taps are constant within one OFDM block, so
intra-block Doppler (ICI) is not modelled.  Doppler appears as block-to-block
tap variation.

The receiver estimates the channel from comb pilots (least squares, then
periodic linear interpolation), equalises with zero forcing, and forms LLRs
with the per-subcarrier post-equalisation variance.  This is a
*channel-specific* simulation path: :class:`modcode.pipeline.Link` cannot use
it, because its LLRs assume one scalar variance for the whole block.
"""

from __future__ import annotations

from dataclasses import dataclass
from numbers import Integral

import numpy as np
from numpy.typing import NDArray

from modcode.channel.awgn import AWGNChannel
from modcode.channel.multipath import TapProcess, subcarrier_gains
from modcode.fec.base import FEC
from modcode.modulation.base import BitArray, Modulation


@dataclass(frozen=True)
class OFDMResult:
    """Captured bits and channel estimates from one OFDM link simulation."""

    source_bits: BitArray
    received_bits: BitArray
    channel_nmse: float
    gains: NDArray[np.complex128]
    estimates: NDArray[np.complex128]
    decoding: str = "hard"
    estimation: str = "ls"

    @property
    def bit_errors(self) -> int:
        return int(np.count_nonzero(self.source_bits != self.received_bits))

    @property
    def ber(self) -> float:
        if self.source_bits.size == 0:
            return 0.0
        return self.bit_errors / self.source_bits.size


@dataclass(frozen=True)
class OFDMLink:
    """An OFDM link: modulation, multipath tap process, pilots, and FEC."""

    modulation: Modulation
    channel: TapProcess
    delays: tuple[int, ...]
    fec: FEC | None = None
    n_fft: int = 64
    pilot_spacing: int = 4
    estimation: str = "ls"

    def __post_init__(self) -> None:
        if (
            isinstance(self.n_fft, bool)
            or not isinstance(self.n_fft, Integral)
            or self.n_fft < 2
        ):
            raise ValueError("n_fft must be an integer of at least 2")
        if self.estimation not in {"ls", "genie"}:
            raise ValueError("estimation must be 'ls' or 'genie'")
        if len(self.delays) != self.channel.taps_count:
            raise ValueError("delays must match the tap process length")
        pilots = self.pilot_subcarriers
        if pilots.size >= self.n_fft:
            raise ValueError("at least one data subcarrier is required")

    @property
    def pilot_subcarriers(self) -> NDArray[np.intp]:
        """Return the comb-pilot subcarrier indices (starting at 0)."""
        if (
            isinstance(self.pilot_spacing, bool)
            or not isinstance(self.pilot_spacing, Integral)
            or self.pilot_spacing < 1
        ):
            raise ValueError("pilot_spacing must be a positive integer")
        return np.arange(0, self.n_fft, int(self.pilot_spacing), dtype=np.intp)

    @property
    def data_subcarriers(self) -> NDArray[np.intp]:
        """Return the data subcarrier indices, i.e. the complement of pilots."""
        return np.setdiff1d(
            np.arange(self.n_fft, dtype=np.intp), self.pilot_subcarriers
        )

    @property
    def effective_rate(self) -> float:
        """Coded bits per subcarrier, counting the pilot overhead."""
        rate = self.fec.rate if self.fec is not None else 1.0
        return rate * self.data_subcarriers.size / self.n_fft


def interpolate_channel(
    pilot_gains: NDArray[np.complexfloating],
    pilot_indices: NDArray[np.integer],
    n_fft: int,
) -> NDArray[np.complex128]:
    """Interpolate pilot gains over all subcarriers, periodically in frequency."""
    values = np.asarray(pilot_gains, dtype=np.complex128)
    indices = np.asarray(pilot_indices, dtype=np.intp)
    if values.ndim != 1 or values.size != indices.size or values.size < 1:
        raise ValueError("pilot_gains must match pilot_indices")
    if indices.size == 1:
        return np.full(int(n_fft), values[0], dtype=np.complex128)
    folded = np.concatenate(
        (indices - int(n_fft), indices, indices + int(n_fft))
    )
    repeated = np.tile(values, 3)
    grid = np.arange(int(n_fft))
    return np.interp(grid, folded, repeated.real) + 1j * np.interp(
        grid, folded, repeated.imag
    )


def run_ofdm_link(
    link: OFDMLink,
    n_bits: int = 20_000,
    ebn0_db: float = 10.0,
    rng: np.random.Generator | None = None,
) -> OFDMResult:
    """Simulate one OFDM link with pilot estimation and zero forcing."""
    if n_bits <= 0:
        raise ValueError("n_bits must be positive")
    if not np.isfinite(ebn0_db):
        raise ValueError("ebn0_db must be finite")
    generator = rng if rng is not None else np.random.default_rng()
    pilots = link.pilot_subcarriers
    data = link.data_subcarriers
    bits_per_symbol = link.modulation.bits_per_symbol
    coded_per_block = data.size * bits_per_symbol
    variance = AWGNChannel.noise_variance(
        ebn0_db, bits_per_symbol, link.effective_rate
    )

    source = generator.integers(0, 2, n_bits, dtype=np.uint8)
    if link.fec is None:
        encoded = np.concatenate(
            (source, np.zeros((-n_bits) % coded_per_block, dtype=np.uint8))
        )
        stride = coded_per_block
        decoding = "hard"
    else:
        data_bits = link.fec.data_bits_per_block
        padded = np.concatenate(
            (source, np.zeros((-n_bits) % data_bits, dtype=np.uint8))
        )
        encoded = link.fec.encode(padded)
        stride = int(np.lcm(coded_per_block, link.fec.encoded_bits_per_block))
        decoding = "soft"
    encoded = np.concatenate(
        (encoded, np.zeros((-encoded.size) % stride, dtype=np.uint8))
    )

    blocks = encoded.size // coded_per_block
    gains = np.empty((blocks, link.n_fft), dtype=np.complex128)
    estimates = np.empty((blocks, link.n_fft), dtype=np.complex128)
    llrs: list[NDArray[np.float64]] = []
    for block in range(blocks):
        taps = link.channel.taps(1)[0]
        gains[block] = subcarrier_gains(taps, link.delays, link.n_fft)
        bits = encoded[block * coded_per_block : (block + 1) * coded_per_block]
        transmitted = np.zeros(link.n_fft, dtype=np.complex128)
        transmitted[pilots] = 1.0 + 0.0j
        transmitted[data] = link.modulation.modulate(bits)
        received = gains[block] * transmitted + np.sqrt(variance) * (
            generator.standard_normal(link.n_fft)
            + 1j * generator.standard_normal(link.n_fft)
        )

        if link.estimation == "genie":
            estimate = gains[block]
        else:
            pilot_gains = received[pilots] / transmitted[pilots]
            estimate = interpolate_channel(pilot_gains, pilots, link.n_fft)
        estimates[block] = estimate

        equalised = received[data] / estimate[data]
        subcarrier_variance = variance / np.abs(estimate[data]) ** 2
        base = link.modulation.soft_demodulate(equalised, 1.0)
        base = base.reshape(data.size, bits_per_symbol)
        llrs.append((base / subcarrier_variance[:, None]).reshape(-1))

    llr_stream = np.concatenate(llrs)
    if link.fec is None:
        received_bits = (llr_stream < 0.0).astype(np.uint8)
    else:
        received_bits = link.fec.decode_soft(llr_stream)
    nmse = float(
        np.sum(np.abs(estimates - gains) ** 2) / np.sum(np.abs(gains) ** 2)
    )
    return OFDMResult(
        source_bits=source,
        received_bits=received_bits[:n_bits],
        channel_nmse=nmse,
        gains=gains,
        estimates=estimates,
        decoding=decoding,
        estimation=link.estimation,
    )


__all__ = ["OFDMLink", "OFDMResult", "interpolate_channel", "run_ofdm_link"]
