# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Multipath and Doppler channels: a symbol-spaced tap-delay line.

The channel is the discrete-time model

$$
y[k] = \\sum_{\\ell} h_{\\ell}[k]\\,x[k-d_{\\ell}] + n[k],
$$

where ``d_l`` are integer symbol delays and the tap gains are unit-average-power
Rayleigh processes.  Doppler is modelled by a first-order autoregressive
process whose one-step correlation is the classical Jakes coefficient

$$
\\rho = \\frac{1}{N}\\sum_{n=1}^{N}\\cos(2\\pi f_d T\\cos\\theta_n),
\\qquad \\theta_n = 2\\pi n/N,
$$

with ``T`` the update interval and ``f_d`` the maximum Doppler.  ``rho`` is one
for a static channel, decreases with ``f_d T``, and is zero for independent
updates.  The AR(1) model matches the Jakes correlation at one step and keeps
the average tap power exact; it is not the full Jakes spectrum.

This is the *time-domain* multipath model.  It produces inter-symbol
interference, so it is **not** a drop-in replacement for the AWGN channel in
:class:`modcode.pipeline.Link`: the scalar-noise LLR path ignores the ISI.  The
equalised, estimation-aware path lives in :mod:`modcode.pipeline.ofdm`.
"""

from __future__ import annotations

from collections.abc import Sequence
from numbers import Integral

import numpy as np
from numpy.typing import NDArray

from modcode.channel.awgn import AWGNChannel
from modcode.channel.base import Channel


class TapProcess:
    """Unit-average-power Rayleigh taps with an AR(1) Doppler correlation.

    ``taps(n)`` returns the next ``n`` consecutive tap vectors, each of shape
    ``(L,)``.  Each block has ``E|h_l|^2 = powers_l``; consecutive blocks have
    correlation ``rho`` from the Jakes one-step coefficient, so ``rho = 1`` is a
    static channel and ``rho = 0`` gives independent blocks.
    """

    def __init__(
        self,
        powers: Sequence[float],
        *,
        doppler_hz: float = 0.0,
        symbol_rate_hz: float = 1.0,
        block_symbols: int = 1,
        sinusoids: int = 16,
        seed: int | None = 0,
    ) -> None:
        taps = np.asarray(powers, dtype=np.float64)
        if taps.ndim != 1 or taps.size == 0:
            raise ValueError("powers must be a non-empty one-dimensional sequence")
        if not np.all(np.isfinite(taps)) or np.any(taps <= 0.0):
            raise ValueError("powers must be finite and strictly positive")
        if not np.isfinite(doppler_hz) or doppler_hz < 0.0:
            raise ValueError("doppler_hz must be finite and non-negative")
        if not np.isfinite(symbol_rate_hz) or symbol_rate_hz <= 0.0:
            raise ValueError("symbol_rate_hz must be finite and positive")
        if (
            isinstance(block_symbols, bool)
            or not isinstance(block_symbols, Integral)
            or block_symbols < 1
        ):
            raise ValueError("block_symbols must be a positive integer")
        if (
            isinstance(sinusoids, bool)
            or not isinstance(sinusoids, Integral)
            or sinusoids < 1
        ):
            raise ValueError("sinusoids must be a positive integer")
        self.powers = taps / taps.sum()
        self.taps_count = self.powers.size
        self.doppler_hz = float(doppler_hz)
        self.symbol_rate_hz = float(symbol_rate_hz)
        self.block_symbols = int(block_symbols)
        self.sinusoids = int(sinusoids)
        self.rng = np.random.default_rng(seed)
        self._state: NDArray[np.complex128] | None = None
        self._rho = self._jakes_correlation(
            self.normalized_doppler * self.block_symbols
        )

    def _jakes_correlation(self, argument: float) -> float:
        if argument <= 0.0:
            return 1.0
        angles = 2.0 * np.pi * np.arange(1, self.sinusoids + 1) / (
            self.sinusoids + 1
        )
        values = np.cos(2.0 * np.pi * argument * np.cos(angles))
        return float(np.clip(np.mean(values), 0.0, 1.0))

    @property
    def normalized_doppler(self) -> float:
        """Maximum Doppler in cycles per symbol."""
        return self.doppler_hz / self.symbol_rate_hz

    @property
    def correlation(self) -> float:
        """One-step tap correlation, from the Jakes coefficient."""
        return self._rho

    def reset(self, seed: int | None = None) -> None:
        """Forget the current state and optionally reseed the generator."""
        self._state = None
        if seed is not None:
            self.rng = np.random.default_rng(seed)

    def taps(self, blocks: int = 1) -> NDArray[np.complex128]:
        """Return the next ``blocks`` tap vectors, shape ``(blocks, L)``."""
        if (
            isinstance(blocks, bool)
            or not isinstance(blocks, Integral)
            or blocks < 1
        ):
            raise ValueError("blocks must be a positive integer")
        innovations = np.sqrt(self.powers / 2.0) * (
            self.rng.standard_normal((int(blocks), self.taps_count))
            + 1j * self.rng.standard_normal((int(blocks), self.taps_count))
        )
        values = np.empty((int(blocks), self.taps_count), dtype=np.complex128)
        state = self._state
        if state is None:
            state = innovations[0].copy()
            values[0] = state
            start = 1
        else:
            start = 0
        scale = np.sqrt(max(1.0 - self._rho**2, 0.0))
        for index in range(start, int(blocks)):
            state = self._rho * state + scale * innovations[index]
            values[index] = state
        self._state = state
        return values


class MultipathChannel(Channel):
    """Symbol-spaced multipath with optional Doppler; no equalisation.

    ``transmit`` applies the tap-delay line to the symbol stream and adds
    circular complex AWGN.  The average channel power is one, so
    ``noise_variance`` matches :class:`AWGNChannel`.  ``last_taps`` stores the
    gains used by the most recent call, with shape ``(samples, taps)``.
    """

    def __init__(
        self,
        powers: Sequence[float] = (0.7, 0.3),
        delays: Sequence[int] = (0, 1),
        *,
        doppler_hz: float = 0.0,
        symbol_rate_hz: float = 1.0,
        sinusoids: int = 16,
        seed: int | None = 0,
    ) -> None:
        delay_values = np.asarray(delays)
        if delay_values.ndim != 1 or delay_values.size == 0:
            raise ValueError("delays must be a non-empty one-dimensional sequence")
        if not np.all(delay_values >= 0):
            raise ValueError("delays must be non-negative")
        process = TapProcess(
            powers,
            doppler_hz=doppler_hz,
            symbol_rate_hz=symbol_rate_hz,
            block_symbols=1,
            sinusoids=sinusoids,
            seed=seed,
        )
        if delay_values.size != process.taps_count:
            raise ValueError("delays and powers must have the same length")
        self.delays = tuple(int(delay) for delay in delay_values)
        self.process = process
        self.last_taps: NDArray[np.complex128] | None = None

    @property
    def powers(self) -> NDArray[np.float64]:
        return self.process.powers

    @staticmethod
    def noise_variance(
        ebn0_db: float, bits_per_symbol: int, code_rate: float = 1.0
    ) -> float:
        return AWGNChannel.noise_variance(ebn0_db, bits_per_symbol, code_rate)

    def transmit(
        self,
        symbols: NDArray[np.complexfloating],
        noise_variance: float,
        rng: np.random.Generator | None = None,
    ) -> NDArray[np.complex128]:
        values = np.asarray(symbols, dtype=np.complex128)
        if values.ndim != 1:
            raise ValueError("symbols must be a one-dimensional array")
        if not np.all(np.isfinite(values)):
            raise ValueError("symbols must contain only finite values")
        variance = float(noise_variance)
        if not np.isfinite(variance) or variance < 0.0:
            raise ValueError("noise_variance must be finite and non-negative")
        generator = rng if rng is not None else np.random.default_rng()
        taps = self.process.taps(values.size)
        received = np.zeros(values.size, dtype=np.complex128)
        for index, delay in enumerate(self.delays):
            if delay >= values.size:
                continue
            received[delay:] += taps[delay:, index] * values[: values.size - delay]
        received += np.sqrt(variance) * (
            generator.standard_normal(values.size)
            + 1j * generator.standard_normal(values.size)
        )
        self.last_taps = taps
        return received


def subcarrier_gains(
    taps: NDArray[np.complexfloating],
    delays: Sequence[int],
    n_fft: int,
) -> NDArray[np.complex128]:
    """Return the DFT of the tap-delay line at the ``n_fft`` subcarriers.

    ``taps`` may carry leading dimensions; the result has shape
    ``taps.shape[:-1] + (n_fft,)`` and ``H[m] = sum_l h_l e^{-j 2 pi m d_l /
    n_fft}``.
    """
    gains = np.asarray(taps, dtype=np.complex128)
    if gains.ndim == 0:
        raise ValueError("taps must have at least one dimension")
    if isinstance(n_fft, bool) or not isinstance(n_fft, Integral) or n_fft < 1:
        raise ValueError("n_fft must be a positive integer")
    delay_array = np.asarray(delays)
    if delay_array.ndim != 1 or delay_array.size != gains.shape[-1]:
        raise ValueError("delays must match the last taps dimension")
    if not np.all(delay_array >= 0):
        raise ValueError("delays must be non-negative")
    subcarriers = np.arange(int(n_fft))
    phase = np.exp(
        -2j * np.pi * subcarriers[None, :] * delay_array[:, None] / int(n_fft)
    )
    return gains @ phase


__all__ = ["MultipathChannel", "TapProcess", "subcarrier_gains"]
