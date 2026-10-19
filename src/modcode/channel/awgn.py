# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Additive white Gaussian noise channel."""

import numpy as np
from numpy.typing import NDArray

from modcode.channel.base import Channel


class AWGNChannel(Channel):
    """Complex AWGN with Eb/N0 normalized per information bit."""

    @staticmethod
    def noise_variance(
        ebn0_db: float, bits_per_symbol: int, code_rate: float = 1.0
    ) -> float:
        """Return noise variance per real component for unit symbol energy."""
        if not np.isfinite(ebn0_db):
            raise ValueError("ebn0_db must be finite")
        if bits_per_symbol <= 0:
            raise ValueError("bits_per_symbol must be positive")
        if not 0.0 < code_rate <= 1.0:
            raise ValueError("code_rate must be in the interval (0, 1]")
        ebn0 = 10.0 ** (ebn0_db / 10.0)
        return 1.0 / (2.0 * bits_per_symbol * code_rate * ebn0)

    def transmit(
        self,
        symbols: NDArray[np.complexfloating],
        ebn0_db: float,
        bits_per_symbol: int,
        code_rate: float = 1.0,
        rng: np.random.Generator | None = None,
    ) -> NDArray[np.complex128]:
        values = np.asarray(symbols, dtype=np.complex128)
        if values.ndim != 1:
            raise ValueError("symbols must be a one-dimensional array")
        if not np.all(np.isfinite(values)):
            raise ValueError("symbols must contain only finite values")
        variance = self.noise_variance(ebn0_db, bits_per_symbol, code_rate)
        generator = rng if rng is not None else np.random.default_rng()
        noise = np.sqrt(variance) * (
            generator.standard_normal(values.size)
            + 1j * generator.standard_normal(values.size)
        )
        return values + noise
