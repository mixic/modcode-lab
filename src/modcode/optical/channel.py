# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Simple symbol-rate impairment models for optical IM/DD experiments."""

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray


@dataclass(frozen=True)
class IMDDChannel:
    """Apply selected memoryless and symbol-rate ISI impairments.

    ``bandwidth`` is a normalized three-tap ISI proxy, not a physical optical
    filter model. ``nonlinearity`` uses a memoryless compression curve.
    """

    noise_std: float = 0.0
    amplitude_noise_std: float = 0.0
    nonlinearity: float = 0.0
    bandwidth: float = 0.0

    def __post_init__(self) -> None:
        parameters = {
            "noise_std": self.noise_std,
            "amplitude_noise_std": self.amplitude_noise_std,
            "nonlinearity": self.nonlinearity,
            "bandwidth": self.bandwidth,
        }
        for name, value in parameters.items():
            if not np.isfinite(value) or value < 0:
                raise ValueError(f"{name} must be finite and nonnegative")
        if self.bandwidth > 0.25:
            raise ValueError("bandwidth proxy must be in the interval [0, 0.25]")

    def transmit(
        self,
        symbols: NDArray[np.floating],
        rng: np.random.Generator | None = None,
    ) -> NDArray[np.float64]:
        """Apply amplitude fluctuation, compression, bandwidth proxy, and noise."""
        values = np.asarray(symbols, dtype=np.float64)
        if values.ndim != 1:
            raise ValueError("symbols must be a one-dimensional array")
        if not np.all(np.isfinite(values)) or np.any(values < 0):
            raise ValueError("intensity symbols must be finite and nonnegative")

        generator = rng if rng is not None else np.random.default_rng()
        received = values.copy()
        if self.amplitude_noise_std:
            gain = 1.0 + generator.normal(0.0, self.amplitude_noise_std, values.size)
            received *= np.maximum(gain, 0.0)
        if self.nonlinearity:
            received = received / (1.0 + self.nonlinearity * received)
        if self.bandwidth and received.size > 2:
            edge = self.bandwidth
            received = np.convolve(
                received, [edge, 1.0 - 2.0 * edge, edge], mode="same"
            )
        if self.noise_std:
            received += generator.normal(0.0, self.noise_std, values.size)
        return received
