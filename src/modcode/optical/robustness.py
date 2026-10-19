# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""BER experiments and candidate level optimization for optical IM/DD."""

from collections.abc import Mapping

import numpy as np
from numpy.typing import NDArray

from modcode.optical.channel import IMDDChannel
from modcode.optical.modulation import IntensityPAM, NonUniformPAM


def optimize_pam_levels(
    noise_std: float,
    amplitude_noise_std: float = 0.0,
    grid_points: int = 2001,
) -> NonUniformPAM:
    """Search for 4-PAM levels maximizing worst adjacent noise margin.

    The search fixes the outer levels at 0 and 2 to preserve unit mean power
    and peak intensity. Its objective is the minimum adjacent separation
    divided by the sum of the two levels' modeled noise standard deviations.
    This is an experimental design criterion, not a claim of global
    optimality or universal BER improvement.
    """
    if not np.isfinite(noise_std) or noise_std <= 0:
        raise ValueError("noise_std must be finite and positive")
    if not np.isfinite(amplitude_noise_std) or amplitude_noise_std < 0:
        raise ValueError("amplitude_noise_std must be finite and nonnegative")
    if isinstance(grid_points, bool) or not isinstance(grid_points, int):
        raise ValueError("grid_points must be an integer")
    if grid_points < 101:
        raise ValueError("grid_points must be at least 101")

    lower = np.linspace(0.05, 1.0, grid_points)
    upper = 2.0 - lower
    candidates = np.column_stack(
        (np.zeros_like(lower), lower, upper, np.full_like(lower, 2.0))
    )
    spacings = np.diff(candidates, axis=1)
    deviations = np.sqrt(noise_std**2 + (amplitude_noise_std * candidates) ** 2)
    margins = spacings / (deviations[:, :-1] + deviations[:, 1:])
    scores = margins.min(axis=1)
    selected = candidates[int(np.argmax(scores))]
    return NonUniformPAM(selected)


def ber_sweep_imdd(
    modulations: Mapping[str, IntensityPAM],
    channels: Mapping[str, IMDDChannel],
    n_bits: int,
    seed: int | None = None,
) -> dict[str, dict[str, NDArray[np.float64]]]:
    """Measure hard-decision BER for each modulation/channel combination."""
    if not modulations:
        raise ValueError("at least one modulation is required")
    if not channels:
        raise ValueError("at least one channel is required")
    if n_bits <= 0:
        raise ValueError("n_bits must be positive")

    generator = np.random.default_rng(seed)
    results: dict[str, dict[str, NDArray[np.float64]]] = {}
    for modulation_name, modulation in modulations.items():
        bit_count = (
            (n_bits + modulation.bits_per_symbol - 1)
            // modulation.bits_per_symbol
            * modulation.bits_per_symbol
        )
        bits = generator.integers(0, 2, bit_count, dtype=np.uint8)
        symbols = modulation.modulate(bits)
        channel_results: dict[str, NDArray[np.float64]] = {}
        for channel_name, channel in channels.items():
            received = channel.transmit(symbols, rng=generator)
            decoded = modulation.demodulate(received)
            errors = decoded[:n_bits] != bits[:n_bits]
            channel_results[channel_name] = np.asarray(
                [np.mean(errors)], dtype=np.float64
            )
        results[modulation_name] = channel_results
    return results
