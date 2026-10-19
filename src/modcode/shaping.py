# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Small utilities for probabilistic constellation and sequence experiments."""

import numpy as np
from numpy.typing import NDArray


def maxwell_boltzmann_probabilities(
    constellation: NDArray[np.complexfloating],
    shaping_parameter: float,
) -> NDArray[np.float64]:
    """Return normalized energy-biased probabilities for complex symbols."""
    points = np.asarray(constellation, dtype=np.complex128)
    if points.ndim != 1 or points.size < 2:
        raise ValueError("constellation must be a one-dimensional array of points")
    if not np.all(np.isfinite(points)):
        raise ValueError("constellation points must be finite")
    if not np.isfinite(shaping_parameter) or shaping_parameter < 0:
        raise ValueError("shaping_parameter must be finite and nonnegative")
    log_weights = -shaping_parameter * np.abs(points) ** 2
    log_weights -= np.max(log_weights)
    weights = np.exp(log_weights)
    return weights / weights.sum()


def sample_shaped_symbols(
    constellation: NDArray[np.complexfloating],
    probabilities: NDArray[np.floating],
    count: int,
    rng: np.random.Generator | None = None,
) -> NDArray[np.complex128]:
    """Sample symbols according to an explicit probability mass function."""
    points = np.asarray(constellation, dtype=np.complex128)
    probs = np.asarray(probabilities, dtype=np.float64)
    if points.ndim != 1 or points.size < 2 or not np.all(np.isfinite(points)):
        raise ValueError("constellation must contain at least two finite points")
    if probs.shape != points.shape or not np.all(np.isfinite(probs)):
        raise ValueError("probabilities must be finite and match the constellation")
    if np.any(probs < 0) or not np.isclose(probs.sum(), 1.0):
        raise ValueError("probabilities must be nonnegative and sum to one")
    if isinstance(count, bool) or not isinstance(count, int) or count < 0:
        raise ValueError("count must be a nonnegative integer")
    generator = rng if rng is not None else np.random.default_rng()
    indexes = generator.choice(points.size, size=count, p=probs)
    return points[indexes]


def periodic_autocorrelation(
    sequence: NDArray[np.complexfloating],
) -> NDArray[np.complex128]:
    """Return the normalized periodic autocorrelation for every cyclic lag."""
    values = np.asarray(sequence, dtype=np.complex128)
    if values.ndim != 1 or values.size == 0:
        raise ValueError("sequence must be a nonempty one-dimensional array")
    if not np.all(np.isfinite(values)):
        raise ValueError("sequence values must be finite")
    return np.asarray(
        [
            np.vdot(values, np.roll(values, lag)) / values.size
            for lag in range(values.size)
        ],
        dtype=np.complex128,
    )
