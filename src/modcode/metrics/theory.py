# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Closed-form uncoded BER references for AWGN."""

from __future__ import annotations

import math

import numpy as np
from numpy.typing import ArrayLike, NDArray

from modcode.modulation.base import Modulation
from modcode.modulation.bpsk import BPSK
from modcode.modulation.qam import QAM
from modcode.modulation.qpsk import QPSK


def q_function(values: ArrayLike) -> NDArray[np.float64]:
    """Return the Gaussian tail probability ``Q(x)`` elementwise."""
    array = np.asarray(values, dtype=np.float64)
    if not np.all(np.isfinite(array)):
        raise ValueError("Q-function arguments must be finite")
    flat = np.ravel(array)
    result = np.asarray(
        [0.5 * math.erfc(float(value) / math.sqrt(2.0)) for value in flat],
        dtype=np.float64,
    )
    return result.reshape(array.shape)


def theoretical_ber(
    modulation: Modulation, ebn0_db: ArrayLike
) -> NDArray[np.float64]:
    """Return the uncoded AWGN BER reference for a supported modulation.

    BPSK and Gray-coded QPSK have the same bit error probability,
    ``Q(sqrt(2 Eb/N0))``.  Square Gray-coded QAM uses the standard
    nearest-neighbour approximation.  Other constellations are deliberately
    rejected rather than given a misleading closed form.
    """
    points = np.asarray(ebn0_db, dtype=np.float64)
    if points.ndim != 1 or points.size == 0:
        raise ValueError("ebn0_db must be a non-empty one-dimensional sequence")
    if not np.all(np.isfinite(points)):
        raise ValueError("ebn0_db values must be finite")
    ebn0 = 10.0 ** (points / 10.0)

    if isinstance(modulation, (BPSK, QPSK)):
        return q_function(np.sqrt(2.0 * ebn0))

    if isinstance(modulation, QAM):
        order = modulation.order
        bits = math.log2(order)
        side = math.sqrt(order)
        argument = np.sqrt(3.0 * bits * ebn0 / (order - 1.0))
        return (
            4.0
            / bits
            * (1.0 - 1.0 / side)
            * q_function(argument)
        )

    raise ValueError(
        f"no theoretical BER reference for {type(modulation).__name__}"
    )


def theoretical_rayleigh_ber(
    modulation: Modulation, ebn0_db: ArrayLike
) -> NDArray[np.float64]:
    """Return the uncoded flat-Rayleigh BER with known CSI.

    BPSK and Gray-coded QPSK have the same per-bit expression,
    ``0.5 * (1 - sqrt(gamma / (1 + gamma)))`` with ``gamma = Eb/N0``.
    """
    points = np.asarray(ebn0_db, dtype=np.float64)
    if points.ndim != 1 or points.size == 0:
        raise ValueError("ebn0_db must be a non-empty one-dimensional sequence")
    if not np.all(np.isfinite(points)):
        raise ValueError("ebn0_db values must be finite")
    if not isinstance(modulation, (BPSK, QPSK)):
        raise ValueError(
            f"no Rayleigh reference for {type(modulation).__name__}"
        )
    ebn0 = 10.0 ** (points / 10.0)
    return 0.5 * (1.0 - np.sqrt(ebn0 / (1.0 + ebn0)))


__all__ = ["q_function", "theoretical_ber", "theoretical_rayleigh_ber"]
