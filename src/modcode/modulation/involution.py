# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Hypothesis tests for the tau (Q-axis reflection) ambiguity.

At ``phi = pi`` the I/Q phase-imbalance channel maps ``x`` to ``conj(x)``.
The two hypotheses are therefore ``H0: y = x + n`` and
``H1: y = conj(x) + n``.  This module supplies the likelihood ratio, a
block-level decision, and a pilot-based decision.

The blind likelihood ratio is deliberately allowed to be identically zero.
For a constellation and data prior invariant under conjugation, the marginal
likelihood of the received vector is the same under both hypotheses for every
observation, so no blind threshold can do better than chance.  The experiment
in ``examples/phase_imbalance_lab.py`` demonstrates that limit; a pilot or an
asymmetric code constraint is what breaks it.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from modcode.modulation.base import Modulation


def _constellation_points(modulation: Modulation) -> NDArray[np.complex128]:
    points = np.asarray(
        list(modulation.constellation.values()), dtype=np.complex128
    )
    if points.size == 0:
        raise ValueError("modulation constellation must not be empty")
    if not np.all(np.isfinite(points)):
        raise ValueError("constellation points must be finite")
    return points


def _logsumexp(values: NDArray[np.float64]) -> NDArray[np.float64]:
    peak = np.max(values, axis=1, keepdims=True)
    return np.squeeze(peak, axis=1) + np.log(
        np.sum(np.exp(values - peak), axis=1)
    )


def tau_log_likelihood_ratio(
    modulation: Modulation,
    symbols: NDArray[np.complexfloating],
    *,
    noise_variance: float = 1.0,
) -> NDArray[np.float64]:
    """Return ``log p(y | H0) - log p(y | H1)`` per received symbol.

    ``H0`` is the ideal channel and ``H1`` is the ``phi = pi`` reflection.
    The constellation prior is uniform.  The noise is circular complex AWGN
    with variance ``noise_variance`` per real component.  Both hypotheses are
    unitary, so the reflection does not colour the noise.
    """
    values = np.asarray(symbols, dtype=np.complex128)
    if values.ndim != 1:
        raise ValueError("symbols must be a one-dimensional array")
    if not np.all(np.isfinite(values)):
        raise ValueError("symbols must contain only finite values")
    variance = float(noise_variance)
    if not np.isfinite(variance) or variance <= 0.0:
        raise ValueError("noise_variance must be finite and positive")
    if values.size == 0:
        return np.empty(0, dtype=np.float64)

    points = _constellation_points(modulation)

    def log_likelihood(observations: NDArray[np.complex128]) -> NDArray[np.float64]:
        distances = np.abs(
            observations[:, None] - points[None, :]
        ) ** 2
        return _logsumexp(-distances / (2.0 * variance))

    return log_likelihood(values) - log_likelihood(np.conj(values))


def block_tau_hypothesis(
    modulation: Modulation,
    symbols: NDArray[np.complexfloating],
    *,
    noise_variance: float = 1.0,
    tie_tolerance: float = 1e-12,
) -> int:
    """Choose ``0`` or ``1`` for a block with one constant phase hypothesis.

    The per-symbol log-likelihood ratios are summed because the hypothesis is
    constant over the block.  An exact tie returns ``0``.  A tie is not a
    defect: for an invariant constellation and prior it is the only possible
    outcome, which is why a blind detector cannot resolve the ambiguity.
    ``tie_tolerance`` is a per-symbol numerical threshold; a total whose
    magnitude is at or below ``tie_tolerance * len(symbols)`` is treated as a
    tie rather than as evidence.
    """
    if not np.isfinite(tie_tolerance) or tie_tolerance < 0.0:
        raise ValueError("tie_tolerance must be finite and non-negative")
    llr = tau_log_likelihood_ratio(
        modulation, symbols, noise_variance=noise_variance
    )
    if llr.size == 0:
        raise ValueError("at least one symbol is required")
    total = float(np.sum(llr))
    if abs(total) <= tie_tolerance * llr.size:
        return 0
    return int(total < 0.0)


def pilot_tau_hypothesis(
    received: NDArray[np.complexfloating],
    transmitted: NDArray[np.complexfloating],
) -> int:
    """Choose the tau hypothesis from known pilot symbols.

    Returns ``0`` for the ideal channel and ``1`` for the ``phi = pi``
    reflection.  A tie returns ``0``; a pilot lying on the I axis is fixed by
    the reflection and therefore carries no information about the hypothesis.
    """
    rx = np.asarray(received, dtype=np.complex128)
    tx = np.asarray(transmitted, dtype=np.complex128)
    if rx.ndim != 1 or tx.ndim != 1:
        raise ValueError("received and transmitted pilots must be 1-D")
    if rx.shape != tx.shape:
        raise ValueError("received and transmitted pilots must have the same shape")
    if rx.size == 0:
        raise ValueError("at least one pilot symbol is required")
    if not np.all(np.isfinite(rx)) or not np.all(np.isfinite(tx)):
        raise ValueError("pilot symbols must be finite")
    ideal_distance = float(np.sum(np.abs(rx - tx) ** 2))
    reflected_distance = float(np.sum(np.abs(rx - np.conj(tx)) ** 2))
    return 0 if ideal_distance <= reflected_distance else 1


def correct_tau(
    symbols: NDArray[np.complexfloating], hypothesis: int
) -> NDArray[np.complex128]:
    """Apply the inverse tau operation for a known hypothesis."""
    values = np.asarray(symbols, dtype=np.complex128)
    if values.ndim != 1:
        raise ValueError("symbols must be a one-dimensional array")
    if not np.all(np.isfinite(values)):
        raise ValueError("symbols must contain only finite values")
    if isinstance(hypothesis, bool) or hypothesis not in (0, 1):
        raise ValueError("hypothesis must be 0 or 1")
    if hypothesis == 0:
        return values.copy()
    return np.conj(values)


__all__ = [
    "block_tau_hypothesis",
    "correct_tau",
    "pilot_tau_hypothesis",
    "tau_log_likelihood_ratio",
]
