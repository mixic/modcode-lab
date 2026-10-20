# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Security-oriented measurements for the toy QC-MDPC demonstration.

These functions are diagnostic experiments, not a security proof. They make
three facts visible:

* the bit-flipping decoder has an error-weight threshold;
* the syndrome-decoding search space grows exponentially with error weight;
* the dual of a QC-MDPC code contains low-weight words, which is exactly what
  key-recovery attacks search for when parameters are too small.
"""

from __future__ import annotations

from math import lgamma, log

import numpy as np
from numpy.typing import ArrayLike, NDArray

from modcode.pqc.qc_mdpc import QCMDPC


def decoder_success_rates(
    code: QCMDPC,
    error_weights: ArrayLike,
    *,
    trials: int = 20,
    seed: int = 0,
) -> tuple[NDArray[np.int64], NDArray[np.float64]]:
    """Estimate bit-flipping success for each error weight.

    The toy decoder's exhaustive fallback is deliberately not used here, so
    the curve shows the iterative decoder's threshold rather than a search
    that would not scale.
    """
    weights = np.asarray(error_weights, dtype=np.int64)
    if weights.ndim != 1 or weights.size == 0:
        raise ValueError("error_weights must be a nonempty one-dimensional array")
    if np.any(weights < 0) or np.any(weights > code.code_length):
        raise ValueError("error weights must be in [0, code_length]")
    if isinstance(trials, bool) or not isinstance(trials, int) or trials < 1:
        raise ValueError("trials must be a positive integer")

    rng = np.random.default_rng(seed)
    rates = np.empty(weights.size, dtype=np.float64)
    for index, weight in enumerate(weights):
        successes = 0
        for _ in range(trials):
            message = rng.integers(
                0, 2, size=code.dimension, dtype=np.uint8
            )
            ciphertext, _ = code.encrypt(
                message,
                int(weight),
                seed=int(rng.integers(0, 2**31 - 1)),
            )
            _, _, _, success = code.bit_flip_decode(ciphertext)
            successes += int(success)
        rates[index] = successes / trials
    return weights, rates


def syndrome_search_space_bits(
    code_length: int, error_weight: int
) -> float:
    """Return ``log2(C(code_length, error_weight))``.

    This is the number of error patterns of that weight; an unstructured
    syndrome decoder must search a space of this size. For a real code the
    decoding problem is not solved by brute force, but the exponential growth
    is the reason the parameters cannot simply be scaled down.
    """
    if isinstance(code_length, bool) or not isinstance(code_length, int):
        raise ValueError("code_length must be an integer")
    if code_length < 1:
        raise ValueError("code_length must be positive")
    if isinstance(error_weight, bool) or not isinstance(error_weight, int):
        raise ValueError("error_weight must be an integer")
    if not 0 <= error_weight <= code_length:
        raise ValueError("error_weight must be in [0, code_length]")
    if error_weight == 0:
        return 0.0
    value = (
        lgamma(code_length + 1)
        - lgamma(error_weight + 1)
        - lgamma(code_length - error_weight + 1)
    ) / log(2.0)
    return float(value)


def dual_weight_distribution(
    code: QCMDPC, *, max_block_length: int = 18
) -> NDArray[np.int64]:
    """Enumerate the weight distribution of the dual code.

    The dual code is the row space of ``H``. Enumerating all ``2**n`` linear
    combinations is only feasible for the tiny blocks used here, but it shows
    why the sparse secret rows are special: the dual contains words of weight
    ``2 * row_weight`` while a random word has weight near ``n``.

    This is an inspection tool for toy parameters, not a scalable attack.
    """
    if code.block_length > max_block_length:
        raise ValueError(
            f"block_length must be at most {max_block_length} for enumeration"
        )
    rows = code.parity_check
    size = code.block_length
    counts = np.zeros(code.code_length + 1, dtype=np.int64)
    current = np.zeros(code.code_length, dtype=np.uint8)
    counts[0] = 1
    # Gray-code enumeration: exactly one row is toggled at each step.
    for index in range(1, 1 << size):
        bit = (index & -index).bit_length() - 1
        current ^= rows[bit]
        counts[int(current.sum())] += 1
    return counts


__all__ = [
    "decoder_success_rates",
    "dual_weight_distribution",
    "syndrome_search_space_bits",
]
