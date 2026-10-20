# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Deterministic bit interleavers for the link pipeline."""

from __future__ import annotations

from numbers import Integral

import numpy as np
from numpy.typing import ArrayLike, NDArray


class PermutationInterleaver:
    """A deterministic permutation of an arbitrary-length bit block.

    The permutation is regenerated from ``seed`` and the block length, so
    :meth:`interleave` and :meth:`deinterleave` agree without storing state or
    transmitting the permutation.  The interleaver is applied to encoded bits
    before modulation and inverted before FEC decoding.

    A random permutation is enough for the lab's purpose: it spreads the bits
    of one codeword across different symbols, so a QAM symbol error does not
    dump two adjacent bits into the same Hamming or BCH codeword.  It is not a
    standards-compliant convolutional or block interleaver.
    """

    def __init__(self, seed: int = 0) -> None:
        if isinstance(seed, bool) or not isinstance(seed, Integral):
            raise ValueError("seed must be an integer")
        self.seed = int(seed)

    def _permutation(self, length: int) -> NDArray[np.intp]:
        if length <= 0:
            raise ValueError("interleaver block length must be positive")
        return np.random.default_rng(self.seed).permutation(length)

    def interleave(self, values: ArrayLike) -> NDArray:
        """Permute a one-dimensional block."""
        array = np.asarray(values)
        if array.ndim != 1:
            raise ValueError("interleaver input must be one-dimensional")
        if array.size == 0:
            raise ValueError("interleaver input must be nonempty")
        return array[self._permutation(array.size)]

    def deinterleave(self, values: ArrayLike) -> NDArray:
        """Invert :meth:`interleave` for the same block length."""
        array = np.asarray(values)
        if array.ndim != 1:
            raise ValueError("interleaver input must be one-dimensional")
        if array.size == 0:
            raise ValueError("interleaver input must be nonempty")
        inverse = np.argsort(self._permutation(array.size), kind="stable")
        return array[inverse]


__all__ = ["PermutationInterleaver"]
