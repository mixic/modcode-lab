# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Shared interface for block forward-error-correction codes."""

from abc import ABC, abstractmethod

import numpy as np
from numpy.typing import NDArray

BitArray = NDArray[np.uint8]


class FEC(ABC):
    """Encode data bits and decode received codeword bits."""

    @property
    @abstractmethod
    def data_bits_per_block(self) -> int:
        """Number of uncoded input bits in each codeword."""

    @property
    @abstractmethod
    def encoded_bits_per_block(self) -> int:
        """Number of encoded bits in each codeword."""

    @property
    @abstractmethod
    def rate(self) -> float:
        """Number of data bits per encoded bit."""

    @abstractmethod
    def encode(self, bits: NDArray[np.integer]) -> BitArray:
        """Encode a one-dimensional data bit sequence."""

    @abstractmethod
    def decode(self, bits: NDArray[np.integer]) -> BitArray:
        """Decode a one-dimensional sequence of encoded bits."""

    def decode_soft(self, llrs: NDArray[np.floating]) -> BitArray:
        """Decode from bit LLRs, positive for bit 0.

        The default is unavailable so that a hard-decision-only code fails
        loudly if a caller asks for soft decoding.  ``run_link`` falls back to
        hard decisions when this method is not overridden.
        """
        raise NotImplementedError(
            f"{type(self).__name__} does not implement soft decoding"
        )

    @staticmethod
    def validated_bits(bits: NDArray[np.integer], name: str = "bits") -> BitArray:
        values = np.asarray(bits)
        if values.ndim != 1:
            raise ValueError(f"{name} must be a one-dimensional array")
        if not np.all((values == 0) | (values == 1)):
            raise ValueError(f"{name} must contain only 0 and 1")
        return values.astype(np.uint8, copy=False)
