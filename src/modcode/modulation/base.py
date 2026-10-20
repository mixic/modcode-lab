# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Shared interface and input checks for digital modulation schemes."""

from abc import ABC, abstractmethod

import numpy as np
from numpy.typing import NDArray

BitArray = NDArray[np.uint8]
SymbolArray = NDArray[np.complex128]


class Modulation(ABC):
    """Map binary data to constellation symbols and back."""

    @property
    @abstractmethod
    def bits_per_symbol(self) -> int:
        """Number of binary digits represented by one symbol."""

    @property
    @abstractmethod
    def constellation(self) -> dict[tuple[int, ...], complex]:
        """Return the mapping from bit labels to complex symbols."""

    @abstractmethod
    def modulate(self, bits: NDArray[np.integer]) -> SymbolArray:
        """Map a one-dimensional bit sequence to complex symbols."""

    @abstractmethod
    def demodulate(self, symbols: NDArray[np.complexfloating]) -> BitArray:
        """Make hard decisions and return the corresponding bits."""

    def demodulate_phase_imbalance(
        self, symbols: NDArray[np.complexfloating], phi: float
    ) -> BitArray:
        """Demodulate after inverting a known I/Q phase mismatch ``phi``.

        The forward channel is ``y = I + j Q exp(j phi)``.  Its inverse is

            I = Re(y) + Im(y) tan(phi),
            Q = Im(y) / cos(phi).

        The inverse exists except at ``phi = pi/2 (mod pi)``, where the
        forward map is rank one and the Q component is lost.  This is a
        known-phase correction; it does not estimate ``phi`` and does not
        resolve the blind ``phi = 0`` versus ``phi = pi`` ambiguity.
        """
        values = self._validated_symbols(symbols)
        angle = float(phi)
        if not np.isfinite(angle):
            raise ValueError("phi must be finite")
        cosine = float(np.cos(angle))
        if abs(cosine) < 1e-12:
            raise ValueError(
                "phase imbalance is singular at phi = pi/2 (mod pi)"
            )
        corrected = (values.real + values.imag * np.tan(angle)) + 1j * (
            values.imag / cosine
        )
        return self.demodulate(corrected)

    def demodulate_tau_aware(
        self, symbols: NDArray[np.complexfloating], *, hypothesis: int
    ) -> BitArray:
        """Demodulate under one of the two tau hypotheses.

        ``hypothesis=0`` assumes the ideal channel; ``hypothesis=1`` assumes
        the ``phi = pi`` Q-axis reflection, so the received symbols are
        conjugated before slicing.  The hypothesis is an input, not an
        estimate: for a constellation and data prior invariant under the
        reflection, the two hypotheses have identical likelihoods and no
        blind detector can choose between them.
        """
        values = self._validated_symbols(symbols)
        if isinstance(hypothesis, bool) or hypothesis not in (0, 1):
            raise ValueError("hypothesis must be 0 or 1")
        if hypothesis == 0:
            return self.demodulate(values)
        return self.demodulate(np.conj(values))

    def _validated_bits(self, bits: NDArray[np.integer]) -> BitArray:
        values = np.asarray(bits)
        if values.ndim != 1:
            raise ValueError("bits must be a one-dimensional array")
        if not np.all((values == 0) | (values == 1)):
            raise ValueError("bits must contain only 0 and 1")
        return values.astype(np.uint8, copy=False)

    def _validated_symbols(
        self, symbols: NDArray[np.complexfloating]
    ) -> NDArray[np.complex128]:
        values = np.asarray(symbols, dtype=np.complex128)
        if values.ndim != 1:
            raise ValueError("symbols must be a one-dimensional array")
        if not np.all(np.isfinite(values)):
            raise ValueError("symbols must contain only finite values")
        return values
