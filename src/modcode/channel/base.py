# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Shared interface for noisy communication channels."""

from abc import ABC, abstractmethod

import numpy as np
from numpy.typing import NDArray


class Channel(ABC):
    """Transform transmitted symbols according to channel conditions."""

    @abstractmethod
    def transmit(
        self,
        symbols: NDArray[np.complexfloating],
        ebn0_db: float,
        bits_per_symbol: int,
        code_rate: float = 1.0,
        rng: np.random.Generator | None = None,
    ) -> NDArray[np.complex128]:
        """Return the received complex symbols."""
