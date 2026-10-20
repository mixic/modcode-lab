# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Toy quasi-cyclic code-based cryptography demonstrations."""

from modcode.pqc.analysis import (
    decoder_success_rates,
    dual_weight_distribution,
    syndrome_search_space_bits,
)
from modcode.pqc.qc_mdpc import (
    QCMDPC,
    DecodeResult,
    circulant,
    is_circulant,
)

__all__ = [
    "QCMDPC",
    "DecodeResult",
    "circulant",
    "decoder_success_rates",
    "dual_weight_distribution",
    "is_circulant",
    "syndrome_search_space_bits",
]
