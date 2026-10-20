# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Forward error-correction codes."""

from modcode.fec.base import FEC
from modcode.fec.bch import BCH15
from modcode.fec.convolutional import ConvolutionalCode
from modcode.fec.crc import CRC
from modcode.fec.hamming import Hamming74
from modcode.fec.ldpc import LDPC36
from modcode.fec.polar import PolarAWGN, PolarBEC
from modcode.fec.reed_solomon import ReedSolomon
from modcode.fec.tau_invariant import TauInvariantCode, rm13

__all__ = [
    "BCH15",
    "CRC",
    "ConvolutionalCode",
    "FEC",
    "Hamming74",
    "LDPC36",
    "PolarAWGN",
    "PolarBEC",
    "ReedSolomon",
    "TauInvariantCode",
    "rm13",
]
