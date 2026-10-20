# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Digital modulation and forward-error-correction simulation tools."""

from modcode.channel.awgn import AWGNChannel
from modcode.channel.burst import BurstChannel
from modcode.channel.phase_imbalance import PhaseImbalanceChannel
from modcode.channel.rayleigh import BlockRayleighChannel
from modcode.fec.bch import BCH15
from modcode.fec.convolutional import ConvolutionalCode
from modcode.fec.hamming import Hamming74
from modcode.fec.ldpc import LDPC36
from modcode.fec.reed_solomon import ReedSolomon
from modcode.modulation.bpsk import BPSK
from modcode.modulation.lattice import LatticeModulation
from modcode.modulation.qam import QAM
from modcode.modulation.qpsk import QPSK
from modcode.modulation.tcm import TCM8PSK
from modcode.pipeline.chain import Link, LinkResult, run_link
from modcode.pipeline.interleaver import PermutationInterleaver

__all__ = [
    "AWGNChannel",
    "BCH15",
    "BPSK",
    "BlockRayleighChannel",
    "BurstChannel",
    "ConvolutionalCode",
    "Hamming74",
    "Link",
    "LinkResult",
    "LatticeModulation",
    "LDPC36",
    "PermutationInterleaver",
    "PhaseImbalanceChannel",
    "QAM",
    "QPSK",
    "ReedSolomon",
    "TCM8PSK",
    "run_link",
]
