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
from modcode.channel.phase_imbalance import PhaseImbalanceChannel
from modcode.fec.hamming import Hamming74
from modcode.modulation.bpsk import BPSK
from modcode.modulation.qam import QAM
from modcode.modulation.qpsk import QPSK
from modcode.pipeline.chain import Link, LinkResult, run_link

__all__ = [
    "AWGNChannel",
    "BPSK",
    "Hamming74",
    "Link",
    "LinkResult",
    "PhaseImbalanceChannel",
    "QAM",
    "QPSK",
    "run_link",
]
