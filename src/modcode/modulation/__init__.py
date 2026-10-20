# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Digital modulation schemes."""

from modcode.modulation.apsk import APSK
from modcode.modulation.base import Modulation
from modcode.modulation.bpsk import BPSK
from modcode.modulation.ct_im import CTIM, ConstantWeightCyclicCode
from modcode.modulation.cyclic_tcm import TailBitingTCM8PSK
from modcode.modulation.involution import (
    block_tau_hypothesis,
    correct_tau,
    pilot_tau_hypothesis,
    tau_log_likelihood_ratio,
)
from modcode.modulation.qam import QAM
from modcode.modulation.qpsk import QPSK

__all__ = [
    "APSK",
    "BPSK",
    "CTIM",
    "ConstantWeightCyclicCode",
    "Modulation",
    "QAM",
    "QPSK",
    "TailBitingTCM8PSK",
    "block_tau_hypothesis",
    "correct_tau",
    "pilot_tau_hypothesis",
    "tau_log_likelihood_ratio",
]
