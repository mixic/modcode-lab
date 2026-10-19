# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Experimental optical intensity-modulation tools."""

from modcode.optical.channel import IMDDChannel
from modcode.optical.modulation import OOK, IntensityPAM, NonUniformPAM
from modcode.optical.robustness import ber_sweep_imdd, optimize_pam_levels

__all__ = [
    "IMDDChannel",
    "IntensityPAM",
    "NonUniformPAM",
    "OOK",
    "ber_sweep_imdd",
    "optimize_pam_levels",
]
