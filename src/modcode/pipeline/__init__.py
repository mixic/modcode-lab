# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""End-to-end communication link simulation."""

from modcode.pipeline.chain import Link, LinkResult, run_link
from modcode.pipeline.interleaver import PermutationInterleaver
from modcode.pipeline.ofdm import OFDMLink, OFDMResult, run_ofdm_link

__all__ = [
    "Link",
    "LinkResult",
    "OFDMLink",
    "OFDMResult",
    "PermutationInterleaver",
    "run_link",
    "run_ofdm_link",
]
