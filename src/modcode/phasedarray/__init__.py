# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Phased-array wavefront synthesis with finite cyclic phase resolution."""

from modcode.phasedarray.array import (
    SPEED_OF_LIGHT,
    TWO_PI,
    LinearArray,
    grating_lobe_limit,
    quantize_phases,
    squint_angle,
)

__all__ = [
    "SPEED_OF_LIGHT",
    "TWO_PI",
    "LinearArray",
    "grating_lobe_limit",
    "quantize_phases",
    "squint_angle",
]
