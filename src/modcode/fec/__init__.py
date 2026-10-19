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
from modcode.fec.hamming import Hamming74

__all__ = ["FEC", "Hamming74"]
