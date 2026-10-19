# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Finite computations for multiplicative orders and repeating fractions."""

from modcode.number_theory.full_reptend import (
    cyclic_multiples,
    full_reptend_cycle,
    full_reptend_primes,
    is_full_reptend_prime,
    is_prime,
    multiplicative_order,
)

__all__ = [
    "cyclic_multiples",
    "full_reptend_cycle",
    "full_reptend_primes",
    "is_full_reptend_prime",
    "is_prime",
    "multiplicative_order",
]
