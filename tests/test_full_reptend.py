# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Finite arithmetic checks for full-reptend primes."""

import pytest

from modcode.number_theory.full_reptend import (
    cyclic_multiples,
    full_reptend_cycle,
    full_reptend_primes,
    is_full_reptend_prime,
    is_prime,
    multiplicative_order,
)


def test_multiplicative_order_matches_modular_powers() -> None:
    assert multiplicative_order(10, 7) == 6
    assert all(pow(10, exponent, 7) != 1 for exponent in range(1, 6))
    assert pow(10, 6, 7) == 1


def test_seven_generates_the_142857_cyclic_number() -> None:
    assert is_full_reptend_prime(7)
    assert full_reptend_cycle(7) == "142857"
    assert cyclic_multiples(7) == (
        "142857",
        "285714",
        "428571",
        "571428",
        "714285",
        "857142",
    )


def test_full_reptend_cycle_preserves_leading_zero() -> None:
    assert full_reptend_cycle(17) == "0588235294117647"


def test_full_reptend_prime_search_is_finite_and_order_based() -> None:
    assert full_reptend_primes(30) == (7, 17, 19, 23, 29)
    assert not is_full_reptend_prime(11)
    assert not is_full_reptend_prime(21)


def test_order_requires_a_unit() -> None:
    with pytest.raises(ValueError):
        multiplicative_order(2, 8)


def test_cycle_rejects_non_full_reptend_primes() -> None:
    with pytest.raises(ValueError):
        full_reptend_cycle(11)


def test_prime_predicate_handles_small_values() -> None:
    assert not is_prime(-1)
    assert not is_prime(1)
    assert is_prime(2)
