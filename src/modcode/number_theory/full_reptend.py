# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Explore full-reptend primes and their repeating fraction cycles.

These functions perform finite computations. They do not establish
Artin's primitive-root conjecture or any Riemann-hypothesis statement.
"""

from numbers import Integral

_DIGITS = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"


def multiplicative_order(value: int, modulus: int) -> int:
    """Return the order of ``value`` in the units modulo ``modulus``."""
    if isinstance(value, bool) or not isinstance(value, Integral):
        raise ValueError("value must be an integer")
    if isinstance(modulus, bool) or not isinstance(modulus, Integral):
        raise ValueError("modulus must be an integer")
    value = int(value)
    modulus = int(modulus)
    if modulus < 2:
        raise ValueError("modulus must be at least 2")
    if not 0 <= value < modulus:
        value %= modulus
    if _gcd(value, modulus) != 1:
        raise ValueError("value must be coprime to modulus")
    current = 1
    for order in range(1, modulus + 1):
        current = (current * value) % modulus
        if current == 1:
            return order
    raise ArithmeticError("multiplicative order search did not terminate")


def is_full_reptend_prime(prime: int, base: int = 10) -> bool:
    """Return whether ``1 / prime`` has period ``prime - 1`` in ``base``."""
    if isinstance(base, bool) or not isinstance(base, Integral):
        raise ValueError("base must be an integer")
    base = int(base)
    if base < 2 or base > len(_DIGITS):
        raise ValueError(f"base must be between 2 and {len(_DIGITS)}")
    if not is_prime(prime) or _gcd(prime, base) != 1:
        return False
    return multiplicative_order(base, prime) == prime - 1


def full_reptend_cycle(prime: int, base: int = 10) -> str:
    """Return the base-``base`` repetend of ``1 / prime``.

    The result is padded to exactly ``prime - 1`` digits and uses uppercase
    letters for digit values 10 through 35.
    """
    if not is_full_reptend_prime(prime, base):
        raise ValueError(f"{prime} is not a full-reptend prime in base {base}")
    remainder = 1
    digits: list[str] = []
    for _ in range(prime - 1):
        remainder *= base
        digit, remainder = divmod(remainder, prime)
        digits.append(_DIGITS[digit])
    if remainder != 1:
        raise ArithmeticError("unexpected remainder while generating repetend")
    return "".join(digits)


def cyclic_multiples(prime: int, base: int = 10) -> tuple[str, ...]:
    """Return the zero-padded base-``base`` multiples of a full-reptend cycle."""
    cycle = full_reptend_cycle(prime, base)
    value = _digits_to_integer(cycle, base)
    width = prime - 1
    return tuple(
        _integer_to_digits(value * multiplier, base, width)
        for multiplier in range(1, prime)
    )


def full_reptend_primes(limit: int, base: int = 10) -> tuple[int, ...]:
    """List full-reptend primes below ``limit`` in the selected base."""
    if isinstance(limit, bool) or not isinstance(limit, Integral):
        raise ValueError("limit must be an integer")
    if isinstance(base, bool) or not isinstance(base, Integral):
        raise ValueError("base must be an integer")
    limit = int(limit)
    base = int(base)
    if limit < 0:
        raise ValueError("limit must be non-negative")
    if base < 2 or base > len(_DIGITS):
        raise ValueError(f"base must be between 2 and {len(_DIGITS)}")
    return tuple(
        prime
        for prime in range(2, limit)
        if is_prime(prime) and is_full_reptend_prime(prime, base)
    )


def is_prime(value: int) -> bool:
    """Return whether ``value`` is prime using trial division."""
    if isinstance(value, bool) or not isinstance(value, Integral):
        raise ValueError("value must be an integer")
    value = int(value)
    if value < 2:
        return False
    if value % 2 == 0:
        return value == 2
    divisor = 3
    while divisor * divisor <= value:
        if value % divisor == 0:
            return False
        divisor += 2
    return True


def _gcd(first: int, second: int) -> int:
    while second:
        first, second = second, first % second
    return abs(first)


def _digits_to_integer(digits: str, base: int) -> int:
    result = 0
    for digit in digits:
        result = result * base + _DIGITS.index(digit)
    return result


def _integer_to_digits(value: int, base: int, width: int) -> str:
    output: list[str] = []
    for _ in range(width):
        value, digit = divmod(value, base)
        output.append(_DIGITS[digit])
    if value:
        raise ArithmeticError("multiple does not fit the repetend width")
    return "".join(reversed(output))
