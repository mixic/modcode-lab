# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Densities of primes with a prescribed residual index.

For a prime ``p`` not dividing the base ``g``, define the residual index

    c_p = (p - 1) / ord_p(g).

Artin's primitive-root conjecture is the case ``c_p = 1``. The case
``c_p = m`` is a half-period or higher-index condition. Under GRH, the density
is given by the Hooley--Lenstra field-degree sum

    delta(g, m) = sum_n mu(n) / [Q(zeta_{mn}, g^(1/(mn))) : Q].

For squarefree ``g`` the degree is ``mn * phi(mn)`` except when the quadratic
field ``Q(sqrt(g))`` is already contained in the cyclotomic field, in which
case it drops by a factor of two. For ``g = 10`` this happens exactly when
``40 | mn``.

This is a numerical/experimental module. The density formula is conditional on
GRH; the finite prime counts are unconditional but converge slowly.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

import numpy as np
from numpy.typing import NDArray

from modcode.number_theory.full_reptend import is_prime


def artin_constant(terms: int = 1000) -> float:
    """Return the Artin constant from the first ``terms`` primes.

    ``A = product_q (1 - 1/(q(q-1)))``. The product converges quickly; 1,000
    primes already give far more digits than the plots need.
    """
    if isinstance(terms, bool) or not isinstance(terms, int) or terms < 1:
        raise ValueError("terms must be a positive integer")
    product = 1.0
    found = 0
    candidate = 2
    while found < terms:
        if is_prime(candidate):
            product *= 1.0 - 1.0 / (candidate * (candidate - 1))
            found += 1
        candidate += 1
    return product


def _mobius_sieve(limit: int) -> NDArray[np.int8]:
    mu = np.ones(limit + 1, dtype=np.int8)
    mu[0] = 0
    composite = bytearray(limit + 1)
    for prime in range(2, limit + 1):
        if composite[prime]:
            continue
        for multiple in range(prime, limit + 1, prime):
            composite[multiple] = 1
            mu[multiple] *= -1
        square = prime * prime
        if square <= limit:
            for multiple in range(square, limit + 1, square):
                mu[multiple] = 0
    return mu


def _totient_sieve(limit: int) -> NDArray[np.int64]:
    phi = np.arange(limit + 1, dtype=np.int64)
    if limit >= 1:
        phi[1] = 1
    for value in range(2, limit + 1):
        if phi[value] == value:
            for multiple in range(value, limit + 1, value):
                phi[multiple] -= phi[multiple] // value
    return phi


def residual_index_density(
    m: int, *, base: int = 10, terms: int = 200_000
) -> float:
    """Return the GRH density of primes with residual index exactly ``m``.

    Only ``base=10`` is supported, because the degree correction uses the
    conductor of ``Q(sqrt(10))``. ``terms`` controls the truncation of the
    Möbius sum.
    """
    if isinstance(m, bool) or not isinstance(m, int) or m < 1:
        raise ValueError("m must be a positive integer")
    if base != 10:
        raise ValueError("only base 10 is supported by this implementation")
    if isinstance(terms, bool) or not isinstance(terms, int) or terms < 1:
        raise ValueError("terms must be a positive integer")

    mu = _mobius_sieve(terms)
    phi = _totient_sieve(m * terms)
    density = 0.0
    for n in range(1, terms + 1):
        if mu[n] == 0:
            continue
        multiple = m * n
        degree = multiple * int(phi[multiple])
        if multiple % 40 == 0:
            degree //= 2
        density += float(mu[n]) / degree
    return density


def known_decimal_density_multiple(m: int) -> Fraction | None:
    """Return the exact rational multiple of Artin's constant for base 10.

    Values are available for the cases derived in the module docstring:
    ``m=1`` (Artin), ``m=2`` (half period), and ``m=4``. The ``m=4`` value is
    ``29/152``, not the sometimes-quoted ``29/52``.
    """
    if isinstance(m, bool) or not isinstance(m, int) or m < 1:
        raise ValueError("m must be a positive integer")
    return {1: Fraction(1, 1), 2: Fraction(3, 4), 4: Fraction(29, 152)}.get(m)


@dataclass(frozen=True)
class ResidualIndexCounts:
    """Finite count of primes by residual index."""

    limit: int
    counts: dict[int, int]
    total: int

    def ratio(self, m: int) -> float:
        """Return the observed fraction of primes with residual index ``m``."""
        return self.counts.get(m, 0) / self.total if self.total else 0.0


def _prime_factors(value: int) -> list[int]:
    factors: list[int] = []
    if value % 2 == 0:
        factors.append(2)
        while value % 2 == 0:
            value //= 2
    divisor = 3
    while divisor * divisor <= value:
        if value % divisor == 0:
            factors.append(divisor)
            while value % divisor == 0:
                value //= divisor
        divisor += 2
    if value > 1:
        factors.append(value)
    return factors


def _multiplicative_order(base: int, prime: int, factors: list[int]) -> int:
    order = prime - 1
    for factor in factors:
        while order % factor == 0 and pow(base, order // factor, prime) == 1:
            order //= factor
    return order


def count_residual_indices(
    limit: int, *, base: int = 10
) -> ResidualIndexCounts:
    """Count primes ``p <= limit`` by the residual index ``(p-1)/ord_p(base)``."""
    if isinstance(limit, bool) or not isinstance(limit, int) or limit < 3:
        raise ValueError("limit must be an integer at least 3")
    if base != 10:
        raise ValueError("only base 10 is supported by this implementation")

    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[0:2] = b"\x00\x00"
    for value in range(2, int(limit**0.5) + 1):
        if sieve[value]:
            start = value * value
            sieve[start : limit + 1 : value] = b"\x00" * (
                ((limit - start) // value) + 1
            )

    counts: dict[int, int] = {}
    total = 0
    for prime in range(2, limit + 1):
        if not sieve[prime] or prime in (2, 5):
            continue
        total += 1
        factors = _prime_factors(prime - 1)
        order = _multiplicative_order(base, prime, factors)
        index = (prime - 1) // order
        counts[index] = counts.get(index, 0) + 1
    return ResidualIndexCounts(limit=limit, counts=counts, total=total)


def residual_index_convergence(
    limits: list[int], *, base: int = 10
) -> dict[int, NDArray[np.float64]]:
    """Return observed ratios for indices 1, 2, and 4 at each limit."""
    if not limits:
        raise ValueError("limits must be nonempty")
    series: dict[int, list[float]] = {1: [], 2: [], 4: []}
    for limit in limits:
        counts = count_residual_indices(limit, base=base)
        for index in series:
            series[index].append(counts.ratio(index))
    return {index: np.asarray(values) for index, values in series.items()}


__all__ = [
    "ResidualIndexCounts",
    "artin_constant",
    "count_residual_indices",
    "known_decimal_density_multiple",
    "residual_index_convergence",
    "residual_index_density",
]
