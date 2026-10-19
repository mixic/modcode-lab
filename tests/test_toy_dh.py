# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Tests of the deliberately insecure, small-parameter DH demonstration."""

import pytest

from modcode.cryptography.toy_dh import ToyDiffieHellman


def test_two_parties_derive_the_same_shared_group_element() -> None:
    dh = ToyDiffieHellman(prime=23, generator=5)
    alice_private, bob_private = 6, 15
    alice_public = dh.public_key(alice_private)
    bob_public = dh.public_key(bob_private)
    assert dh.shared_secret(alice_private, bob_public) == pow(5, 90, 23)
    assert dh.shared_secret(bob_private, alice_public) == pow(5, 90, 23)


@pytest.mark.parametrize(
    ("prime", "generator"),
    [(21, 2), (23, 1), (23, 4)],
)
def test_parameters_require_a_prime_and_primitive_root(
    prime: int, generator: int
) -> None:
    with pytest.raises(ValueError):
        ToyDiffieHellman(prime, generator)


def test_private_and_peer_values_are_validated() -> None:
    dh = ToyDiffieHellman(23, 5)
    with pytest.raises(ValueError):
        dh.public_key(0)
    with pytest.raises(ValueError):
        dh.shared_secret(3, 1)


def test_private_key_generator_returns_valid_exponents() -> None:
    dh = ToyDiffieHellman(23, 5)
    assert 1 <= dh.generate_private_key() <= 21
