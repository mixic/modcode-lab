# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Tiny, deliberately insecure Diffie-Hellman demonstration.

The small parameter sizes and lack of protocol protections make this module
unsuitable for protecting real data. Use a reviewed cryptographic library for
applications.
"""

import secrets

from modcode.number_theory.full_reptend import is_prime, multiplicative_order


class ToyDiffieHellman:
    """Demonstrate Diffie-Hellman in the multiplicative group modulo a prime."""

    def __init__(self, prime: int, generator: int) -> None:
        if isinstance(prime, bool) or not isinstance(prime, int):
            raise ValueError("prime must be an integer")
        if isinstance(generator, bool) or not isinstance(generator, int):
            raise ValueError("generator must be an integer")
        if not is_prime(prime):
            raise ValueError("prime must be a prime integer")
        if not 2 <= generator < prime:
            raise ValueError("generator must be in the range [2, prime)")
        if multiplicative_order(generator, prime) != prime - 1:
            raise ValueError("generator must generate the full multiplicative group")
        self.prime = prime
        self.generator = generator

    def generate_private_key(self) -> int:
        """Generate an exponent in the nonzero range for the toy group."""
        if self.prime <= 4:
            raise ValueError("prime must be at least 5 for private-key generation")
        return secrets.randbelow(self.prime - 2) + 1

    def public_key(self, private_key: int) -> int:
        """Compute ``generator ** private_key mod prime``."""
        self._validate_private_key(private_key)
        return pow(self.generator, private_key, self.prime)

    def shared_secret(self, private_key: int, peer_public_key: int) -> int:
        """Compute the shared integer after validating the peer's group element."""
        self._validate_private_key(private_key)
        if isinstance(peer_public_key, bool) or not isinstance(peer_public_key, int):
            raise ValueError("peer public key must be an integer")
        if not 2 <= peer_public_key <= self.prime - 2:
            raise ValueError("peer public key is outside the accepted range")
        return pow(peer_public_key, private_key, self.prime)

    def _validate_private_key(self, private_key: int) -> None:
        if isinstance(private_key, bool) or not isinstance(private_key, int):
            raise ValueError("private key must be an integer")
        if not 1 <= private_key <= self.prime - 2:
            raise ValueError("private key must be in the range [1, prime - 2]")
