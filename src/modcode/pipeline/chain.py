# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Composition of FEC, modulation, and channel into a simulated link."""

import math
from dataclasses import dataclass

import numpy as np

from modcode.channel.base import Channel
from modcode.fec.base import FEC
from modcode.modulation.base import BitArray, Modulation, SymbolArray


@dataclass(frozen=True)
class Link:
    """A transmit/receive chain for one modulation, FEC, and channel."""

    modulation: Modulation
    channel: Channel
    fec: FEC | None = None


@dataclass(frozen=True)
class LinkResult:
    """Bit and symbol arrays captured during a single link simulation."""

    source_bits: BitArray
    received_bits: BitArray
    transmitted_bits: BitArray
    transmitted_symbols: SymbolArray
    received_symbols: SymbolArray

    @property
    def bit_errors(self) -> int:
        return int(np.count_nonzero(self.source_bits != self.received_bits))

    @property
    def ber(self) -> float:
        if self.source_bits.size == 0:
            return 0.0
        return self.bit_errors / self.source_bits.size


def run_link(
    link: Link,
    n_bits: int = 10_000,
    ebn0_db: float = 5.0,
    rng: np.random.Generator | None = None,
) -> LinkResult:
    """Generate random data and pass it through the configured link."""
    if n_bits <= 0:
        raise ValueError("n_bits must be positive")
    if not np.isfinite(ebn0_db):
        raise ValueError("ebn0_db must be finite")
    generator = rng if rng is not None else np.random.default_rng()
    source = generator.integers(0, 2, size=n_bits, dtype=np.uint8)

    block_size = link.modulation.bits_per_symbol
    if link.fec is not None:
        blocks_per_symbol_multiple = block_size // math.gcd(
            link.fec.encoded_bits_per_block, block_size
        )
        block_size = link.fec.data_bits_per_block * blocks_per_symbol_multiple
    padded_size = math.ceil(n_bits / block_size) * block_size
    padded = np.pad(source, (0, padded_size - n_bits))
    transmitted = link.fec.encode(padded) if link.fec else padded
    tx_symbols = link.modulation.modulate(transmitted)
    rx_symbols = link.channel.transmit(
        tx_symbols,
        ebn0_db,
        link.modulation.bits_per_symbol,
        link.fec.rate if link.fec else 1.0,
        generator,
    )
    demodulated = link.modulation.demodulate(rx_symbols)
    decoded = link.fec.decode(demodulated) if link.fec else demodulated
    return LinkResult(
        source_bits=source,
        received_bits=decoded[:n_bits],
        transmitted_bits=transmitted,
        transmitted_symbols=tx_symbols,
        received_symbols=rx_symbols,
    )
