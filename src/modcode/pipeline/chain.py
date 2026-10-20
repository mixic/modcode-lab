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
    decoding: str = "auto"

    def __post_init__(self) -> None:
        if self.decoding not in {"auto", "hard", "soft"}:
            raise ValueError("decoding must be 'auto', 'hard', or 'soft'")


@dataclass(frozen=True)
class LinkResult:
    """Bit and symbol arrays captured during a single link simulation."""

    source_bits: BitArray
    received_bits: BitArray
    transmitted_bits: BitArray
    transmitted_symbols: SymbolArray
    received_symbols: SymbolArray
    decoding: str = "hard"

    @property
    def bit_errors(self) -> int:
        return int(np.count_nonzero(self.source_bits != self.received_bits))

    @property
    def ber(self) -> float:
        if self.source_bits.size == 0:
            return 0.0
        return self.bit_errors / self.source_bits.size


def _decode(
    link: Link, received_symbols: SymbolArray, ebn0_db: float
) -> tuple[BitArray, str]:
    """Decode with soft LLRs when both the channel and FEC support them."""
    if link.fec is None:
        return link.modulation.demodulate(received_symbols), "hard"

    if link.decoding == "hard":
        return (
            link.fec.decode(link.modulation.demodulate(received_symbols)),
            "hard",
        )

    try:
        variance = link.channel.noise_variance(
            ebn0_db,
            link.modulation.bits_per_symbol,
            link.fec.rate,
        )
    except NotImplementedError:
        if link.decoding == "soft":
            raise ValueError(
                f"{type(link.channel).__name__} does not expose a noise "
                "variance; soft decoding is unavailable"
            ) from None
        return (
            link.fec.decode(link.modulation.demodulate(received_symbols)),
            "hard",
        )

    try:
        llrs = link.modulation.soft_demodulate(received_symbols, variance)
        return link.fec.decode_soft(llrs), "soft"
    except NotImplementedError:
        if link.decoding == "soft":
            raise ValueError(
                f"{type(link.fec).__name__} does not implement soft decoding"
            ) from None
        return (
            link.fec.decode(link.modulation.demodulate(received_symbols)),
            "hard",
        )


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
    decoded, decoding = _decode(link, rx_symbols, ebn0_db)
    return LinkResult(
        source_bits=source,
        received_bits=decoded[:n_bits],
        transmitted_bits=transmitted,
        transmitted_symbols=tx_symbols,
        received_symbols=rx_symbols,
        decoding=decoding,
    )
