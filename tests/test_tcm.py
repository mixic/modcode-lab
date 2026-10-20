# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Tests for the four-state 8-PSK TCM example."""

import numpy as np
import pytest

from modcode.channel.awgn import AWGNChannel
from modcode.modulation.qpsk import QPSK
from modcode.modulation.tcm import TCM8PSK
from modcode.pipeline.chain import Link, run_link


def test_tcm_noiseless_round_trip() -> None:
    tcm = TCM8PSK()
    rng = np.random.default_rng(1)
    for pairs in (1, 2, 10, 200):
        bits = rng.integers(0, 2, 2 * pairs, dtype=np.uint8)
        np.testing.assert_array_equal(tcm.decode(tcm.encode(bits)), bits)


def test_tcm_beats_uncoded_qpsk_at_6db() -> None:
    tcm = TCM8PSK()
    rng = np.random.default_rng(2)
    bits = rng.integers(0, 2, 40_000, dtype=np.uint8)
    symbols = tcm.encode(bits)
    variance = AWGNChannel.noise_variance(6.0, 2)
    received = AWGNChannel().transmit(
        symbols, variance, np.random.default_rng(3)
    )
    tcm_ber = float(np.mean(tcm.decode(received) != bits))
    qpsk_ber = run_link(
        Link(QPSK(), AWGNChannel()),
        n_bits=40_000,
        ebn0_db=6.0,
        rng=np.random.default_rng(3),
    ).ber
    assert tcm_ber < qpsk_ber


def test_tcm_rejects_invalid_input() -> None:
    with pytest.raises(ValueError, match="even"):
        TCM8PSK().encode(np.asarray([0, 1, 1]))
    with pytest.raises(ValueError, match="nonempty"):
        TCM8PSK().decode(np.empty(0, dtype=complex))
    with pytest.raises(ValueError, match="finite"):
        TCM8PSK().decode(np.asarray([np.inf]))
