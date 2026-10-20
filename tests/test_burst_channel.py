# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Tests for the synthetic burst channel."""

import numpy as np
import pytest

from modcode.channel.burst import BurstChannel


def test_burst_channel_is_identity_without_noise_or_burst() -> None:
    channel = BurstChannel(burst_length=3, burst_probability=0.0)
    symbols = np.asarray([1 + 2j, -1 + 0j])
    received = channel.transmit(symbols, 0.0, rng=np.random.default_rng(1))
    np.testing.assert_allclose(received, symbols)
    assert channel.last_burst is None


def test_burst_erases_a_contiguous_interval() -> None:
    channel = BurstChannel(burst_length=4, burst_probability=1.0)
    received = channel.transmit(
        np.ones(20, dtype=complex), 0.0, rng=np.random.default_rng(2)
    )
    assert channel.last_burst is not None
    start, stop = channel.last_burst
    assert stop - start <= 4
    np.testing.assert_allclose(received[start:stop], 0.0)
    assert not np.allclose(received, 0.0)


def test_burst_channel_validation() -> None:
    with pytest.raises(ValueError, match="burst_length"):
        BurstChannel(burst_length=0)
    with pytest.raises(ValueError, match="burst_probability"):
        BurstChannel(burst_probability=1.5)
