# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Tests for the deterministic bit interleaver and its Link integration."""

import numpy as np
import pytest

from modcode.channel.awgn import AWGNChannel
from modcode.fec.hamming import Hamming74
from modcode.modulation.qpsk import QPSK
from modcode.pipeline.chain import Link, run_link
from modcode.pipeline.interleaver import PermutationInterleaver


def test_interleaver_is_a_permutation_and_inverts() -> None:
    interleaver = PermutationInterleaver(seed=7)
    values = np.arange(101)
    permuted = interleaver.interleave(values)
    assert np.array_equal(np.sort(permuted), values)
    assert np.array_equal(interleaver.deinterleave(permuted), values)


def test_interleaver_is_deterministic() -> None:
    first = PermutationInterleaver(seed=3)
    second = PermutationInterleaver(seed=3)
    values = np.arange(64)
    assert np.array_equal(first.interleave(values), second.interleave(values))
    assert not np.array_equal(
        first.interleave(values),
        PermutationInterleaver(seed=4).interleave(values),
    )


def test_interleaved_link_decodes_with_soft_llrs() -> None:
    result = run_link(
        Link(
            QPSK(),
            AWGNChannel(),
            Hamming74(),
            decoding="soft",
            interleaver=PermutationInterleaver(seed=11),
        ),
        n_bits=20_000,
        ebn0_db=8.0,
        rng=np.random.default_rng(2),
    )
    assert result.decoding == "soft"
    assert result.ber < 0.01


def test_interleaver_requires_fec() -> None:
    with pytest.raises(ValueError, match="interleaver"):
        Link(QPSK(), AWGNChannel(), interleaver=PermutationInterleaver())


def test_interleaver_rejects_non_vector_input() -> None:
    interleaver = PermutationInterleaver()
    with pytest.raises(ValueError, match="one-dimensional"):
        interleaver.interleave(np.zeros((2, 2)))
    with pytest.raises(ValueError, match="nonempty"):
        interleaver.deinterleave([])
