# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Tests for the explicitly tail-biting 8-PSK trellis-coded example."""

import numpy as np
import pytest

from modcode.modulation.cyclic_tcm import TailBitingTCM8PSK


def test_tail_biting_encoder_is_equivariant_to_pairwise_cyclic_shifts() -> None:
    encoder = TailBitingTCM8PSK()
    bits = np.array([0, 0, 1, 0, 1, 1, 0, 1], dtype=np.uint8)
    pairs = bits.reshape(-1, 2)
    shifted = np.roll(pairs, 1, axis=0).reshape(-1)
    assert np.allclose(encoder.encode(shifted), np.roll(encoder.encode(bits), 1))


@pytest.mark.parametrize("n_pairs", [1, 2, 5, 12])
def test_tail_biting_viterbi_decodes_noiseless_sequence(n_pairs: int) -> None:
    encoder = TailBitingTCM8PSK()
    bits = np.random.default_rng(n_pairs).integers(
        0, 2, size=2 * n_pairs, dtype=np.uint8
    )
    assert np.array_equal(encoder.decode(encoder.encode(bits)), bits)


def test_tail_biting_viterbi_corrects_small_perturbations() -> None:
    encoder = TailBitingTCM8PSK()
    bits = np.array([1, 0, 0, 1, 1, 1, 0, 0], dtype=np.uint8)
    symbols = encoder.encode(bits)
    noisy = symbols + 0.01 * (
        np.random.default_rng(3).normal(size=symbols.size)
        + 1j * np.random.default_rng(4).normal(size=symbols.size)
    )
    assert np.array_equal(encoder.decode(noisy), bits)


def test_encoder_requires_complete_nonempty_input_pairs() -> None:
    encoder = TailBitingTCM8PSK()
    with pytest.raises(ValueError):
        encoder.encode(np.array([], dtype=np.uint8))
    with pytest.raises(ValueError):
        encoder.encode(np.array([1], dtype=np.uint8))
