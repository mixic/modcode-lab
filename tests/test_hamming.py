# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Tests for cyclic Hamming(7,4)."""

import numpy as np
import pytest

from modcode.fec.hamming import Hamming74


def test_known_systematic_codeword() -> None:
    code = Hamming74()
    assert np.array_equal(code.encode(np.array([1, 0, 1, 0])), [1, 0, 1, 0, 0, 1, 1])


def test_corrects_each_single_bit_error() -> None:
    code = Hamming74()
    data = np.array([1, 0, 1, 0], dtype=np.uint8)
    codeword = code.encode(data)
    for position in range(7):
        corrupted = codeword.copy()
        corrupted[position] ^= 1
        assert np.array_equal(code.decode(corrupted), data)


def test_round_trips_multiple_blocks() -> None:
    code = Hamming74()
    data = np.tile(np.array([1, 0, 1, 0], dtype=np.uint8), 3)
    assert np.array_equal(code.decode(code.encode(data)), data)


def test_syndrome_table_covers_every_nonzero_syndrome() -> None:
    code = Hamming74()
    assert np.all(code._syndrome_to_position[1:] > 0)
    assert code._syndrome_to_position[0] == 0


@pytest.mark.parametrize("bits", [[0, 1], [0, 1, 0, 1, 0]])
def test_encode_requires_complete_blocks(bits: list[int]) -> None:
    with pytest.raises(ValueError):
        Hamming74().encode(np.asarray(bits))
