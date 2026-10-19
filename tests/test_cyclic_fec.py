# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Cyclic-shift closure checks for the Hamming code and a shortened variant."""

import itertools

import numpy as np

from modcode.fec.cyclic import cyclic_shift
from modcode.fec.hamming import Hamming74


def test_every_hamming74_codeword_is_closed_under_cyclic_shifts() -> None:
    code = Hamming74()
    data_words = np.asarray(list(itertools.product((0, 1), repeat=4)), dtype=np.uint8)
    codewords = {tuple(code.encode(word).tolist()) for word in data_words}
    assert len(codewords) == 16
    for word in codewords:
        for shift in range(7):
            assert tuple(cyclic_shift(np.asarray(word), shift).tolist()) in codewords


def test_shortening_hamming_code_breaks_six_bit_cyclic_closure() -> None:
    code = Hamming74()
    data_words = np.asarray(list(itertools.product((0, 1), repeat=3)), dtype=np.uint8)
    shortened_words = {
        tuple(code.encode(np.concatenate(([0], word)))[1:].tolist())
        for word in data_words
    }
    assert len(shortened_words) == 8
    assert any(
        tuple(cyclic_shift(np.asarray(word)).tolist()) not in shortened_words
        for word in shortened_words
    )
