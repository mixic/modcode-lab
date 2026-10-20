# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Tests for CRC encoding, detection, and division traces."""

import numpy as np
import pytest

from modcode.fec.crc import CRC, CRC8


def _bits_of(text: str) -> np.ndarray:
    return np.unpackbits(np.frombuffer(text.encode("ascii"), dtype=np.uint8))


def test_crc8_check_value_for_123456789() -> None:
    assert CRC8.remainder(_bits_of("123456789")) == 0xF4


def test_encode_appends_remainder_and_zero_check_passes() -> None:
    message = np.asarray([1, 0, 1, 1, 0, 1, 0, 0], dtype=np.uint8)
    codeword = CRC8.encode(message)
    assert codeword.shape == (16,)
    assert np.array_equal(codeword[: message.size], message)
    assert CRC8.check(codeword)
    assert CRC8.syndrome(codeword) == 0


def test_every_single_bit_error_is_detected() -> None:
    message = np.asarray([1, 0, 1, 1, 0, 1, 0, 0], dtype=np.uint8)
    codeword = CRC8.encode(message)
    for position in range(codeword.size):
        corrupted = codeword.copy()
        corrupted[position] ^= 1
        assert not CRC8.check(corrupted)
        assert CRC8.syndrome(corrupted) != 0


@pytest.mark.parametrize("length", [1, 2, 8])
def test_bursts_up_to_width_bits_are_detected(length: int) -> None:
    message = np.asarray([1, 0, 1, 1, 0, 1, 0, 0], dtype=np.uint8)
    codeword = CRC8.encode(message)
    for start in range(codeword.size - length + 1):
        corrupted = codeword.copy()
        corrupted[start : start + length] ^= 1
        assert not CRC8.check(corrupted)


def test_division_trace_shape_and_final_remainder() -> None:
    message = np.asarray([1, 0, 1, 1, 0, 1], dtype=np.uint8)
    trace = CRC8.division_trace(message)
    assert trace.shape == (message.size, 8)
    assert np.array_equal(trace[-1], CRC8.encode(message)[-8:])


def test_invalid_crc_parameters_and_bits() -> None:
    with pytest.raises(ValueError):
        CRC(width=0, polynomial=1)
    with pytest.raises(ValueError):
        CRC(width=8, polynomial=0x100)
    with pytest.raises(ValueError):
        CRC(width=8, polynomial=0x06)
    with pytest.raises(ValueError):
        CRC8.encode([0, 2])
    with pytest.raises(ValueError):
        CRC8.check([0] * 7)
