# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Tests for the toy QC-MDPC McEliece-style demonstration."""

import numpy as np
import pytest

from modcode.pqc import QCMDPC, circulant, is_circulant


def test_circulant_rows_are_rotations_of_the_first_row() -> None:
    first_row = np.asarray([1, 0, 1, 1, 0], dtype=np.uint8)
    matrix = circulant(first_row)
    assert matrix.shape == (5, 5)
    assert is_circulant(matrix)
    for index in range(5):
        assert np.array_equal(matrix[index], np.roll(first_row, index))


def test_qc_blocks_and_generator_are_circulant() -> None:
    code = QCMDPC(block_length=31, row_weight=3, seed=5)
    assert is_circulant(code.h0)
    assert is_circulant(code.h1)
    assert is_circulant(code.q_block)
    assert code.parity_check.shape == (31, 62)
    assert code.generator.shape == (31, 62)


def test_parity_check_annihilates_the_generator() -> None:
    code = QCMDPC(block_length=31, row_weight=3, seed=5)
    product = (
        code.parity_check.astype(np.int64)
        @ code.generator.T.astype(np.int64)
    ) % 2
    assert not np.any(product)


def test_round_trip_recovers_a_low_weight_error() -> None:
    code = QCMDPC(block_length=31, row_weight=3, seed=5)
    rng = np.random.default_rng(9)
    message = rng.integers(0, 2, size=code.dimension, dtype=np.uint8)
    ciphertext, error = code.encrypt(message, error_weight=2, seed=11)

    result = code.decode(ciphertext, error_weight=2)
    assert np.array_equal(result.message, message)
    assert np.array_equal(result.corrected, code.encode(message))
    assert int(error.sum()) == 2


def test_key_size_is_linear_in_the_block_length() -> None:
    code = QCMDPC(block_length=31, row_weight=3, seed=5)
    sizes = code.key_size_bits
    assert sizes["dense_generator"] == 31 * 31
    assert sizes["qc_generator"] == 31
    assert code.compression_ratio == 31.0


def test_single_error_has_nonzero_syndrome() -> None:
    code = QCMDPC(block_length=31, row_weight=3, seed=5)
    codeword = code.encode(np.zeros(code.dimension, dtype=np.uint8))
    corrupted = codeword.copy()
    corrupted[7] ^= 1
    assert not np.any(code.syndrome(codeword))
    assert np.any(code.syndrome(corrupted))


def test_invalid_parameters_are_rejected() -> None:
    with pytest.raises(ValueError):
        QCMDPC(block_length=2)
    with pytest.raises(ValueError):
        QCMDPC(block_length=31, row_weight=0)
    with pytest.raises(ValueError):
        QCMDPC(block_length=31, row_weight=2)
    with pytest.raises(ValueError):
        QCMDPC(block_length=31, row_weight=31)
    with pytest.raises(ValueError):
        circulant([1, 2])
    assert not is_circulant(np.zeros((3, 2)))

