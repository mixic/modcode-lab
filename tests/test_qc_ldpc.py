# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Tests that expose the circulant identity blocks in QC-LDPC matrices."""

import numpy as np
import pytest

from modcode.fec.qc_ldpc import circulant_identity, expand_qc_base_matrix


def test_shift_value_expands_to_cyclically_shifted_identity() -> None:
    identity = np.eye(4, dtype=np.uint8)
    assert np.array_equal(circulant_identity(4, 0), identity)
    assert np.array_equal(circulant_identity(4, 1), np.roll(identity, 1, axis=1))
    assert np.array_equal(circulant_identity(4, 5), np.roll(identity, 1, axis=1))


def test_five_by_five_shift_two_matches_documented_example() -> None:
    assert np.array_equal(
        circulant_identity(5, 2),
        np.array(
            [
                [0, 0, 1, 0, 0],
                [0, 0, 0, 1, 0],
                [0, 0, 0, 0, 1],
                [1, 0, 0, 0, 0],
                [0, 1, 0, 0, 0],
            ],
            dtype=np.uint8,
        ),
    )


def test_qc_base_matrix_expands_to_expected_block_matrix() -> None:
    base = np.array([[0, 1, -1], [2, -1, 0]], dtype=np.int64)
    expanded = expand_qc_base_matrix(base, lifting_factor=3)
    identity = np.eye(3, dtype=np.uint8)
    expected = np.block(
        [
            [identity, np.roll(identity, 1, axis=1), np.zeros((3, 3), dtype=np.uint8)],
            [np.roll(identity, 2, axis=1), np.zeros((3, 3), dtype=np.uint8), identity],
        ]
    )
    assert expanded.shape == (6, 9)
    assert np.array_equal(expanded, expected)


def test_each_nonzero_block_is_a_circulant_permutation_matrix() -> None:
    base = np.array([[0, 2], [1, -1]], dtype=np.int64)
    lifting_factor = 4
    matrix = expand_qc_base_matrix(base, lifting_factor)
    for block_row in range(base.shape[0]):
        for block_column in range(base.shape[1]):
            block = matrix[
                block_row * lifting_factor : (block_row + 1) * lifting_factor,
                block_column * lifting_factor : (block_column + 1) * lifting_factor,
            ]
            if base[block_row, block_column] == -1:
                assert not np.any(block)
            else:
                assert np.all(block.sum(axis=0) == 1)
                assert np.all(block.sum(axis=1) == 1)


@pytest.mark.parametrize(
    ("base", "lifting_factor"),
    [
        (np.array([[0, 3]], dtype=np.int64), 3),
        (np.array([[-2]], dtype=np.int64), 3),
        (np.array([[0.0]]), 3),
    ],
)
def test_qc_expansion_rejects_invalid_shift_matrices(
    base: np.ndarray, lifting_factor: int
) -> None:
    with pytest.raises(ValueError):
        expand_qc_base_matrix(base, lifting_factor)
