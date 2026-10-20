# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Tests for the CT-IM constant-weight index code and mapper."""

import numpy as np
import pytest

from modcode.modulation.ct_im import CTIM, ConstantWeightCyclicCode


def test_codebook_has_seven_weight_three_words() -> None:
    code = ConstantWeightCyclicCode()
    assert len(code) == 7
    assert code.patterns.shape == (7, 7)
    assert np.all(code.patterns.sum(axis=1) == 3)


def test_codebook_is_one_cyclic_orbit() -> None:
    code = ConstantWeightCyclicCode()
    for label in range(len(code)):
        shifted = {
            tuple(np.roll(code.encode(label), steps).tolist())
            for steps in range(7)
        }
        assert shifted == {tuple(word.tolist()) for word in code.patterns}


def test_every_pair_of_words_has_distance_four() -> None:
    code = ConstantWeightCyclicCode()
    matrix = code.distance_matrix()
    off_diagonal = matrix[~np.eye(7, dtype=bool)]
    assert np.all(off_diagonal == 4)


def test_every_single_index_error_is_corrected() -> None:
    code = ConstantWeightCyclicCode()
    for label in range(len(code)):
        word = code.encode(label)
        for position in range(7):
            corrupted = word.copy()
            corrupted[position] ^= 1
            decoded, distance = code.decode(corrupted)
            assert distance <= 1
            assert decoded == label


def test_ctim_round_trip_without_noise() -> None:
    ct = CTIM(M=4)
    rng = np.random.default_rng(2026)
    for _ in range(20):
        bits = rng.integers(0, 2, size=ct.bits_per_block, dtype=np.uint8)
        label = (bits[0] << 2) | (bits[1] << 1) | bits[2]
        if label == 7:
            bits[:3] = [0, 1, 1]
        recovered = ct.decode(ct.encode(bits))
        assert np.array_equal(recovered, bits)


def test_single_subcarrier_error_does_not_break_the_index() -> None:
    ct = CTIM(M=4)
    rng = np.random.default_rng(7)
    bits = rng.integers(0, 2, size=ct.bits_per_block, dtype=np.uint8)
    vector = ct.encode(bits)
    corrupted = vector.copy()
    position = int(np.argmax(np.abs(corrupted)))
    corrupted[position] *= np.exp(1j * np.pi)
    label, distance = ct.decode_index(corrupted)
    original_label = (bits[0] << 2) | (bits[1] << 1) | bits[2]
    assert distance <= 1
    assert label == original_label


def test_single_activity_error_is_corrected_in_the_index() -> None:
    ct = CTIM(M=4)
    rng = np.random.default_rng(13)
    bits = rng.integers(0, 2, size=ct.bits_per_block, dtype=np.uint8)
    bits[:3] = [0, 0, 0]
    vector = ct.encode(bits)
    original_label = int((bits[0] << 2) | (bits[1] << 1) | bits[2])
    pattern = ct.index_code.encode(original_label)

    silent = int(np.flatnonzero(pattern == 0)[0])
    switched_on = vector.copy()
    switched_on[silent] = 1.0 + 0.0j
    label_on, distance_on = ct.decode_index(switched_on)
    assert distance_on == 1
    assert label_on == original_label

    active = int(np.flatnonzero(pattern == 1)[0])
    switched_off = vector.copy()
    switched_off[active] = 0.0
    label_off, distance_off = ct.decode_index(switched_off)
    assert distance_off == 1
    assert label_off == original_label


def test_spectral_efficiency_matches_the_rate_budget() -> None:
    ct = CTIM(M=4)
    assert ct.bits_per_block == 9
    assert ct.spectral_efficiency == pytest.approx(9 / 7)


def test_shift_selection_lowers_or_preserves_papr() -> None:
    ct = CTIM(M=4, cyclic_shift=True)
    rng = np.random.default_rng(11)
    bits = rng.integers(0, 2, size=ct.bits_per_block, dtype=np.uint8)
    bits[:3] = [0, 0, 0]
    vector, shift = ct.encode(bits, return_shift=True)
    unshifted = np.roll(vector, -shift)
    assert ct.papr_db(vector) <= ct.papr_db(unshifted) + 1e-12


def test_all_cyclic_shifts_are_valid_codewords() -> None:
    code = ConstantWeightCyclicCode()
    for label in range(len(code)):
        word = code.encode(label)
        for steps in range(7):
            assert code.is_codeword(np.roll(word, steps))


def test_invalid_ctim_inputs() -> None:
    ct = CTIM(M=4)
    with pytest.raises(ValueError):
        ct.encode(np.zeros(ct.bits_per_block - 1, dtype=np.uint8))
    with pytest.raises(ValueError):
        ct.decode(np.zeros(6, dtype=complex))
    with pytest.raises(ValueError):
        ConstantWeightCyclicCode().encode(7)
    with pytest.raises(ValueError):
        ConstantWeightCyclicCode().decode([0] * 6)
