# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Tests for soft-decision LLRs, Hamming soft decoding, and BER theory."""

import numpy as np
import pytest

from modcode.channel.awgn import AWGNChannel
from modcode.fec.hamming import Hamming74
from modcode.metrics.theory import q_function, theoretical_ber
from modcode.modulation.apsk import APSK
from modcode.modulation.bpsk import BPSK
from modcode.modulation.qam import QAM
from modcode.modulation.qpsk import QPSK
from modcode.pipeline.chain import Link, run_link


def test_bpsk_llr_sign_and_scale() -> None:
    symbols = np.asarray([1.0 + 0.0j, -1.0 + 0.0j])
    llr = BPSK().soft_demodulate(symbols, noise_variance=0.5)
    np.testing.assert_allclose(llr, [4.0, -4.0])


def test_qpsk_llr_separates_i_and_q() -> None:
    qpsk = QPSK()
    point = qpsk.constellation[(0, 1)]
    llr = qpsk.soft_demodulate(np.asarray([point]), noise_variance=0.5)
    assert llr.shape == (2,)
    assert llr[0] > 0.0
    assert llr[1] < 0.0


def test_empty_soft_demodulation_returns_empty() -> None:
    assert BPSK().soft_demodulate(np.empty(0, dtype=complex), 1.0).size == 0


def test_hamming_soft_decoding_corrects_every_single_bit_error() -> None:
    code = Hamming74()
    rng = np.random.default_rng(5)
    data = rng.integers(0, 2, 16 * 4, dtype=np.uint8)
    codeword = code.encode(data)
    for position in range(codeword.size):
        corrupted = codeword.copy()
        corrupted[position] ^= 1
        llr = 1.0 - 2.0 * corrupted.astype(np.float64)
        np.testing.assert_array_equal(code.decode_soft(llr), data)


def test_hamming_soft_decoding_uses_unreliable_bits() -> None:
    code = Hamming74()
    data = np.asarray([1, 0, 1, 1], dtype=np.uint8)
    codeword = code.encode(data)
    llr = 1.0 - 2.0 * codeword.astype(np.float64)
    llr[0] = -0.05
    llr[1] = 8.0
    np.testing.assert_array_equal(code.decode_soft(llr), data)


def test_link_reports_soft_decoding_for_hamming() -> None:
    result = run_link(
        Link(BPSK(), AWGNChannel(), Hamming74()),
        n_bits=4_000,
        ebn0_db=8.0,
        rng=np.random.default_rng(3),
    )
    assert result.decoding == "soft"
    assert result.ber < 0.02


def test_link_hard_decoding_returns_data_bits() -> None:
    result = run_link(
        Link(BPSK(), AWGNChannel(), Hamming74(), decoding="hard"),
        n_bits=4_000,
        ebn0_db=8.0,
        rng=np.random.default_rng(3),
    )
    assert result.decoding == "hard"
    assert result.received_bits.size == result.source_bits.size
    assert result.ber < 0.02


def test_link_rejects_unknown_decoding_mode() -> None:
    with pytest.raises(ValueError, match="decoding"):
        Link(BPSK(), AWGNChannel(), Hamming74(), decoding="turbo")


def test_theoretical_ber_for_bpsk_and_qpsk() -> None:
    value = theoretical_ber(BPSK(), [0.0])
    assert np.isclose(value[0], 0.07864960352514257, rtol=1e-12)
    np.testing.assert_allclose(theoretical_ber(QPSK(), [0.0]), value)


def test_theoretical_ber_for_16qam_is_monotone() -> None:
    values = theoretical_ber(QAM(16), [0.0, 4.0, 8.0, 12.0])
    assert np.all(np.diff(values) < 0.0)
    assert np.all((values > 0.0) & (values < 1.0))


def test_q_function_values() -> None:
    np.testing.assert_allclose(
        q_function([0.0, 1.0, 2.0]),
        [0.5, 0.15865525393145707, 0.02275013194817921],
        rtol=1e-12,
    )


def test_theory_rejects_unsupported_modulation() -> None:
    with pytest.raises(ValueError, match="no theoretical"):
        theoretical_ber(APSK(), [0.0])
