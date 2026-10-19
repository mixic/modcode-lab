# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Inspect a lifted QC-LDPC matrix and Polar-code BEC channel polarization."""

import numpy as np

from modcode.fec.polar import PolarBEC, bec_capacity
from modcode.fec.qc_ldpc import expand_qc_base_matrix


def main() -> None:
    base = np.array([[0, 1, -1], [2, -1, 0]], dtype=np.int64)
    parity_check = expand_qc_base_matrix(base, lifting_factor=3)
    print("QC-LDPC base matrix (-1 means a zero block):")
    print(base)
    print("Expanded binary parity-check matrix:")
    print(parity_check)
    print("Column sums:", parity_check.sum(axis=0))

    erasure_probability = 0.4
    code = PolarBEC(
        block_length=16,
        information_bits=9,
        erasure_probability=erasure_probability,
    )
    probabilities = code.synthetic_erasure_probabilities
    print(f"\nBEC capacity: {bec_capacity(erasure_probability):.3f} bits/use")
    print(f"Selected Polar rate: {code.rate:.3f}")
    print("Synthetic-channel erasure probabilities:")
    print(np.round(probabilities, 4))
    print("Selected information-channel indices:", code.information_indices)
    print(
        "Mean synthetic-channel capacity:",
        f"{np.mean(1.0 - probabilities):.3f} bits/use",
    )

    rng = np.random.default_rng(2026)
    message = rng.integers(0, 2, code.information_bits, dtype=np.uint8)
    received = code.encode(message).astype(np.int8)
    received[rng.random(code.block_length) < erasure_probability] = -1
    decoded = code.decode_bec(received)
    print("\nOne reproducible BEC decoding trial:")
    print("Information bits:", message)
    print("Received bits (-1 means erased):", received)
    print("SC decoded bits:", decoded)
    print("Decoded correctly:", np.array_equal(message, decoded))


if __name__ == "__main__":
    main()
