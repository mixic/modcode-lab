#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Visualise CRC detection and the circulating blocks of a toy QC-MDPC code.

Run from the repository root:

    python examples/cyclic_codes_pqc_lab.py

Writes three figures to ``examples/figures/`` and prints a numeric report.
The QC-MDPC section is a structural demonstration only: the parameters are
tiny, and the module is not cryptography.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from modcode.fec.crc import CRC8
from modcode.pqc import (
    QCMDPC,
    decoder_success_rates,
    dual_weight_distribution,
    syndrome_search_space_bits,
)
from modcode.viz.crc import (
    plot_crc_division,
    plot_crc_frame,
    plot_crc_syndromes,
)
from modcode.viz.pqc import (
    plot_bit_vector,
    plot_circulant_block,
    plot_decoder_success,
    plot_dual_weight_histogram,
    plot_key_scaling,
    plot_parity_check,
    plot_search_space,
    plot_syndrome_history,
)

plt.switch_backend("Agg")

OUTPUT_DIR = Path(__file__).resolve().parent / "figures"
RULE = "=" * 78


def _save(
    figure: plt.Figure,
    path: Path,
    *,
    rect: tuple[float, ...] | None = None,
) -> None:
    figure.tight_layout(rect=rect)
    figure.savefig(path, dpi=140)
    plt.close(figure)
    print(f"  figure -> {path}")


def _crc_figure() -> None:
    message = np.asarray(
        [1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1],
        dtype=np.uint8,
    )
    parity = CRC8.encode(message)[-CRC8.width :]
    codeword = CRC8.encode(message)
    corrupted = codeword.copy()
    corrupted[3] ^= 1
    syndromes = []
    for position in range(codeword.size):
        candidate = codeword.copy()
        candidate[position] ^= 1
        syndromes.append(CRC8.remainder(candidate))

    figure, axes = plt.subplots(1, 3, figsize=(16.0, 4.2))
    figure.suptitle(
        "CRC-8: systematic cyclic encoding, shift-register division, and detection",
        fontsize=12,
    )
    plot_crc_frame(
        axes[0],
        message,
        parity,
        error_positions=[3],
        title="codeword with one bit error",
    )
    plot_crc_division(
        axes[1],
        CRC8.division_trace(message),
        title="remainder register after each message bit",
    )
    plot_crc_syndromes(
        axes[2],
        range(codeword.size),
        syndromes,
        title="single-bit errors all leave a nonzero remainder",
    )
    _save(figure, OUTPUT_DIR / "crc_encoding.png", rect=(0, 0, 1, 0.94))

    print()
    print(RULE)
    print("CRC-8 DETECTION")
    print(RULE)
    print(f"message bits          : {message.size}")
    print(f"CRC width             : {CRC8.width}")
    print(f"codeword bits         : {codeword.size}")
    print(f"generator             : 0x{CRC8.polynomial:02X} (x^8 + x^2 + x + 1)")
    print("corrupted position    : 3")
    print(f"syndrome after error  : 0x{CRC8.remainder(corrupted):02X}")
    print(f"detected              : {not CRC8.check(corrupted)}")
    print("  Every single-bit error in the 24-bit codeword leaves a nonzero")
    print("  remainder; the right-hand panel plots those 24 syndromes.")


def _qc_figure() -> None:
    code = QCMDPC(block_length=47, row_weight=3, seed=7)
    rng = np.random.default_rng(2026)
    message = rng.integers(0, 2, size=code.dimension, dtype=np.uint8)
    ciphertext, error = code.encrypt(message, error_weight=2, seed=11)
    result = code.decode(ciphertext, error_weight=2)
    syndrome_before = code.syndrome(ciphertext)

    figure, axes = plt.subplots(1, 3, figsize=(16.0, 5.0))
    figure.suptitle(
        "QC-MDPC McEliece-style toy: one sparse row defines each circulant block",
        fontsize=12,
    )
    plot_parity_check(axes[0], code, title="H = [H0 | H1]")
    plot_circulant_block(axes[1], code.h0, title="H0")
    sizes = np.asarray([31, 47, 63, 95, 127])
    plot_key_scaling(
        axes[2],
        sizes,
        sizes.astype(float) ** 2,
        sizes.astype(float),
        title="circulant Q: n bits, not n^2",
    )
    _save(figure, OUTPUT_DIR / "qc_mdpc_matrix.png", rect=(0, 0, 1, 0.93))

    figure, axes = plt.subplots(1, 3, figsize=(16.0, 4.2))
    figure.suptitle(
        "Toy QC-MDPC decoding: low-weight error, syndrome, and bit-flipping progress",
        fontsize=12,
    )
    plot_bit_vector(axes[0], error, title="injected error (weight 2)", color="tab:red")
    plot_bit_vector(
        axes[1],
        syndrome_before,
        title="syndrome before decoding",
        color="tab:purple",
    )
    plot_syndrome_history(axes[2], result.syndrome_history, title="syndrome weight")
    _save(figure, OUTPUT_DIR / "qc_mdpc_decoding.png", rect=(0, 0, 1, 0.93))

    sizes = code.key_size_bits
    print()
    print(RULE)
    print("TOY QC-MDPC MCELIECE STRUCTURE (NOT SECURE)")
    print(RULE)
    print(f"block length n        : {code.block_length}")
    print(f"code length 2n        : {code.code_length}")
    print(f"dimension k           : {code.dimension}")
    print(f"rate                  : {code.rate:.3f}")
    print(f"row weight per block  : {code.row_weight}")
    print(f"column weight of H    : {2 * code.row_weight}")
    print(f"dense Q bits          : {sizes['dense_generator']:,}")
    print(f"circulant Q bits      : {sizes['qc_generator']:,}")
    print(f"compression ratio     : {code.compression_ratio:.0f}x")
    print(f"injected error weight : {int(error.sum())}")
    print(f"decoder iterations    : {result.iterations}")
    print(f"exhaustive fallback   : {result.used_search}")
    print(f"message recovered     : {np.array_equal(result.message, message)}")
    print()
    print("  The same circulant-block idea appears in QC-LDPC FEC. Here it is")
    print("  used for key-size reduction in a McEliece-style construction. The")
    print("  security assumption is the hardness of decoding a random")
    print("  quasi-cyclic code; the toy parameters and decoder are not secure.")


def _qc_security_figure() -> None:
    code = QCMDPC(block_length=47, row_weight=3, seed=7)
    weights, success = decoder_success_rates(
        code, np.arange(1, 7), trials=15, seed=42
    )
    search_weights = np.arange(1, 16)
    search_bits = np.asarray(
        [
            syndrome_search_space_bits(code.code_length, int(weight))
            for weight in search_weights
        ]
    )

    small = QCMDPC(block_length=15, row_weight=3, seed=3)
    counts = dual_weight_distribution(small)
    minimum = int(np.flatnonzero(counts[1:])[0] + 1)
    secret_weight = 2 * small.row_weight

    figure, axes = plt.subplots(1, 3, figsize=(16.0, 4.4))
    figure.suptitle(
        "QC-MDPC security experiments: decoder threshold, search space, "
        "and low-weight dual words",
        fontsize=12,
    )
    plot_decoder_success(
        axes[0],
        weights,
        success,
        title="bit-flipping success (n=47)",
    )
    plot_search_space(
        axes[1],
        search_weights,
        search_bits,
        title="log2 error patterns (2n=94)",
    )
    plot_dual_weight_histogram(
        axes[2],
        counts,
        secret_weight=secret_weight,
        title="dual weights (n=15 toy)",
    )
    _save(figure, OUTPUT_DIR / "qc_mdpc_security.png", rect=(0, 0, 1, 0.92))

    print()
    print(RULE)
    print("QC-MDPC SECURITY EXPERIMENTS (TOY)")
    print(RULE)
    summary = ", ".join(
        f"{int(weight)}:{rate:.2f}"
        for weight, rate in zip(weights, success, strict=True)
    )
    print(f"decoder success by weight : {summary}")
    print(f"log2 patterns at t=2      : {search_bits[1]:.1f} bits")
    print(f"log2 patterns at t=6      : {search_bits[5]:.1f} bits")
    print(f"small toy minimum dual wt : {minimum}")
    print(f"secret row weight         : {secret_weight}")
    print()
    print("  The success curve falls as the error weight crosses the decoder's")
    print("  correction capability. The search-space panel shows why the")
    print("  syndrome-decoding problem does not become easy at toy scale. The")
    print("  right panel exposes the structural attack surface: the dual code")
    print("  contains the sparse secret rows, so key recovery is a low-weight")
    print("  codeword search. Real parameters must make that search infeasible.")


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    print(RULE)
    print("CYCLIC CODES + POST-QUANTUM STRUCTURE LAB")
    print(f"Python {sys.version_info.major}.{sys.version_info.minor}."
          f"{sys.version_info.micro} | numpy {np.__version__}")
    print(RULE)
    _crc_figure()
    _qc_figure()
    _qc_security_figure()
    print()
    print(RULE)
    print(f"done -- figures in {OUTPUT_DIR}")
    print(RULE)


if __name__ == "__main__":
    main()
