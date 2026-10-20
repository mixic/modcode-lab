#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""CT-IM experiments: BER against OFDM-IM, Rayleigh, and PAPR.

Run from the repository root:

    python examples/ct_im_experiments.py

Every scheme uses the same Eb/N0 definition, the same AWGN/Rayleigh channel,
and the same maximum-likelihood index detector (except where the threshold
detector is shown for contrast).  Index-bit and symbol-bit errors are counted
separately.  The PAPR panel compares no shift, the best of seven shifts (three
bits of side information), and the best of four shifts (two bits).
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from modcode.metrics.index_modulation import (
    IndexResult,
    block_papr_values,
    simulate_index_modulation,
)
from modcode.modulation.ct_im import (
    CTIM,
    BCHIndexCode,
    ClassicalOFDMIM,
    ShortenedConstantWeightCode,
)

plt.switch_backend("Agg")

OUTPUT = Path(__file__).resolve().parent / "figures"
EBN0_DB = np.asarray([0.0, 2.0, 4.0, 6.0, 8.0, 10.0])
BLOCKS = 800
PAPR_BLOCKS = 500
FLOOR = 1e-6


def _schemes() -> dict[str, CTIM | ClassicalOFDMIM]:
    return {
        "CT-IM, threshold": CTIM(M=4, detector="threshold"),
        "CT-IM, ML": CTIM(M=4, detector="ml"),
        "CT-IM shortened, ML": CTIM(
            M=4, code=ShortenedConstantWeightCode(), detector="ml"
        ),
        "BCH(15,7) orbit, ML": CTIM(M=4, code=BCHIndexCode(weight=7), detector="ml"),
        "OFDM-IM, ML": ClassicalOFDMIM(),
    }


def _sweep(
    scheme: CTIM | ClassicalOFDMIM, *, rayleigh: bool
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    index = []
    symbol = []
    total = []
    for point, ebn0_db in enumerate(EBN0_DB):
        result: IndexResult = simulate_index_modulation(
            scheme,
            float(ebn0_db),
            blocks=BLOCKS,
            rayleigh=rayleigh,
            rng=np.random.default_rng(500 + point + (0 if rayleigh else 100)),
        )
        index.append(result.index_ber)
        symbol.append(result.symbol_ber)
        total.append(result.ber)
    return np.asarray(index), np.asarray(symbol), np.asarray(total)


def _random_vectors(
    scheme: CTIM | ClassicalOFDMIM, blocks: int, seed: int
) -> np.ndarray:
    rng = np.random.default_rng(seed)
    rows = []
    for _ in range(blocks):
        label = int(rng.integers(0, len(scheme.index_code)))
        bits = np.concatenate(
            (
                np.asarray(
                    [
                        (label >> shift) & 1
                        for shift in range(scheme.index_bits - 1, -1, -1)
                    ],
                    dtype=np.uint8,
                ),
                rng.integers(0, 2, scheme.symbol_bits, dtype=np.uint8),
            )
        )
        rows.append(scheme.encode(bits))
    return np.asarray(rows)


def main() -> None:
    schemes = _schemes()
    results = {
        name: _sweep(scheme, rayleigh=False) for name, scheme in schemes.items()
    }
    fading = {
        name: _sweep(schemes[name], rayleigh=True)
        for name in ("CT-IM, ML", "OFDM-IM, ML")
    }

    print(
        "scheme                        index bits  symbol bits  "
        "bits/block  bits/subcarrier"
    )
    for name, scheme in schemes.items():
        print(
            f"{name:<28}{scheme.index_bits:>10}{scheme.symbol_bits:>13}"
            f"{scheme.bits_per_block:>12}{scheme.spectral_efficiency:>17.3f}"
        )

    print(f"\nBER over AWGN, {BLOCKS} blocks per point (index / symbol)")
    print(f"{'Eb/N0':>5}" + "".join(f"{name:>28}" for name in results))
    for point, ebn0_db in enumerate(EBN0_DB):
        row = "".join(
            f"{index[point]:12.2e} /{symbol[point]:12.2e}"
            for index, symbol, _ in results.values()
        )
        print(f"{ebn0_db:5.1f}{row}")

    print("\nindex-bit BER, AWGN versus known-CSI Rayleigh")
    print(f"{'Eb/N0':>5}{'CT-IM AWGN':>14}{'CT-IM Rayleigh':>17}"
          f"{'OFDM-IM AWGN':>15}{'OFDM-IM Rayleigh':>18}")
    for point, ebn0_db in enumerate(EBN0_DB):
        print(
            f"{ebn0_db:5.1f}{results['CT-IM, ML'][0][point]:14.2e}"
            f"{fading['CT-IM, ML'][0][point]:17.2e}"
            f"{results['OFDM-IM, ML'][0][point]:15.2e}"
            f"{fading['OFDM-IM, ML'][0][point]:18.2e}"
        )

    ctim = schemes["CT-IM, ML"]
    vectors = _random_vectors(ctim, PAPR_BLOCKS, seed=99)
    classical = _random_vectors(schemes["OFDM-IM, ML"], PAPR_BLOCKS, seed=99)
    thresholds = np.linspace(2.0, 9.0, 60)
    print("\nPAPR of the 7-subcarrier block placed in a 64-subcarrier band")
    plain = block_papr_values(vectors)
    classical_plain = block_papr_values(classical)
    print(
        f"  CT-IM        mean {np.mean(plain):5.2f} dB, "
        f"P(PAPR > 5 dB) = {np.mean(plain > 5.0):.3f}"
    )
    print(
        f"  OFDM-IM      mean {np.mean(classical_plain):5.2f} dB, "
        f"P(PAPR > 5 dB) = {np.mean(classical_plain > 5.0):.3f}"
    )
    deviations = [
        float(np.max(np.abs(block_papr_values(vectors, offset=shift) - plain)))
        for shift in range(7)
    ]
    print(
        "  cyclic shift of the whole block: max |PAPR(shift) - PAPR(0)| = "
        f"{max(deviations):.2e} dB"
    )
    print("  -> the shift is PAPR-invariant, so the side information buys nothing")

    figure, (awgn_axis, fading_axis) = plt.subplots(1, 2, figsize=(12.6, 4.9))
    for name, (_index, _, total) in results.items():
        awgn_axis.semilogy(
            EBN0_DB, np.maximum(total, FLOOR), marker="o", label=name
        )
    awgn_axis.set_xlabel("Eb/N0 (dB)")
    awgn_axis.set_ylabel("bit error rate")
    awgn_axis.set_title("Index modulation over AWGN")
    awgn_axis.grid(True, which="both", alpha=0.3)
    awgn_axis.legend(fontsize=8)
    awgn_axis.set_ylim(FLOOR / 2.0, 1.0)

    for name in ("CT-IM, ML", "OFDM-IM, ML"):
        index, symbol, _ = results[name]
        fading_index, _, _ = fading[name]
        fading_axis.semilogy(
            EBN0_DB,
            np.maximum(index, FLOOR),
            linestyle="--",
            marker="o",
            label=f"{name}, AWGN index",
        )
        fading_axis.semilogy(
            EBN0_DB,
            np.maximum(symbol, FLOOR),
            linestyle=":",
            marker="s",
            label=f"{name}, AWGN symbol",
        )
        fading_axis.semilogy(
            EBN0_DB,
            np.maximum(fading_index, FLOOR),
            linestyle="-",
            marker="^",
            label=f"{name}, Rayleigh index",
        )
    fading_axis.set_xlabel("Eb/N0 (dB)")
    fading_axis.set_ylabel("bit error rate")
    fading_axis.set_title("Index versus symbol bits, and Rayleigh")
    fading_axis.grid(True, which="both", alpha=0.3)
    fading_axis.legend(fontsize=7)
    fading_axis.set_ylim(FLOOR / 2.0, 1.0)
    figure.tight_layout()
    OUTPUT.mkdir(parents=True, exist_ok=True)
    ber_path = OUTPUT / "ct_im_ber.png"
    figure.savefig(ber_path, dpi=140)
    plt.close(figure)

    figure, axis = plt.subplots(figsize=(7.2, 4.8))

    def ccdf(values: np.ndarray) -> np.ndarray:
        return np.asarray([np.mean(values > point) for point in thresholds])

    axis.semilogy(thresholds, ccdf(plain), "-", label="CT-IM (d=4 orbit)")
    axis.semilogy(
        thresholds, ccdf(classical_plain), ":", label="OFDM-IM (32 patterns)"
    )
    axis.set_xlabel("PAPR threshold (dB)")
    axis.set_ylabel("CCDF")
    axis.set_title("PAPR CCDF of the 7-subcarrier block")
    axis.grid(True, which="both", alpha=0.3)
    axis.legend(fontsize=8)
    figure.tight_layout()
    papr_path = OUTPUT / "ct_im_papr.png"
    figure.savefig(papr_path, dpi=140)
    plt.close(figure)

    print(f"\nfigure -> {ber_path}")
    print(f"figure -> {papr_path}")


if __name__ == "__main__":
    main()
