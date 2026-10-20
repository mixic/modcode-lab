#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Multipath/Doppler OFDM link with pilot-based channel estimation.

Run from the repository root:

    python examples/multipath_lab.py

The channel is frequency-selective and block-fading: one Rayleigh tap
realisation per OFDM block, with an AR(1) Doppler correlation between blocks.
The receiver estimates the channel from comb pilots (least squares plus linear
interpolation), equalises with zero forcing, and uses the per-subcarrier
post-equalisation variance for the LLRs.  The "genie" curves use the true gains
to separate channel estimation loss from the fading penalty.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from modcode.channel.multipath import TapProcess, subcarrier_gains
from modcode.metrics.theory import theoretical_ber, theoretical_rayleigh_ber
from modcode.modulation.qpsk import QPSK
from modcode.pipeline.ofdm import OFDMLink, run_ofdm_link
from modcode.viz.ber_plot import plot_ber

plt.switch_backend("Agg")

OUTPUT = Path(__file__).resolve().parent / "figures" / "multipath.png"
EBN0_DB = np.asarray([0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0])
N_BITS = 300_000
FLOOR = 1e-6
N_FFT = 64
PILOT_SPACING = 4
TWO_TAP = ((0.7, 0.3), (0, 1))
FOUR_TAP = ((0.5, 0.25, 0.15, 0.1), (0, 1, 3, 5))


def _process(powers, *, seed: int) -> TapProcess:
    return TapProcess(
        powers, doppler_hz=1.0, symbol_rate_hz=1.0, block_symbols=N_FFT, seed=seed
    )


def _link(powers, delays, estimation: str, *, seed: int, spacing: int = PILOT_SPACING):
    return OFDMLink(
        QPSK(),
        _process(powers, seed=seed),
        delays,
        n_fft=N_FFT,
        pilot_spacing=spacing,
        estimation=estimation,
    )


def _sweep(powers, delays, estimation: str) -> np.ndarray:
    values = []
    for index, ebn0_db in enumerate(EBN0_DB):
        result = run_ofdm_link(
            _link(powers, delays, estimation, seed=11 + index),
            n_bits=N_BITS,
            ebn0_db=float(ebn0_db),
            rng=np.random.default_rng(2027 + index),
        )
        values.append(result.ber)
    return np.asarray(values)


def _measure(powers, delays, estimation: str, ebn0_db: float, seed: int):
    return run_ofdm_link(
        _link(powers, delays, estimation, seed=seed),
        n_bits=N_BITS,
        ebn0_db=ebn0_db,
        rng=np.random.default_rng(seed),
    )


def main() -> None:
    flat_genie = _sweep((1.0,), (0,), "genie")
    two_genie = _sweep(TWO_TAP[0], TWO_TAP[1], "genie")
    two_ls = _sweep(TWO_TAP[0], TWO_TAP[1], "ls")
    four_ls = _sweep(FOUR_TAP[0], FOUR_TAP[1], "ls")

    reference = OFDMLink(QPSK(), _process((1.0,), seed=1), (0,))
    print(
        f"OFDM: n_fft={N_FFT}, pilots every {PILOT_SPACING} subcarriers, "
        f"effective rate {reference.effective_rate:.3f}"
    )
    print(f"{N_BITS} bits per point, fast block fading (independent blocks)\n")
    print(
        f"{'Eb/N0':>5}{'flat genie':>13}{'2-tap genie':>13}"
        f"{'2-tap LS':>11}{'4-tap LS':>11}"
    )
    for index, ebn0_db in enumerate(EBN0_DB):
        print(
            f"{ebn0_db:5.1f}{flat_genie[index]:13.2e}{two_genie[index]:13.2e}"
            f"{two_ls[index]:11.2e}{four_ls[index]:11.2e}"
        )

    print("\nchannel-estimation loss at 20 dB (2 taps)")
    for spacing in (2, 4, 8):
        result = run_ofdm_link(
            _link(TWO_TAP[0], TWO_TAP[1], "ls", seed=31, spacing=spacing),
            n_bits=N_BITS,
            ebn0_db=20.0,
            rng=np.random.default_rng(31),
        )
        rate = OFDMLink(
            QPSK(), _process((1.0,), seed=1), (0,), pilot_spacing=spacing
        ).effective_rate
        print(
            f"  pilot spacing {spacing}: BER {result.ber:.2e}, "
            f"NMSE {result.channel_nmse:.2e}, effective rate {rate:.3f}"
        )
    genie = _measure(TWO_TAP[0], TWO_TAP[1], "genie", 20.0, 32)
    print(f"  genie CSI        : BER {genie.ber:.2e}, NMSE {genie.channel_nmse:.2e}")

    print("\nDoppler correlation across one OFDM block (64 symbols)")
    for doppler in (0.0, 0.002, 0.005, 0.02, 1.0):
        process = TapProcess(
            (1.0,),
            doppler_hz=doppler,
            symbol_rate_hz=1.0,
            block_symbols=N_FFT,
            seed=0,
        )
        print(
            f"  f_d = {doppler:6.3f}: rho = {process.correlation:.3f} "
            f"(argument {process.normalized_doppler * N_FFT:.3f})"
        )

    figure, (ber_axis, gain_axis) = plt.subplots(1, 2, figsize=(12.6, 4.9))
    plot_ber(
        EBN0_DB,
        {
            "OFDM flat, genie": np.maximum(flat_genie, FLOOR),
            "OFDM 2-tap, genie": np.maximum(two_genie, FLOOR),
            "OFDM 2-tap, LS": np.maximum(two_ls, FLOOR),
            "OFDM 4-tap, LS": np.maximum(four_ls, FLOOR),
        },
        ber_axis,
        theory={
            "uncoded QPSK, AWGN": theoretical_ber(QPSK(), EBN0_DB),
            "flat Rayleigh, pilot-adjusted": theoretical_rayleigh_ber(
                QPSK(), EBN0_DB + 10.0 * np.log10(0.75)
            ),
        },
    )
    ber_axis.set_ylim(FLOOR / 2.0, 1.0)
    ber_axis.set_title("OFDM over multipath with estimated CSI")
    ber_axis.legend(fontsize=8)

    result = _measure(TWO_TAP[0], TWO_TAP[1], "ls", 10.0, 41)
    four_taps = np.asarray(FOUR_TAP[0], dtype=complex) ** 0.5 * np.exp(
        1j * np.asarray([0.3, 1.1, 2.0, 2.6])
    )
    four_gains = subcarrier_gains(four_taps, FOUR_TAP[1], N_FFT)
    subcarriers = np.arange(N_FFT)
    gain_axis.plot(
        subcarriers, np.abs(result.gains[0]), label="2-tap, true |H| (one block)"
    )
    gain_axis.plot(
        subcarriers,
        np.abs(result.estimates[0]),
        marker="o",
        markersize=3,
        linewidth=1.0,
        label="2-tap, LS estimate (10 dB)",
    )
    gain_axis.plot(subcarriers, np.abs(four_gains), label="4-tap, true |H|")
    gain_axis.set_xlabel("subcarrier")
    gain_axis.set_ylabel("|H(f)|")
    gain_axis.set_title("Frequency selectivity and estimation")
    gain_axis.grid(True, alpha=0.3)
    gain_axis.legend(fontsize=8)

    figure.tight_layout()
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(OUTPUT, dpi=140)
    plt.close(figure)
    print(f"\nfigure -> {OUTPUT}")


if __name__ == "__main__":
    main()
