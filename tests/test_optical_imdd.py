# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Tests for optical intensity modulation and the candidate design search."""

import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg")

from modcode.optical import (
    OOK,
    IMDDChannel,
    IntensityPAM,
    NonUniformPAM,
    ber_sweep_imdd,
    optimize_pam_levels,
)
from modcode.viz.optical import (
    plot_intensity_constellations,
    plot_received_levels,
    plot_robustness,
)


def test_ook_and_pam_levels_are_nonnegative_with_unit_average() -> None:
    assert np.array_equal(OOK().constellation, [0.0, 2.0])
    modulation = IntensityPAM(4)
    assert np.all(modulation.constellation >= 0)
    assert np.isclose(modulation.constellation.mean(), 1.0)


def test_intensity_pam_round_trip() -> None:
    modulation = IntensityPAM(4)
    bits = np.array([0, 0, 0, 1, 1, 0, 1, 1], dtype=np.uint8)
    assert np.array_equal(modulation.demodulate(modulation.modulate(bits)), bits)


def test_neighboring_pam_levels_have_gray_labels() -> None:
    modulation = IntensityPAM(4)
    decoded = modulation.demodulate(modulation.constellation)
    labels = decoded.reshape(-1, modulation.bits_per_symbol)
    assert all(
        np.count_nonzero(first != second) == 1
        for first, second in zip(labels, labels[1:], strict=False)
    )


def test_nonuniform_levels_require_physical_ordered_unit_mean_levels() -> None:
    with pytest.raises(ValueError):
        NonUniformPAM(np.array([0.0, 0.5, 0.4, 3.1]))
    with pytest.raises(ValueError):
        NonUniformPAM(np.array([-0.1, 0.5, 1.5, 2.1]))


def test_optimizer_returns_ordered_unit_mean_candidate() -> None:
    modulation = optimize_pam_levels(
        noise_std=0.05,
        amplitude_noise_std=0.1,
    )
    assert modulation.order == 4
    assert modulation.levels[0] == 0.0
    assert modulation.levels[-1] == 2.0
    assert np.all(np.diff(modulation.levels) > 0)
    assert np.isclose(modulation.levels.mean(), 1.0)


def test_noiseless_channel_preserves_intensity_and_seed_is_reproducible() -> None:
    symbols = np.array([0.0, 2 / 3, 4 / 3, 2.0])
    channel = IMDDChannel()
    first = channel.transmit(symbols, rng=np.random.default_rng(3))
    second = channel.transmit(symbols, rng=np.random.default_rng(3))
    assert np.array_equal(first, symbols)
    assert np.array_equal(first, second)


def test_channel_applies_noise_amplitude_nonlinearity_and_bandwidth_proxy() -> None:
    symbols = np.tile(np.array([0.0, 2.0]), 1000)
    channel = IMDDChannel(
        noise_std=0.02,
        amplitude_noise_std=0.1,
        nonlinearity=0.1,
        bandwidth=0.1,
    )
    received = channel.transmit(symbols, rng=np.random.default_rng(1))
    assert received.shape == symbols.shape
    assert not np.array_equal(received, symbols)


def test_optical_ber_sweep_compares_candidate_formats_and_conditions() -> None:
    candidate = optimize_pam_levels(noise_std=0.1, amplitude_noise_std=0.05)
    results = ber_sweep_imdd(
        {"OOK": OOK(), "PAM4": IntensityPAM(4), "candidate": candidate},
        {
            "noise": IMDDChannel(noise_std=0.1),
            "combined": IMDDChannel(
                noise_std=0.1,
                amplitude_noise_std=0.05,
                nonlinearity=0.05,
                bandwidth=0.05,
            ),
        },
        n_bits=20_000,
        seed=4,
    )
    assert set(results) == {"OOK", "PAM4", "candidate"}
    assert all(
        0.0 <= float(ber[0]) <= 1.0
        for conditions in results.values()
        for ber in conditions.values()
    )


def test_optical_plots_return_matplotlib_figures() -> None:
    modulation = IntensityPAM(4)
    results = ber_sweep_imdd(
        {"PAM4": modulation},
        {"noise": IMDDChannel(noise_std=0.1)},
        n_bits=1000,
        seed=5,
    )
    figures = [
        plot_intensity_constellations({"OOK": OOK(), "PAM4": modulation})[0],
        plot_received_levels(modulation, IMDDChannel(noise_std=0.1))[0],
        plot_robustness(results)[0],
    ]
    for figure in figures:
        figure.canvas.draw()
        matplotlib.pyplot.close(figure)
