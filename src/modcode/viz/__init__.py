# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Plotting helpers for modulation and coding experiments."""

from modcode.viz.ber_plot import plot_ber
from modcode.viz.constellation import plot_constellation
from modcode.viz.crc import (
    plot_crc_division,
    plot_crc_frame,
    plot_crc_syndromes,
)
from modcode.viz.number_theory import (
    plot_residual_index_convergence,
    plot_residual_index_density,
    plot_residual_index_histogram,
)
from modcode.viz.phasedarray import (
    plot_pattern,
    plot_pattern_family,
    plot_phase_ramp,
    plot_squint,
    plot_wavefront,
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
from modcode.viz.qfc import (
    plot_noise_budget,
    plot_squeezing_survival,
    plot_thermal_occupation,
    plot_transport_budget,
    plot_wigner,
)

__all__ = [
    "plot_ber",
    "plot_bit_vector",
    "plot_circulant_block",
    "plot_constellation",
    "plot_crc_division",
    "plot_crc_frame",
    "plot_crc_syndromes",
    "plot_decoder_success",
    "plot_dual_weight_histogram",
    "plot_key_scaling",
    "plot_noise_budget",
    "plot_parity_check",
    "plot_pattern",
    "plot_pattern_family",
    "plot_phase_ramp",
    "plot_residual_index_convergence",
    "plot_residual_index_density",
    "plot_residual_index_histogram",
    "plot_search_space",
    "plot_squeezing_survival",
    "plot_squint",
    "plot_syndrome_history",
    "plot_thermal_occupation",
    "plot_transport_budget",
    "plot_wavefront",
    "plot_wigner",
]
