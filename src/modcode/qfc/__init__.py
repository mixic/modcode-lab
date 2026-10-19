# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Photon budgets, transduction noise, and transport loss for QFC links."""

from modcode.qfc.link import (
    TransportPath,
    aperture_gain_db,
    coax_loss,
    fiber_loss,
    free_space_link_loss,
    free_space_path_loss,
    pointing_loss_db,
    transport_budget,
)
from modcode.qfc.photons import (
    BOLTZMANN,
    PLANCK,
    SPEED_OF_LIGHT,
    Transducer,
    covariance_squeezing_db,
    photon_energy,
    squeezed_vacuum_covariance,
    thermal_occupation,
    thermal_squeezed_covariance,
    wigner_function,
)

__all__ = [
    "BOLTZMANN",
    "PLANCK",
    "SPEED_OF_LIGHT",
    "Transducer",
    "TransportPath",
    "aperture_gain_db",
    "coax_loss",
    "covariance_squeezing_db",
    "fiber_loss",
    "free_space_link_loss",
    "free_space_path_loss",
    "photon_energy",
    "pointing_loss_db",
    "squeezed_vacuum_covariance",
    "thermal_occupation",
    "thermal_squeezed_covariance",
    "transport_budget",
    "wigner_function",
]
