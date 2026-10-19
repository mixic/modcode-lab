#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Wavefront synthesis and quantum frequency conversion, end to end.

Run from the repository root:

    python examples/wavefront_qfc_lab.py

The script writes five figures to ``examples/figures/`` and prints a numeric
report. It is the integrated counterpart of the two scratch scripts under
``src/modcode/cyclic-5g/``: the same questions, but built on the importable
``modcode.phasedarray`` and ``modcode.qfc`` modules and checked by the test
suite instead of by prose.

Part A -- wavefronts
    A1  Field synthesis: exact phase ramps, verified against a full
        cylindrical-wave sum rather than the array-factor formula.
    A2  Cyclic phase resolution: per-element rounding versus quantising the
        inter-element step once, at the awkward angle of 17 degrees.
    A3  Beamwidth and the grating-lobe spacing limit.
    A4  Wideband squint: a phase-steered beam moves towards broadside as the
        drive frequency rises.

Part B -- quantum frequency conversion
    B1  What 100 km of coax, waveguide, free space, and fibre cost.
    B2  What conversion costs a quantum state: efficiency versus noise.
    B3  The thermal floor at the microwave end.
    B4  Whether a squeezed state survives the converter and the bath.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from modcode.phasedarray import (
    SPEED_OF_LIGHT,
    LinearArray,
    grating_lobe_limit,
    squint_angle,
)
from modcode.qfc import (
    Transducer,
    covariance_squeezing_db,
    free_space_link_loss,
    photon_energy,
    squeezed_vacuum_covariance,
    thermal_occupation,
    transport_budget,
    wigner_function,
)
from modcode.viz.phasedarray import (
    plot_pattern_family,
    plot_phase_ramp,
    plot_squint,
    plot_wavefront,
)
from modcode.viz.qfc import (
    plot_noise_budget,
    plot_squeezing_survival,
    plot_thermal_occupation,
    plot_transport_budget,
    plot_wigner,
)

plt.switch_backend("Agg")

RF_FREQUENCY = 10.0e9
OPTICAL_WAVELENGTH = 1550e-9
OUTPUT_DIR = Path(__file__).resolve().parent / "figures"
RULE = "=" * 78


def _section(title: str) -> None:
    print()
    print(RULE)
    print(title)
    print(RULE)


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


def _wavefront_figure(array: LinearArray) -> None:
    angles = (0.0, 15.0, 30.0, 45.0)
    x = np.linspace(-0.36, 0.36, 520)
    z = np.linspace(0.004, 0.72, 580)
    figure, axes = plt.subplots(2, 4, figsize=(17.5, 8.2))
    figure.suptitle(
        "Wavefront shaping by cyclic phase ramps   "
        f"(f={array.frequency / 1e9:.0f} GHz, lambda={array.wavelength * 100:.1f} cm, "
        f"N={array.elements}, d=lambda/2, aperture={array.aperture * 100:.1f} cm)",
        fontsize=13,
    )

    _section("A1. FIELD SYNTHESIS: full-wave peak versus commanded angle")
    print(
        f"{'commanded':>10}{'dphi':>9}{'array factor':>14}"
        f"{'field sum':>12}{'error':>9}"
    )
    print("-" * 78)
    for column, theta in enumerate(angles):
        phases = array.exact_phases(theta)
        peak = array.far_field_peak(phases)
        print(
            f"{theta:>9.0f}deg"
            f"{np.degrees(array.steering_step(theta)):>8.1f}deg"
            f"{array.peak_angle(phases):>13.2f}deg"
            f"{peak:>11.2f}deg"
            f"{peak - theta:>+8.2f}"
        )
        plot_wavefront(
            axes[0, column],
            array,
            x,
            z,
            phases,
            theta,
        )
        plot_wavefront(
            axes[1, column],
            array,
            x,
            z,
            phases,
            theta,
            component="phase",
            title=f"phase, theta0={theta:.0f} deg",
        )
        axes[0, column].set_title(
            f"theta0={theta:.0f} deg\n"
            f"dphi={np.degrees(array.steering_step(theta)):.1f} deg",
            fontsize=9,
        )
        if column == 0:
            axes[0, column].set_ylabel("range z (m)", fontsize=8)
            axes[1, column].set_ylabel("range z (m)", fontsize=8)

    figure.text(
        0.5,
        0.012,
        "Dashed line = commanded direction. The fringes are wavefronts and stay "
        f"curved out to the Fresnel distance L^2/lambda = "
        f"{array.fresnel_distance:.2f} m.",
        ha="center",
        fontsize=9,
    )
    _save(figure, OUTPUT_DIR / "wavefront_steering.png", rect=(0, 0.03, 1, 0.95))
    print()
    print("  The field sum adds real cylindrical waves on a 10 m circle and")
    print("  agrees with the array factor to the 0.01 deg sampling step.")


def _quantization_figure(array: LinearArray) -> None:
    theta0 = 17.0
    variants = [("exact", None), ("C4", 4), ("C8", 8), ("C16", 16)]
    x = np.linspace(-0.36, 0.36, 520)
    z = np.linspace(0.004, 0.72, 580)
    figure, axes = plt.subplots(3, 4, figsize=(17.5, 10.4))
    figure.suptitle(
        f"Commanded to {theta0:.0f} deg: two ways to spend a cyclic phase budget "
        f"(ideal step {np.degrees(array.steering_step(theta0)):.1f} deg)",
        fontsize=13,
    )

    _section("A2. CYCLIC PHASE RESOLUTION: per-element rounding versus step rounding")
    print(
        f"{'variant':<8}{'quantum':>9}{'elem peak':>11}{'loss@cmd':>10}"
        f"{'rms(rad)':>10}{'strehl':>9}{'step peak':>11}{'step@cmd':>10}"
    )
    print("-" * 78)

    for column, (label, order) in enumerate(variants):
        if order is None:
            phases = array.exact_phases(theta0)
            step_phases = phases
            peak = array.peak_angle(phases)
            loss = array.scan_loss_db(theta0, phases)
            rms = array.aberration_rms(phases, theta0)
            strehl_db = 10.0 * np.log10(array.strehl_estimate(phases, theta0))
            step_peak_text = "--"
            step_loss_text = "--"
            title = "exact ramp"
        else:
            phases = array.element_quantized_phases(theta0, order)
            step_phases = array.step_quantized_phases(theta0, order)
            peak = array.peak_angle(phases)
            loss = array.scan_loss_db(theta0, phases)
            rms = array.aberration_rms(phases, theta0)
            strehl_db = 10.0 * np.log10(array.strehl_estimate(phases, theta0))
            step_peak = array.peak_angle(step_phases)
            step_loss = array.scan_loss_db(theta0, step_phases)
            step_peak_text = f"{step_peak:>10.2f}"
            step_loss_text = f"{step_loss:>+9.2f}"
            title = f"{label}, quantum {360.0 / order:.1f} deg"

        quantum_text = "--" if order is None else f"{360.0 / order:.1f}"
        print(
            f"{label:<8}{quantum_text:>9}{peak:>11.2f}{loss:>+10.2f}"
            f"{rms:>10.3f}{strehl_db:>+9.2f}{step_peak_text:>11}{step_loss_text:>10}"
        )

        plot_phase_ramp(
            axes[0, column],
            array,
            theta0,
            phases,
            label=label,
            title=f"applied phase: {title}",
        )
        plot_wavefront(
            axes[1, column],
            array,
            x,
            z,
            phases,
            theta0,
            title=f"Re E: {title}",
        )
        patterns = {label: phases}
        if order is not None:
            patterns["uniform step"] = step_phases
        pattern_title = (
            "pattern: exact" if order is None else "pattern: element vs uniform step"
        )
        plot_pattern_family(
            axes[2, column],
            array,
            patterns,
            commanded_deg=theta0,
            floor_db=-45.0,
            title=pattern_title,
        )
        if column == 0:
            axes[1, column].set_ylabel("range z (m)", fontsize=8)
            axes[2, column].set_ylabel("array factor (dB)", fontsize=8)

    figure.text(
        0.5,
        0.012,
        "Per-element rounding keeps the beam near 17 deg and turns the excess into "
        "wavefront ripple. Rounding the step once keeps the wavefront flat but aims "
        "at the quantised angle (30 deg for C4, 14.5 deg for C8/C16).",
        ha="center",
        fontsize=9,
    )
    _save(figure, OUTPUT_DIR / "phase_quantization.png", rect=(0, 0.03, 1, 0.955))
    print()
    print("  Element rounding: the average ramp slope survives, so pointing")
    print("  barely moves; the staircase is an aberration. The Marechal")
    print("  estimate exp(-rms^2) predicts the scan loss to within 0.02 dB.")
    print("  Uniform-step rounding: the wavefront is perfect, so the beam has")
    print("  full gain -- aimed at the quantised angle, not the commanded one.")


def _beam_report(array: LinearArray) -> None:
    _section("A3. BEAMWIDTH AND GRATING LOBES")
    big = LinearArray(elements=100, spacing_over_lambda=0.5, frequency=RF_FREQUENCY)
    width = big.half_power_width_deg(30.0, big.exact_phases(30.0), span_deg=5.0)
    theory = np.degrees(0.886 / (big.elements * big.spacing_over_lambda))
    theory_scanned = theory / np.cos(np.radians(30.0))
    print(f"100-element array at 30 deg: -3 dB width {width:.3f} deg")
    print(f"  theory 0.886*lambda/(N*d) broadside = {theory:.3f} deg")
    print(f"  broadened by 1/cos(30 deg)        = {theory_scanned:.3f} deg")
    print()
    limit = grating_lobe_limit(30.0)
    print(f"No-grating-lobe bound at 30 deg: d/lambda <= {limit:.3f}")
    print(f"  half-wavelength design sits {(1.0 - 0.5 / limit) * 100.0:.0f}% below it.")
    wide = LinearArray(elements=16, spacing_over_lambda=1.0, frequency=RF_FREQUENCY)
    wide_phases = wide.exact_phases(30.0)
    levels = wide.array_factor([-30.0, 30.0], wide_phases)
    print(
        "  At d=lambda the image lobe is a second full-gain beam: "
        f"AF(-30 deg)={levels[0]:.3f}, AF(+30 deg)={levels[1]:.3f}."
    )
    print("  That, not appearance, is why tactical arrays use d=lambda/2.")


def _squint_figure(array: LinearArray) -> None:
    ratios = np.linspace(0.85, 1.15, 121)
    figure, axes = plt.subplots(1, 2, figsize=(13.5, 4.9))
    figure.suptitle(
        "Wideband squint: a fixed phase ramp is only pointed at its design frequency",
        fontsize=12,
    )

    _section("A4. WIDEBAND SQUINT: sin(theta) = sin(theta0) / (f/f0)")
    print(f"{'theta0':>8}{'f/f0':>7}{'pointing':>11}{'gain@theta0':>14}")
    print("-" * 78)
    for theta0 in (15.0, 30.0, 45.0):
        phases = array.exact_phases(theta0)
        for ratio in (0.9, 1.0, 1.1):
            pointing = squint_angle(theta0, ratio)
            response = float(
                array.array_factor(theta0, phases, frequency_ratio=ratio)[0]
            )
            gain = 20.0 * np.log10(max(response, 1e-6))
            print(
                f"{theta0:>7.0f}deg{ratio:>7.2f}"
                f"{pointing:>10.2f}deg{gain:>+13.2f} dB"
            )
        print()

    plot_squint(
        axes[0],
        (15.0, 30.0, 45.0),
        ratios,
        title="realised pointing vs drive frequency",
    )
    grid = np.arange(-90.0, 90.0 + 0.05, 0.05)
    phases = array.exact_phases(30.0)
    for ratio, color in ((0.9, "tab:blue"), (1.0, "black"), (1.1, "tab:red")):
        axes[1].plot(
            grid,
            array.pattern_db(grid, phases, frequency_ratio=ratio),
            color=color,
            linewidth=1.2,
            label=f"f/f0={ratio:.2f}: {squint_angle(30.0, ratio):.1f} deg",
        )
    axes[1].axvline(30.0, color="0.35", linestyle=":", linewidth=1.0)
    axes[1].set_xlim(-90.0, 90.0)
    axes[1].set_ylim(-40.0, 2.0)
    axes[1].set_xlabel("theta (deg)")
    axes[1].set_ylabel("array factor (dB)")
    axes[1].set_title("same ramp, three frequencies, commanded 30 deg", fontsize=10)
    axes[1].grid(alpha=0.3)
    axes[1].legend(fontsize=7)

    _save(figure, OUTPUT_DIR / "beam_squint.png")
    print("  A 10% frequency shift at 30 deg moves the beam about 3 deg. In a")
    print("  5G wideband link that is a pointing budget, not a rounding error.")
    print("  True-time-delay steering avoids it; phase shifters do not.")


def _transport_figure() -> None:
    paths = transport_budget(100e3, pointing_jitter_rad=1e-6)
    _section("B1. TRANSPORT LOSS OVER 100 km (0.3 m aperture at both ends)")
    print(f"{'carrier / medium':<38}{'net loss':>14}")
    print("-" * 78)
    for path in paths:
        print(f"{path.label:<38}{path.loss_db:>12,.1f} dB")
    print()
    print("  A 10 GHz state in room-temperature coax dies in metres. Free-space")
    print("  RF reaches 100 km at ~93 dB but needs line of sight and both ends")
    print("  aimed. Fibre at 1550 nm costs 22 dB, and 1550 nm is also where")
    print("  low-noise detectors and erbium amplifiers live. That gap is the")
    print("  whole reason to transduce.")
    jittered_optical = free_space_link_loss(
        100e3, OPTICAL_WAVELENGTH, 0.3, 0.3, pointing_jitter_rad=10e-6
    )
    print(
        f"  At 10 urad jitter the optical free-space row becomes "
        f"{jittered_optical:.1f} dB. Sub-microradian pointing is the real"
    )
    print("  tax on free-space optical; fibre does not care.")

    temperatures = np.geomspace(0.01, 300.0, 160)
    microwave = np.array(
        [thermal_occupation(RF_FREQUENCY, value) for value in temperatures]
    )
    optical = np.array(
        [
            thermal_occupation(SPEED_OF_LIGHT / OPTICAL_WAVELENGTH, value)
            for value in temperatures
        ]
    )
    efficiency = np.geomspace(1e-3, 1.0, 200)
    added_noise = 1e-3
    calibrated = added_noise / efficiency
    uncalibrated = (1.0 - efficiency + added_noise) / efficiency

    figure, axes = plt.subplots(1, 3, figsize=(18.0, 5.6))
    figure.suptitle(
        "Quantum frequency conversion: the transport win and the noise bill",
        fontsize=13,
    )
    plot_transport_budget(
        axes[0], paths, title="why the carrier has to change"
    )
    plot_thermal_occupation(
        axes[1],
        temperatures,
        microwave,
        optical,
        title="thermal photons per mode",
    )
    plot_noise_budget(
        axes[2],
        efficiency,
        calibrated,
        uncalibrated,
        title="input-referred excess noise",
    )
    _save(figure, OUTPUT_DIR / "qfc_transport.png")

    _section("B2. WHAT CONVERSION COSTS A QUANTUM STATE (added noise 1e-3 SNU)")
    print(
        f"{'eta':>8}{'calibrated':>13}{'uncalibrated':>15}"
        f"{'conservative verdict':>26}"
    )
    print("-" * 78)
    for value in (1e-3, 1e-2, 0.1, 0.5, 0.9, 0.99):
        device = Transducer(value, added_noise)
        trusted = device.excess_noise(calibrated=True)
        untrusted = device.excess_noise(calibrated=False)
        verdict = (
            "hopeless" if untrusted > 10.0 else
            "rate crushed" if untrusted > 0.1 else
            "usable"
        )
        print(f"{value:>8.3g}{trusted:>13.2g}{untrusted:>15.2g}{verdict:>26}")
    print()
    print("  Calibrated loss is not excess noise: fibre CV-QKD runs over tens")
    print("  of dB because the loss is known. The right-hand column is the")
    print("  uncalibrated worst case, where a drifting converter looks exactly")
    print("  like an eavesdropper with a lossy component. High efficiency is")
    print("  necessary; stable, characterised efficiency is what is sufficient.")
    print("  Note that calibrated efficiency still has to clear the budget on")
    print("  its own: at eta=1e-3 the 1e-3 SNU floor becomes 1 SNU input-referred.")

    _section("B3. THE THERMAL FLOOR AT THE MICROWAVE END (10 GHz)")
    print(f"{'bath':>10}{'thermal photons':>18}{'vs 0.1 SNU':>16}")
    print("-" * 78)
    for temperature in (290.0, 4.0, 1.0, 0.02):
        occupation = thermal_occupation(RF_FREQUENCY, temperature)
        ratio = occupation / 0.1
        factor = f"{ratio:,.0f}x" if ratio >= 1.0 else f"{ratio:.1e}x"
        print(f"{temperature:>8.2f} K{occupation:>18.3g}{factor:>16}")
    print()
    print("  Conversion is linear and photon-number preserving, so the thermal")
    print("  bath converts along with the signal. At 4 K the floor is already")
    print("  78x above a QKD-grade budget; at 290 K it is 6,000x. Only a")
    print("  millikelvin stage clears it.")
    print()
    print("  Architecture consequence: the antenna can live in the sky, but")
    print("  the transducer, pump laser and detector cannot. A tactical")
    print("  microwave quantum link is a cryostat on a mast with a photonic")
    print("  tether -- fine for a fixed node, not for a man-portable one.")
    print("  For classical traffic only the 70 dB transport advantage matters,")
    print("  which is why QFC reaches the field there first.")


def _state_figure() -> None:
    squeezing_in = 10.0
    input_covariance = squeezed_vacuum_covariance(squeezing_in)
    cases = [
        ("input: 10 dB squeezed", input_covariance),
        (
            "eta=0.2, n_add=0.05",
            Transducer(0.2, 0.05).convert_covariance(input_covariance),
        ),
    ]
    q = np.linspace(-8.0, 8.0, 241)

    figure, axes = plt.subplots(2, 2, figsize=(12.5, 10.0))
    figure.suptitle(
        "A fragile state through the converter: covariance, squeezing, and the bath",
        fontsize=13,
    )
    for axis, (title, covariance) in zip(axes[0], cases, strict=True):
        plot_wigner(
            axis,
            q,
            q,
            wigner_function(q, q, covariance),
            title=f"{title}\n{covariance_squeezing_db(covariance):+.1f} dB",
        )

    efficiency = np.linspace(0.01, 1.0, 200)
    curves = {
        f"n_add={noise:g}": [
            Transducer(value, noise).surviving_squeezing_db(squeezing_in)
            for value in efficiency
        ]
        for noise in (0.0, 0.01, 0.1)
    }
    plot_squeezing_survival(
        axes[1, 0],
        efficiency,
        curves,
        title="10 dB in, no bath: efficiency alone is survivable",
    )

    temperatures = np.geomspace(0.05, 20.0, 120)
    for value, color in ((1.0, "tab:green"), (0.5, "tab:blue")):
        curve = [
            Transducer(value).surviving_squeezing_db(
                squeezing_in,
                bath_occupation=thermal_occupation(RF_FREQUENCY, temperature),
            )
            for temperature in temperatures
        ]
        axes[1, 1].semilogx(
            temperatures, curve, color=color, linewidth=1.6, label=f"eta={value:g}"
        )
    axes[1, 1].axhline(0.0, color="0.35", linestyle=":", linewidth=1.0)
    axes[1, 1].set_xlabel("bath temperature (K)", fontsize=8)
    axes[1, 1].set_ylabel("output squeezing (dB)", fontsize=8)
    axes[1, 1].set_title("10 GHz bath: 4 K erases 10 dB of squeezing", fontsize=10)
    axes[1, 1].grid(alpha=0.3)
    axes[1, 1].legend(fontsize=7)
    _save(figure, OUTPUT_DIR / "qfc_state.png")

    _section("B4. DOES THE SQUEEZED STATE SURVIVE?")
    print(
        f"{'eta':>8}{'cold':>12}{'bath 4 K':>12}{'bath 1 K':>12}"
        f"{'bath 0.02 K':>13}"
    )
    print("-" * 78)
    for efficiency_value in (0.2, 0.5, 0.9, 1.0):
        device = Transducer(efficiency_value)
        cold = device.surviving_squeezing_db(squeezing_in)
        row = [device.surviving_squeezing_db(
            squeezing_in,
            bath_occupation=thermal_occupation(RF_FREQUENCY, temperature),
        ) for temperature in (4.0, 1.0, 0.02)]
        print(
            f"{efficiency_value:>8.2f}{cold:>+11.1f} dB"
            f"{row[0]:>+10.1f} dB{row[1]:>+10.1f} dB{row[2]:>+11.1f} dB"
        )
    print()
    print("  Negative = still squeezed. A 4 K bath at 10 GHz carries 7.8")
    print("  photons per mode, so the state is thermal before the converter")
    print("  touches it. The converter is not what kills the microwave")
    print("  squeezing; the bath is. QFC should be argued for the transport")
    print("  and detector advantage, not as a way around cryogenics.")


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    print(RULE)
    print("WAVEFRONT + QFC LAB -- cyclic phase resolution and microwave-to-optics")
    print(
        f"Python {sys.version_info.major}.{sys.version_info.minor}."
        f"{sys.version_info.micro} | numpy {np.__version__}"
    )
    print(RULE)
    array = LinearArray()
    print(
        f"array: N={array.elements}, f={array.frequency / 1e9:.1f} GHz, "
        f"lambda={array.wavelength * 100:.2f} cm, d={array.spacing * 100:.2f} cm, "
        f"aperture={array.aperture * 100:.1f} cm, "
        f"Fresnel distance={array.fresnel_distance:.2f} m"
    )
    optical_frequency = SPEED_OF_LIGHT / OPTICAL_WAVELENGTH
    print(
        f"10 GHz photon energy {photon_energy(RF_FREQUENCY):.3e} J; "
        f"1550 nm photon energy {photon_energy(optical_frequency):.3e} J"
    )
    _wavefront_figure(array)
    _quantization_figure(array)
    _beam_report(array)
    _squint_figure(array)
    _transport_figure()
    _state_figure()
    print()
    print(RULE)
    print(f"done -- figures in {OUTPUT_DIR}")
    print(RULE)


if __name__ == "__main__":
    main()
