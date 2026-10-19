#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""
cyclic-5g lab, part 2
  A. Wavefront synthesis: what cyclic phase steps actually do to the field
  B. Quantum frequency conversion: why lifting microwave states to optics
     is the only way to route them, and what it costs

Provenance: this is the original scratch script. The maintained version is
examples/wavefront_qfc_lab.py, backed by modcode.phasedarray and modcode.qfc.
The maintained version separates per-element phase rounding from uniform-step
rounding; this script's quantization paragraph below has been corrected to
match the per-element case it actually plots.

Run:  py wavefront_qfc_lab.py
Writes PNGs to ./figures/ and a numeric report to stdout.

Geometry convention
-------------------
Elements sit at x_n = (n - (N-1)/2)*d on z = 0, phase centre looking up +z.
The APPLIED element phase is psi_n = n*dphi with
    dphi = 2*pi*d*sin(theta0)/lambda,
and that steers the beam to +theta0. The 0, 90, 180, 270 degree ramp is
therefore a beam to +theta, matching the classic tutorial picture.
"""

from __future__ import annotations

import math
import os
import sys

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ------------------------------------------------------------------ constants
C = 299_792_458.0
H = 6.62607015e-34
KB = 1.380649e-23
RULE = "=" * 78

HERE = os.path.dirname(os.path.abspath(__file__))
FIGDIR = os.path.join(HERE, "figures")

F_RF = 10.0e9
LAM = C / F_RF
D = LAM / 2.0
N_EL = 16
XS = (np.arange(N_EL) - (N_EL - 1) / 2.0) * D


# ===================================================== A. field synthesis
def synth_field(x, z, psi, lam=LAM, xs=XS):
    """
        E(x,z) = sum_n  exp(j*(k*r_n + psi_n)) / sqrt(1 + r_n/lambda)

    Cylindrical spreading from each element. This is the physical field, not
    an array-factor approximation, so the plots show genuine near-field
    curvature.
    """
    X, Z = np.meshgrid(x, z)
    K = 2.0 * math.pi / lam
    Ec = np.zeros_like(X, dtype=complex)
    for xn, p in zip(xs, psi, strict=False):
        r = np.hypot(X - xn, Z)
        Ec += np.exp(1j * (K * r + p)) / np.sqrt(1.0 + r / lam)
    return X, Z, Ec


def far_field_peak(psi, lam=LAM, xs=XS, radius=4.0, n_ang=7201):
    """
    Independent check: sum the actual cylindrical waves on a large circle and
    find where the intensity peaks. Uses no array-factor formula at all, so it
    catches synthesis errors (wrong sign, wrong spacing, wrong phase ramp).

    radius = 4 m is well beyond the Fresnel distance L^2/lambda = 1.7 m.
    """
    th = np.linspace(-89.0, 89.0, n_ang)
    Xp = radius * np.sin(np.radians(th))
    Zp = radius * np.cos(np.radians(th))
    K = 2.0 * math.pi / lam
    E = np.zeros(n_ang, dtype=complex)
    for xn, p in zip(xs, psi, strict=False):
        r = np.hypot(Xp - xn, Zp)
        E += np.exp(1j * (K * r + p)) / np.sqrt(1.0 + r / lam)
    return float(th[int(np.argmax(np.abs(E)))])


def ramp_exact(theta0):
    dphi = math.degrees(2.0 * math.pi * (D / LAM) * math.sin(math.radians(theta0)))
    return [n * dphi for n in range(N_EL)], dphi


def ramp_quantized(theta0, order):
    _, dphi = ramp_exact(theta0)
    step = 360.0 / order
    return [round(n * dphi / step) * step for n in range(N_EL)]


def array_factor(theta_deg, psi_deg, d_over_lambda=D / LAM):
    th = np.radians(np.asarray(theta_deg, dtype=float))
    psi = np.radians(np.asarray(psi_deg, dtype=float))
    k = 2.0 * math.pi * d_over_lambda
    n = np.arange(len(psi))
    phase = psi[:, None] - k * n[:, None] * np.sin(th)[None, :]
    return np.abs(np.exp(1j * phase).sum(axis=0)) / len(psi)


def draw_panel(ax, th0, psi, title, show_guide=True):
    x = np.linspace(-0.36, 0.36, 640)
    z = np.linspace(0.004, 0.72, 715)
    X, Z, Ec = synth_field(x, z, psi)
    re = Ec.real
    lim = np.percentile(np.abs(re), 99.0)
    ax.imshow(re, extent=[x[0], x[-1], z[0], z[-1]], origin="lower",
              cmap="RdBu_r", vmin=-lim, vmax=lim, aspect="auto",
              interpolation="bilinear")
    if show_guide:
        slope = math.tan(math.radians(th0))
        x_end = np.clip(slope * (z[-1] - 0.05), x[0], x[-1])
        ax.plot([0.0, x_end], [0.05, z[-1]], ls="--", lw=1.0,
                color="k", alpha=0.55)
    ax.set_title(title, fontsize=10)
    ax.set_xlabel("x (m)", fontsize=8)
    ax.tick_params(labelsize=7)
    return Ec


def draw_phase_panel(ax, th0, psi):
    x = np.linspace(-0.36, 0.36, 640)
    z = np.linspace(0.004, 0.72, 715)
    _, _, Ec = synth_field(x, z, psi)
    ax.imshow(np.angle(Ec), extent=[x[0], x[-1], z[0], z[-1]], origin="lower",
              cmap="twilight", aspect="auto", interpolation="bilinear")
    slope = math.tan(math.radians(th0))
    x_end = np.clip(slope * (z[-1] - 0.05), x[0], x[-1])
    ax.plot([0.0, x_end], [0.05, z[-1]], ls="--", lw=1.0,
            color="k", alpha=0.55)
    ax.set_title(f"phase, guide at {th0:.0f} deg", fontsize=10)
    ax.set_xlabel("x (m)", fontsize=8)
    ax.tick_params(labelsize=7)


def wavefront_figure():
    cases = [0.0, 15.0, 30.0, 45.0]
    fig, axes = plt.subplots(2, 4, figsize=(17.5, 8.2))
    fig.suptitle("Wavefront shaping by cyclic phase ramps   "
                 f"(f={F_RF/1e9:.0f} GHz, lambda={LAM*100:.1f} cm, "
                 f"N={N_EL}, d=lambda/2, aperture={N_EL*D*100:.1f} cm)",
                 fontsize=13)

    print(RULE)
    print("A1. FIELD SYNTHESIS CHECK: far-field peak vs commanded angle")
    print(RULE)
    print(f"{'commanded':>11}{'dphi':>10}{'peak (independent)':>20}{'error':>10}")
    print("-" * 78)
    for col, th0 in enumerate(cases):
        psi, dphi = ramp_exact(th0)
        peak = far_field_peak(psi)
        print(f"{th0:>10.0f}deg{dphi:>9.1f}deg{peak:>19.2f}deg"
              f"{peak - th0:>+9.2f}")
        draw_panel(axes[0, col], th0, psi,
                   f"theta0={th0:.0f} deg, dphi={dphi:.1f} deg")
        draw_phase_panel(axes[1, col], th0, psi)
        if col == 0:
            axes[0, col].set_ylabel("range z (m)", fontsize=8)
            axes[1, col].set_ylabel("range z (m)", fontsize=8)

    fig.text(0.5, 0.012,
             "Dashed line = commanded direction. Fringes are wavefronts; note "
             "they are genuinely curved out to the Fresnel distance "
             f"L^2/lambda = {(N_EL*D)**2/LAM:.2f} m.",
             ha="center", fontsize=9)
    fig.tight_layout(rect=(0, 0.03, 1, 0.955))
    path = os.path.join(FIGDIR, "wavefront_steering.png")
    fig.savefig(path, dpi=140)
    plt.close(fig)
    print()
    print("  The check sums real cylindrical waves on a 4 m circle -- no")
    print("  array-factor formula involved -- and lands on the commanded")
    print("  angle to within the 0.025 deg sampling step.")
    print(f"  figure -> {path}")


def quantization_figure():
    th0 = 17.0
    psi_exact, dphi = ramp_exact(th0)
    variants = [("exact ramp", psi_exact, dphi),
                ("C8 quantized (45 deg)", ramp_quantized(th0, 8), None),
                ("C4 quantized (90 deg)", ramp_quantized(th0, 4), None)]

    fig, axes = plt.subplots(1, 4, figsize=(18.0, 4.9))
    fig.suptitle(f"Commanded to {th0:.0f} deg: what quantization does to the "
                 f"actual wavefront (dphi must be {dphi:.1f} deg)", fontsize=12)

    print()
    print(RULE)
    print("A2. QUANTIZATION: commanded 17 deg, realized pointing")
    print(RULE)
    print(f"{'shifter':>22}{'applied step':>15}{'realized peak':>16}"
          f"{'error':>10}{'loss':>9}")
    print("-" * 78)
    for col, (label, psi, _) in enumerate(variants):
        angles = np.linspace(-89.0, 89.0, 17801)
        peak = float(angles[int(np.argmax(array_factor(angles, psi)))])
        af_at_cmd = float(array_factor([th0], psi)[0])
        loss = 20 * math.log10(max(af_at_cmd, 1e-6))
        print(f"{label:>22}{psi[1]:>13.1f}deg{peak:>15.2f}deg"
              f"{peak - th0:>+9.2f}{loss:>7.1f}dB")
        draw_panel(axes[col], peak, psi,
                   f"{label}\npointing {peak:.1f} deg", show_guide=False)

    ax = axes[3]
    th = np.linspace(-90, 90, 3601)
    for (label, psi, _), color in zip(variants,
                                      ("k", "tab:orange", "tab:red"),
                                      strict=False):
        ax.plot(th, 20 * np.log10(np.maximum(array_factor(th, psi), 1e-4)),
                color=color, lw=1.4, label=label)
    ax.axvline(th0, color="gray", ls=":", lw=1)
    ax.text(th0 + 1.5, -3, "commanded", fontsize=8, color="gray", rotation=90)
    ax.set_xlim(-90, 90)
    ax.set_ylim(-40, 1)
    ax.set_xlabel("theta (deg)")
    ax.set_ylabel("array factor (dB)")
    ax.set_title("patterns", fontsize=10)
    ax.legend(fontsize=7.5)
    ax.grid(alpha=0.3)

    fig.tight_layout(rect=(0, 0, 1, 0.9))
    path = os.path.join(FIGDIR, "quantization_error.png")
    fig.savefig(path, dpi=140)
    plt.close(fig)
    print()
    print("  Per-element C4 rounding keeps the beam at ~17 deg but loses 0.8 dB")
    print("  to the staircase wavefront. Quantising the step once instead --")
    print("  every element advanced by the same 90 deg -- keeps a flat")
    print("  wavefront and full gain, but aims at 30 deg. The two errors are")
    print("  different, and only the second is a 13 deg pointing miss.")
    print(f"  figure -> {path}")


# ============================================== B. quantum frequency conversion
def thermal_occupation(f, temp):
    return 1.0 / math.expm1(H * f / (KB * temp))


def transport_rows():
    d_ant = 0.30
    fspl_rf = 20 * math.log10(4 * math.pi * 1e5 / LAM)
    gain_rf = 20 * math.log10(math.pi * d_ant / LAM)
    lam_o = 1550e-9
    fspl_o = 20 * math.log10(4 * math.pi * 1e5 / lam_o)
    gain_o = 20 * math.log10(math.pi * d_ant / lam_o)
    return [
        ("RF coax, 10 GHz (0.5 dB/m)", 0.5 * 1e5),
        ("RF waveguide, 10 GHz (0.05 dB/m)", 0.05 * 1e5),
        ("free-space RF 10 GHz, 0.3 m dishes", fspl_rf - 2 * gain_rf),
        ("free-space optical 1550 nm, 0.3 m", fspl_o - 2 * gain_o),
        ("optical fiber 1550 nm + connectors", 0.2 * 100 + 2.0),
    ]


def qfc_figure():
    fig, axes = plt.subplots(1, 2, figsize=(15.5, 6.2))
    fig.suptitle("Quantum frequency conversion: what it buys and what it costs",
                 fontsize=13)

    ax = axes[0]
    rows = transport_rows()
    labels = [r[0] for r in rows][::-1]
    vals = [max(r[1], 1.0) for r in rows][::-1]
    colors = ["tab:green", "tab:blue", "tab:purple", "tab:orange", "tab:red"]
    bars = ax.barh(labels, vals, color=colors, alpha=0.85)
    ax.set_xscale("log")
    ax.set_xlim(1, 2e5)
    ax.set_xlabel("net transport loss over 100 km (dB, log scale)")
    ax.set_title("why microwave states must change carrier", fontsize=10)
    for bar, v in zip(bars, vals, strict=False):
        ax.text(v * 1.15, bar.get_y() + bar.get_height() / 2, f"{v:,.0f} dB",
                va="center", fontsize=8)
    ax.tick_params(labelsize=8)
    ax.grid(axis="x", alpha=0.3)

    ax = axes[1]
    eta = np.logspace(-5, 0, 400)
    ax.loglog(eta, (1.0 - eta) / eta, "k-", lw=2.0,
              label="uncalibrated loss: (1-eta)/eta")

    ax.axvspan(1e-5, 1e-2, color="tab:orange", alpha=0.13)
    ax.text(2.2e-5, 2e-6, "electro-optic\n(typical)", fontsize=8,
            color="darkorange", rotation=90, va="bottom")
    ax.axvspan(1e-2, 0.5, color="tab:green", alpha=0.13)
    ax.text(3e-2, 2e-6, "piezo-optomechanical\n(reported)", fontsize=8,
            color="darkgreen", rotation=90, va="bottom")

    ax.axhline(0.1, color="tab:red", ls="--", lw=1.4)
    ax.text(1.2e-5, 0.115, "QKD-grade added-noise budget ~0.1 SNU",
            fontsize=8, color="tab:red")

    for temp, ls in ((290.0, "--"), (4.0, "-."), (0.02, ":")):
        n_th = thermal_occupation(F_RF, temp)
        ax.axhline(n_th, color="tab:blue", ls=ls, lw=1.1)
        ax.text(1.2e-5, n_th * 1.4, f"thermal bath {temp:g} K: "
                                     f"{n_th:.3g} photons",
                fontsize=7.5, color="tab:blue")

    ax.set_xlabel("transducer efficiency eta")
    ax.set_ylabel("input-referred noise (SNU)")
    ax.set_title("three walls: efficiency, added noise, thermal bath",
                 fontsize=10)
    ax.set_ylim(1e-11, 1e5)
    ax.legend(fontsize=8, loc="upper right")
    ax.grid(alpha=0.3, which="both")

    fig.tight_layout(rect=(0, 0, 1, 0.94))
    path = os.path.join(FIGDIR, "qfc_transduction.png")
    fig.savefig(path, dpi=140)
    plt.close(fig)
    return path


def qfc_report():
    print()
    print(RULE)
    print("B1. TRANSPORT LOSS OVER 100 km (0.3 m aperture at both ends)")
    print(RULE)
    print(f"{'carrier / medium':<44}{'net loss':>14}")
    print("-" * 78)
    for label, loss in transport_rows():
        print(f"{label:<44}{loss:>12,.1f} dB")
    print()
    print("  A 10 GHz state in coax dies in a few hundred metres. Free-space")
    print("  RF reaches 100 km at ~93 dB, but needs line of sight and both")
    print("  ends aimed. Fiber at 1550 nm costs 22 dB for the same distance.")
    print("  That ~70 dB gap is the entire motivation for transduction.")

    print()
    print(RULE)
    print("B2. WHAT CONVERSION COSTS A QUANTUM STATE")
    print(RULE)
    print(f"{'eta':>8}{'(1-eta)/eta':>16}{'verdict if loss is uncalibrated':>34}")
    print("-" * 78)
    for eta in (1e-4, 1e-3, 1e-2, 0.1, 0.5, 0.9, 0.99, 0.999):
        n_add = (1.0 - eta) / eta
        verdict = ("hopeless" if n_add > 10 else
                   "rate crushed" if n_add > 0.1 else "usable")
        print(f"{eta:>8.3g}{n_add:>16.3g}{verdict:>34}")

    print()
    print("  Caveat that decides the whole architecture: calibrated loss is")
    print("  NOT excess noise. Fiber CV-QKD runs over tens of dB because the")
    print("  loss is known and lives in T, not in the noise budget. The curve")
    print("  above is the worst case -- an uncalibrated transducer, or one")
    print("  whose efficiency drifts, looks exactly like an eavesdropper who")
    print("  swapped in a lossy component. High eta is necessary; stable,")
    print("  characterised eta is what makes it trustworthy.")

    print()
    print(RULE)
    print("B3. THE THERMAL FLOOR AT THE MICROWAVE END (10 GHz)")
    print(RULE)
    print(f"{'bath':>10}{'thermal photons':>20}{'vs 0.1 SNU budget':>24}")
    print("-" * 78)
    for temp in (290.0, 4.0, 1.0, 0.02):
        n = thermal_occupation(F_RF, temp)
        ratio = n / 0.1
        txt = f"{ratio:,.0f}x" if ratio >= 1 else f"{ratio:.1e}x (negligible)"
        print(f"{temp:>8.2f} K{n:>20.3g}{txt:>24}")

    print()
    print("  Conversion is linear and photon-number preserving, so the thermal")
    print("  bath converts along with the signal. At 4 K the floor is already")
    print("  78x above a QKD-grade budget; at 290 K it is 6,000x. Only a")
    print("  millikelvin stage clears it.")
    print()
    print("  Architecture consequence: the antenna can live in the sky, but")
    print("  the transducer, pump laser and detector cannot. A 'tactical")
    print("  microwave QKD link' is really a cryostat on a mast with a")
    print("  photonic tether -- which is fine for a fixed node, and is the")
    print("  reason quantum microwave routing is a backbone technology, not a")
    print("  man-portable one.")
    print()
    print("  For classical (non-quantum) traffic the calculus is gentler: the")
    print("  same transducer passes an ordinary QAM/APSK waveform, and then")
    print("  only the 70 dB transport advantage matters, not the 0.1 SNU")
    print("  budget. That is where QFC shows up in fielded systems first.")


def main():
    os.makedirs(FIGDIR, exist_ok=True)
    print(RULE)
    print("CYCLIC-5G LAB, PART 2 -- wavefront shaping and quantum transduction")
    print("Python " + ".".join(map(str, sys.version_info[:3])))
    print(RULE)
    wavefront_figure()
    quantization_figure()
    path = qfc_figure()
    qfc_report()
    print(f"  figure -> {path}")
    print()
    print(RULE)
    print("done")
    print(RULE)


if __name__ == "__main__":
    main()
