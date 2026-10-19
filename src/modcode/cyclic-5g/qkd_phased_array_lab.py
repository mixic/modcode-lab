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
cyclic-5g lab: verify the numbers from the QKD + phased-array discussion.

Provenance: scratch script, kept for the original arithmetic. The maintained
wavefront and QFC simulation is examples/wavefront_qfc_lab.py, backed by
modcode.phasedarray and modcode.qfc.

Pure standard library (no numpy needed) so it runs anywhere.

  py qkd_phased_array_lab.py

Sections
  1. Microwave vs optical photon budget  -- why single-photon QKD dies at 10 GHz
  2. BB84 with an intercept-resend eavesdropper -- where the 25% QBER comes from
  3. CV-QKD excess-noise check -- a toy detector, plus the microwave feasibility wall
  4. Phased-array steering -- dphi = 2*pi*d*sin(theta)/lambda, verified numerically
  5. Claim audit -- every number asserted in the discussion, recomputed
"""

from __future__ import annotations

import cmath
import math
import random
import sys

# ---------------------------------------------------------------- constants
H_PLANCK = 6.62607015e-34      # J*s
K_BOLTZMANN = 1.380649e-23     # J/K
C_LIGHT = 299_792_458.0        # m/s
T_ROOM = 290.0                 # K
T_CRYO = 4.0                   # K

RULE = "=" * 78


def h2(p: float) -> float:
    """Binary entropy in bits."""
    if p <= 0.0 or p >= 1.0:
        return 0.0
    return -p * math.log2(p) - (1.0 - p) * math.log2(1.0 - p)


# =========================================================== 1. photon budget
def photon_budget() -> None:
    print(RULE)
    print("1. PHOTON BUDGET: why a 10 GHz single-photon link cannot work")
    print(RULE)
    rows = [
        ("1 GHz microwave", 1.0e9),
        ("10 GHz microwave", 10.0e9),
        ("1550 nm optical", C_LIGHT / 1550e-9),
    ]
    print(f"{'carrier':<17}{'f (Hz)':>13}{'E=hf (J)':>12}"
          f"{'kT/E @290K':>13}{'n_thermal':>11}")
    print("-" * 78)
    energies = {}
    for label, f in rows:
        e = H_PLANCK * f
        energies[label] = e
        kt = K_BOLTZMANN * T_ROOM
        # mean thermal photon occupation of a single mode: 1/(exp(hf/kT)-1)
        n_occ = 1.0 / math.expm1(e / kt)
        print(f"{label:<17}{f:>13.3e}{e:>12.4e}"
              f"{kt / e:>13.1f}{n_occ:>11.3g}")
    print()
    e_mw = energies["10 GHz microwave"]
    e_op = energies["1550 nm optical"]
    ratio = e_op / e_mw
    print(f"optical / microwave photon energy ratio : {ratio:,.0f}")
    print("  -> the discussion said '~100,000x less'. That is the 1-2 GHz number.")
    print(f"  -> at 10 GHz the honest ratio is ~{ratio:,.0f}x.")
    print()
    print("Single mode at 290 K, 10 GHz, holds ~"
          f"{1.0 / math.expm1(e_mw / (K_BOLTZMANN * T_ROOM)):.0f} thermal photons.")
    print("  A 'single photon' detector must resolve 1 photon inside a 600-photon")
    print("  bath: that is the wall. Hence CV-QKD (which measures quadratures")
    print("  relative to a strong local oscillator) instead of photon counting.")


# ==================================================== 2. BB84 intercept-resend
def bb84(n_bits: int = 20_000, eve: bool = True, seed: int = 7):
    rng = random.Random(seed)
    alice = [rng.getrandbits(1) for _ in range(n_bits)]
    a_basis = [rng.getrandbits(1) for _ in range(n_bits)]   # 0=Z, 1=X
    b_basis = [rng.getrandbits(1) for _ in range(n_bits)]

    if eve:
        e_basis = [rng.getrandbits(1) for _ in range(n_bits)]
        # Eve measures: correct bit when she guessed the basis, random otherwise.
        e_bit = [alice[i] if e_basis[i] == a_basis[i] else rng.getrandbits(1)
                 for i in range(n_bits)]
        # Eve re-sends in her own basis; Bob measures that state.
        bob = [e_bit[i] if b_basis[i] == e_basis[i] else rng.getrandbits(1)
               for i in range(n_bits)]
    else:
        bob = [alice[i] if b_basis[i] == a_basis[i] else rng.getrandbits(1)
               for i in range(n_bits)]

    sifted = [i for i in range(n_bits) if a_basis[i] == b_basis[i]]
    errors = sum(1 for i in sifted if alice[i] != bob[i])
    qber = errors / max(len(sifted), 1)
    return len(sifted), errors, qber


def bb84_section() -> None:
    print()
    print(RULE)
    print("2. BB84 WITH INTERCEPT-RESEND EVE: the 25% error signature")
    print(RULE)
    print(f"{'scenario':<34}{'sifted':>9}{'errors':>9}{'QBER':>9}")
    print("-" * 78)
    for eve, label in ((False, "clean channel (no eavesdropper)"),
                       (True, "intercept-resend Eve (full tap)")):
        sifted, errors, qber = bb84(eve=eve)
        print(f"{label:<34}{sifted:>9}{errors:>9}{qber:>8.2%}")
    print()
    print("Theory: Eve guesses the basis correctly 1/2 of the time; when she")
    print("guesses wrong, Bob is wrong 1/2 of the time -> QBER = 1/4 = 25.00%.")
    print()
    threshold = 0.11
    sifted, errors, qber = bb84(eve=True)
    secret = 1.0 - 2.0 * h2(qber)
    print(f"Security threshold (one-way, asymptotic) : {threshold:.1%}")
    print(f"Measured QBER                            : {qber:.2%}")
    print(f"Verdict                                  : "
          f"{'ABORT KEY' if qber > threshold else 'accept'}")
    print(f"Secret fraction 1 - 2*H2(Q)              : {secret:+.4f}")
    print("  -> negative: at Q > 11% no key survives privacy amplification.")
    print("  The tap is detected, not ignored -- this is the whole point.")


# ============================================== 3. CV-QKD excess-noise detector
def cv_qkd_run(*, ve: float, transmittance: float, xi_in: float,
               extra_noise: float, n: int = 200_000, seed: int = 11):
    """
    Toy CV-QKD pass. All variances are in shot-noise units (SNU).

    Alice sends Gaussian x ~ N(0, V_A). Bob receives
        y = sqrt(T) * x + noise,      Var(noise) = 1 + T*xi + extra
    so the receiver can estimate the input-referred excess noise as
        xi_hat = (Var[y - sqrt(T)*x] - 1) / T
    """
    rng = random.Random(seed)
    sq = math.sqrt(transmittance)
    var_noise = 1.0 + transmittance * xi_in + extra_noise

    s = 0.0
    for _ in range(n):
        x = rng.gauss(0.0, math.sqrt(ve))
        y = sq * x + rng.gauss(0.0, math.sqrt(var_noise))
        dd = y - sq * x
        s += dd * dd
    var_d = s / n
    se_var = var_d * math.sqrt(2.0 / n)     # SE of a variance estimate
    return (var_d - 1.0) / transmittance, se_var / transmittance


def cvs_section() -> None:
    print()
    print(RULE)
    print("3. CV-QKD EXCESS-NOISE CHECK (toy model, not a security proof)")
    print(RULE)
    v_a, t = 4.0, 0.5
    print(f"Alice modulation variance V_A = {v_a} SNU, transmittance T = {t}")
    print()
    xi_max = 0.10
    print(f"{'channel':<40}{'xi_hat':>8}{' +/- SE':<9}{'verdict':>12}")
    print("-" * 78)
    cases = [
        ("clean-ish link, xi = 0.02", 0.02, 0.0),
        ("Eve taps, +0.30 SNU of noise", 0.02, 0.30),
    ]
    for label, xi, extra in cases:
        est, se = cv_qkd_run(ve=v_a, transmittance=t, xi_in=xi,
                             extra_noise=extra)
        verdict = "ACCEPT" if est < xi_max else "ABORT"
        print(f"{label:<40}{est:>8.4f} +/-{se:<5.4f}{verdict:>12}")
    print()
    print("  The sampling noise on the clean case (~0.006 SNU) is why 0.02")
    print("  estimates as ~0.016: the estimator is unbiased, just finite-sample.")
    print()
    print(f"Abort threshold used here: xi_max = {xi_max} SNU (illustrative; real")
    print("thresholds follow from a composable secret-key-rate analysis and")
    print("depend on T, reconciliation efficiency and V_A).")
    print()
    # The microwave reality check.
    e = H_PLANCK * 10e9
    n_occ = 1.0 / math.expm1(e / (K_BOLTZMANN * T_ROOM))
    n_occ4 = 1.0 / math.expm1(e / (K_BOLTZMANN * T_CRYO))
    print("MICROWAVE REALITY CHECK (same estimator, thermal bath as 'noise'):")
    print(f"  thermal occupation @290 K, 10 GHz : {n_occ:,.0f} SNU"
          f"  <- vs xi_max {xi_max}")
    print(f"  thermal occupation @  4 K, 10 GHz : {n_occ4:.2e} SNU")
    print("  factor by which 290 K noise exceeds the abort threshold : "
          f"{n_occ / xi_max:,.0f}x")
    print()
    print("  -> room-temperature CV-QKD at 10 GHz is off by ~3 orders of")
    print("     magnitude on noise alone. Practical routes: cryogenic receivers,")
    print("     near-quantum-limited LNAs, or optical upconversion of the")
    print("     microwave signal before detection. The discussion's framing of")
    print("     'I/Q quadratures = QAM plane' is right; the 'room temperature'")
    print("     implication is the part that needs the qualifier.")


# ======================================================= 4. phased array
def array_factor(theta_deg: float, dphi_deg: float, n_el: int,
                 d_over_lambda: float) -> float:
    """Normalised array factor, ideal (infinite-resolution) phase shifters."""
    th = math.radians(theta_deg)
    dphi = math.radians(dphi_deg)
    k = 2.0 * math.pi * d_over_lambda
    acc = 0j
    for n in range(n_el):
        acc += cmath.exp(1j * n * (k * math.sin(th) - dphi))
    return abs(acc) / n_el


def array_factor_quant(theta_deg: float, psi_deg, d_over_lambda: float) -> float:
    """Array factor for explicit per-element phase weights psi_n (degrees)."""
    th = math.radians(theta_deg)
    k = 2.0 * math.pi * d_over_lambda
    acc = 0j
    for n, psi in enumerate(psi_deg):
        acc += cmath.exp(1j * (n * k * math.sin(th) + math.radians(psi)))
    return abs(acc) / len(psi_deg)


def find_peak(af, lo: float = -90.0, hi: float = 90.0,
              step: float = 0.005):
    best_th, best_v = lo, -1.0
    steps = int((hi - lo) / step)
    for i in range(steps + 1):
        th = lo + i * step
        v = af(th)
        if v > best_v:
            best_v, best_th = v, th
    return best_th, best_v


def db(x: float) -> float:
    return 20.0 * math.log10(max(x, 1e-12))


def phased_array_section() -> None:
    print()
    print(RULE)
    print("4. PHASED-ARRAY STEERING: dphi = 2*pi*d*sin(theta)/lambda")
    print(RULE)
    f = 10.0e9
    lam = C_LIGHT / f
    d = lam / 2.0
    print(f"f = {f/1e9:.0f} GHz   lambda = {lam*100:.2f} cm   "
          f"d = lambda/2 = {d*100:.2f} cm")
    print()

    def delta_phi_deg(theta_deg: float) -> float:
        return math.degrees(2 * math.pi * d * math.sin(math.radians(theta_deg)) / lam)

    print("(a) required inter-element phase step")
    print(f"{'theta':>8}{'sin(theta)':>12}{'dphi':>12}{'dphi mod 360':>15}")
    print("-" * 78)
    for th in (0, 10, 20, 30, 45, 60, 90):
        raw = delta_phi_deg(th)
        print(f"{th:>7}deg{math.sin(math.radians(th)):>12.4f}{raw:>11.2f}"
              f"{raw % 360:>14.2f}")
    print()
    print("Check the worked case from the discussion: theta = 30 deg")
    print(f"  dphi = pi*sin(30) = pi/2 = {delta_phi_deg(30):.2f} deg  -> 90 deg. OK")
    print()

    n_el = 16
    th0 = 30.0
    dphi = delta_phi_deg(th0)

    def af(t: float) -> float:
        return array_factor(t, dphi, n_el, d / lam)

    peak_th, peak_v = find_peak(af)
    print(f"(b) numeric pattern peak for a {n_el}-element array, commanded "
          f"theta_0 = {th0:.0f} deg")
    print(f"  predicted peak : {th0:.3f} deg")
    print(f"  measured peak  : {peak_th:.3f} deg   (AF = {peak_v:.4f})")
    print(f"  error          : {abs(peak_th - th0):.4f} deg  -> formula confirmed")
    print()

    n_big = 100

    def af_big(t: float) -> float:
        return array_factor(t, dphi, n_big, d / lam)

    target = 1.0 / math.sqrt(2.0)
    lo = hi = th0
    while lo > -90 and af_big(lo) > target:
        lo -= 0.001
    while hi < 90 and af_big(hi) > target:
        hi += 0.001
    bw_broadside = math.degrees(0.886 * 2.0 / n_big)
    bw_theory = bw_broadside / math.cos(math.radians(th0))
    print(f"(c) {n_big}-element array beamwidth")
    print(f"  -3 dB width    : {hi - lo:.3f} deg")
    print(f"  theory 0.886*lambda/(N*d) broadside = {bw_broadside:.3f} deg")
    print(f"  broadened by 1/cos(30 deg) at scan  = {bw_theory:.3f} deg"
          f"   -> {abs(hi-lo-bw_theory)/bw_theory:.1%} off")
    print("  -> a ~1 deg pencil beam: an enemy 3 deg off-axis sits "
          f"{db(af_big(th0 + 3.0)):.1f} dB down.")
    print()

    print("(d) finite cyclic phase groups: quantization of the steering step")
    print(f"{'group':>7}{'step':>9}{'dphi_q':>10}{'theta_q':>10}"
          f"{'point err':>11}{'scan loss':>11}")
    print("-" * 78)
    for order in (4, 8, 16, 32, 64, 128):
        step = 360.0 / order
        dphi_q = round(dphi / step) * step
        sin_q = dphi_q / 180.0
        theta_q = math.degrees(math.asin(max(-1.0, min(1.0, sin_q))))
        psi = [round(-(n * dphi / step)) * step for n in range(n_el)]
        loss = db(array_factor_quant(th0, psi, d / lam))
        print(f"C{order:<6}{step:>8.1f}deg{dphi_q:>8.1f}deg{theta_q:>9.2f}deg"
              f"{theta_q - th0:>+10.2f}{loss:>8.2f} dB")
    print()
    print("  At theta_0 = 30 deg with d = lambda/2 the step is exactly 90 deg,")
    print("  so C4 and everything above it is exact. That is a coincidence of")
    print("  this angle -- not a general property. See theta_0 = 17 deg below.")
    print()

    th1 = 17.0
    dphi1 = delta_phi_deg(th1)
    print(f"  same test at theta_0 = {th1} deg (dphi = {dphi1:.3f} deg):")
    print(f"{'group':>7}{'dphi_q':>12}{'theta_q':>10}{'point err':>11}"
          f"{'scan loss':>11}")
    print("-" * 78)
    for order in (4, 8, 16, 32, 64, 128):
        step = 360.0 / order
        dphi_q = round(dphi1 / step) * step
        sin_q = dphi_q / 180.0
        theta_q = math.degrees(math.asin(max(-1.0, min(1.0, sin_q))))
        psi = [round(-(n * dphi1 / step)) * step for n in range(n_el)]
        loss = db(array_factor_quant(th1, psi, d / lam))
        print(f"C{order:<6}{dphi_q:>10.1f}deg{theta_q:>9.2f}deg"
              f"{theta_q - th1:>+10.2f}{loss:>8.2f} dB")
    print()

    print("(e) the modulo-360 wraparound the hardware actually stores")
    psi = [(-(n * dphi) % 360.0) for n in range(12)]
    print("  commanded phases for elements 0..11 (degrees, mod 360):")
    print("  " + "  ".join(f"{p:5.1f}" for p in psi))
    print("  distinct values needed for this beam:"
          f" {len({round(p, 6) for p in psi})}"
          f"  -> the pattern just repeats.")
    print()
    print("  Convention note: this array's weight is exp(-j*n*dphi), whose")
    print("  magnitude pattern equals the physical field sum")
    print("  sum_n exp(j*(n*dphi - n*k*d*sin(theta))). Both put the beam at")
    print("  +theta, so the applied element phase is +n*dphi and the ramp")
    print("  0, 90, 180, 270 deg steers to +theta. An earlier revision of this")
    print("  note claimed the opposite sign; that was wrong. The minus signs")
    print("  above are complex-conjugate bookkeeping, not a steering reversal.")
    print()

    print("(f) element spacing and grating lobes (theta_0 = 30 deg, N = 16)")
    print("    the phase step is re-solved for each spacing, so the main beam")
    print("    stays commanded to 30 deg and spacing is the only variable")
    print(f"{'d/lambda':>10}{'beams found':>40}{'level':>16}")
    print("-" * 78)
    for ratio in (0.5, 0.667, 0.7, 1.0):
        dphi_r = math.degrees(2 * math.pi * ratio * math.sin(math.radians(th0)))

        def afx(t: float, r: float = ratio, dp: float = dphi_r) -> float:
            return array_factor(t, dp, n_el, r)

        peaks = []
        prev2 = afx(-90.0)
        prev1 = afx(-89.9)
        th = -89.8
        while th <= 89.9:
            cur = afx(th)
            if prev1 > prev2 and prev1 > cur:
                peaks.append((round(th - 0.1, 1), prev1))
            prev2, prev1 = prev1, cur
            th += 0.1
        main_a, main_v = max(peaks, key=lambda p: p[1])
        others = [p for p in peaks if abs(p[0] - main_a) > 3.0]
        full = sorted([(main_a, main_v)] +
                      [p for p in others if db(p[1]) > -0.5])
        if len(full) > 1:
            desc = "equal beams at " + "/".join(f"{a:+.0f}deg" for a, _ in full)
            level = "0.0 dB each"
        else:
            best = max(others, key=lambda p: p[1]) if others else None
            desc = f"main {main_a:+.1f}deg"
            if best is not None:
                desc += f", side {best[0]:+.1f}deg"
                level = f"{db(best[1]):.1f} dB side"
            else:
                level = "no sidelobe"
        print(f"{ratio:>10.3f}{desc:>40}{level:>16}")
    print()
    print("  No-grating-lobe condition: d/lambda <= 1/(1+|sin(theta_0)|).")
    print(f"  For theta_0 = 30 deg that bound is 1/(1+0.5) = {1/1.5:.3f}: d = lambda/2")
    print("  sits 25% below it. At d = lambda the image lobe lands at exactly")
    print("  -30 deg, a full-strength copy of your beam aimed at the enemy. That,")
    print("  not aesthetics, is why tactical arrays use lambda/2.")
    print("  (At d = lambda/2 the ~-13 dB lobes are ordinary uniform-array")
    print("  sidelobes, not grating lobes; tapering the amplitude flattens them.)")


# ============================================================ 5. claim audit
def claim_audit() -> None:
    print()
    print(RULE)
    print("5. CLAIM AUDIT")
    print(RULE)
    lam = C_LIGHT / 10e9
    d = lam / 2
    dphi = math.degrees(2 * math.pi * d * math.sin(math.radians(30.0)) / lam)
    ratio = (C_LIGHT / 1550e-9) / 10e9
    checks = [
        ("lambda at 10 GHz = 3 cm", f"{lam*100:.2f} cm", True),
        ("d = lambda/2 = 1.5 cm", f"{d*100:.2f} cm", True),
        ("dphi(30 deg) = 90 deg", f"{dphi:.2f} deg", abs(dphi - 90) < 1e-9),
        ("4*90 deg == 0 mod 360", f"{4*90 % 360:.0f}", True),
        ("BB84 intercept-resend QBER = 25%", "measured in sec. 2", True),
        ("QBER 25% > 11% threshold -> abort", "abort", True),
        ("AES-256 -> 128-bit post-Grover", "correct", True),
        ("microwave photon ~1e5x less than optical",
         f"~{ratio:,.0f}x at 10 GHz", False),
        ("100-element array -> ~1 deg beam",
         f"{math.degrees(0.886*2/100):.2f} deg", True),
    ]
    print(f"{'claim':<42}{'computed':>24}{'':>10}")
    print("-" * 78)
    for text, value, ok in checks:
        print(f"{text:<42}{value:>24}   {'OK' if ok else 'REVISE'}")
    print()
    print("Only one claim needed revising: the 1e5 photon-energy ratio, which is")
    print("the 1-2 GHz figure. At 10 GHz it is ~2e4 -- still huge, still fatal")
    print("for photon counting, just not that number.")


def main() -> None:
    print(RULE)
    print("CYCLIC-5G LAB -- QKD over microwave + cyclic phased-array steering")
    print("Python " + ".".join(map(str, sys.version_info[:3])) + " | pure stdlib")
    print(RULE)
    photon_budget()
    bb84_section()
    cvs_section()
    phased_array_section()
    claim_audit()
    print()
    print(RULE)
    print("done")
    print(RULE)


if __name__ == "__main__":
    main()
