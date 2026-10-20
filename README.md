<div align="center">

# modcode-lab

**A reproducible simulation lab for digital signal modulation and forward error correction.**

*Modulation × Coding — measured, plotted, reproducible — and free as in freedom.*

[![CI](https://github.com/mixic/modcode-lab/actions/workflows/ci.yml/badge.svg)](https://github.com/mixic/modcode-lab/actions/workflows/ci.yml)
[![Python](https://img.shields.io/badge/python-3.10%20%7C%203.11%20%7C%203.12-blue.svg)](https://www.python.org/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Code style: ruff](https://img.shields.io/badge/code%20style-ruff-000000.svg)](https://github.com/astral-sh/ruff)

</div>

---

## Table of Contents

- [Why modcode-lab?](#why-modcode-lab)
- [Features](#features)
- [Installation](#installation)
- [Quickstart](#quickstart)
- [How it works](#how-it-works)
- [The math behind it](#the-math-behind-it)
- [Figures](#figures)
- [Module reference](#module-reference)
- [Documentation](#documentation)
- [Project structure](#project-structure)
- [Testing and reproducibility](#testing-and-reproducibility)
- [Extending modcode-lab](#extending-modcode-lab)
- [Roadmap](#roadmap)
- [Contributing](#contributing)
- [References](#references)
- [License](#license)

---

## Why modcode-lab?

Digital communication is the product of two inseparable disciplines:

1. **Modulation** — how bits become analog waveforms (BPSK, QPSK, QAM, APSK).
2. **Forward error correction (FEC)** — how bits are protected *before* they hit the noisy channel (Hamming, BCH, Reed-Solomon, convolutional, TCM).

Most textbooks treat them separately. Real systems (Wi-Fi, 5G, DVB-S2X) fuse them. `modcode-lab` is a **small, hackable, numerically verified sandbox** that lets you:

- swap any modulation scheme for any other with a one-line change,
- toggle FEC on/off and measure the resulting *coding gain* in dB,
- visualize constellations *before* and *after* the channel,
- reproduce every published BER curve you read in a paper,
- extend the framework with your own modulation or code in <50 lines of Python.

It's the missing bridge between a lecture on cyclic groups and a real link-budget spreadsheet.

---

## Features

### Baseline implementation (M0 + M1)

| Category | What's available |
|---|---|
| **Modulation** | BPSK, QPSK, 16-QAM, 64-QAM (Gray-mapped, unit-power normalized) |
| **FEC** | Hamming(7,4) as a *cyclic* code with generator `g(x) = x³ + x + 1` |
| **Channel** | AWGN with information-bit Eb/N0 normalization |
| **Metrics** | Bit error rate (BER) as a function of Eb/N0 |
| **Visualization** | Constellation diagrams, BER sweeps, wavefront/beam maps, and QFC trade-off plots |
| **Tests** | pytest suite covering modulation, FEC, channels, optics, arrays, QFC, and waveforms |
| **CI** | GitHub Actions on Python 3.10 / 3.11 / 3.12 |

### Experimental research prototypes

| Area | Prototype |
|---|---|
| **Cyclic geometry** | APSK rings with explicit cyclic phase actions |
| **Trellis coding** | Tail-biting rate-2/3 8-PSK example with pair-rotation equivariance |
| **Probabilistic shaping** | Maxwell-Boltzmann probabilities and periodic autocorrelation |
| **Structured FEC** | QC-LDPC base-matrix lifting and Polar coding over the BEC |
| **Phased arrays** | Linear-array wavefront synthesis with cyclic phase quantisation |
| **Quantum transduction** | Photon/thermal budgets, QFC noise, transport-loss comparison |
| **Waveforms** | PSK constellations, Zadoff–Chu generation, periodic autocorrelation |
| **CRC** | Systematic CRC-8/16/32 encoding, detection, and division traces |
| **Post-quantum structure** | Toy QC-MDPC McEliece-style circulant blocks and key-size scaling |
| **Number theory** | Full-reptend primes, cyclic multiples, multiplicative orders |
| **Residual indices** | Densities of primes with prescribed `(p-1)/ord_p(10)` index |
| **CT-IM** | Constant-weight cyclic index modulation with the Hamming(7,4) orbit |
| **Cryptography** | Toy Diffie–Hellman in a prime field |
| **Optical IM/DD** | OOK, nonnegative PAM, impairment experiments |
| **Cyclic-5g lab** | Standalone QKD/phased-array scratch scripts and generated figures |

### Planned

- **M2** — APSK (multi-ring), Rayleigh fading, multipath echo channel
- **M3** — BCH, Reed-Solomon, convolutional codes
- **M4** — Trellis-coded modulation (TCM / C-TCM)
- **M5** — Full documentation site + Jupyter notebooks

---

## Installation

Requires **Python ≥ 3.10**.

```bash
git clone https://github.com/mixic/modcode-lab.git
cd modcode-lab

python -m venv .venv
source .venv/bin/activate           # Windows: .venv\Scripts\activate

pip install -e ".[dev]"
```

## Quickstart

```python
from modcode import AWGNChannel, BPSK, Link, run_link

link = Link(modulation=BPSK(), channel=AWGNChannel())
result = run_link(link, n_bits=100_000, ebn0_db=4.0)
print(f"BER: {result.ber:.4g}")
```

Run the examples from the repository root:

```bash
python examples/hamming_demo.py
python examples/run_ber_sweep.py
python examples/explore_cyclic_structures.py
python examples/explore_qc_polar.py
python examples/qc_ldpc_shift_demo.py
python examples/wavefront_qfc_lab.py
python examples/cyclic_codes_pqc_lab.py
python examples/residual_index_lab.py
python examples/ct_im_lab.py
```

| Script | What it demonstrates |
|---|---|
| `hamming_demo.py` | Encode `1010`, inject one bit error, and correct it with Hamming(7,4) |
| `run_ber_sweep.py` | BER versus Eb/N0 for BPSK, QPSK, and QAM |
| `explore_cyclic_structures.py` | 16-APSK ring actions, shaped symbols, and the tail-biting 8-PSK TCM |
| `explore_qc_polar.py` | QC-LDPC lifting and exact BEC Polar construction |
| `qc_ldpc_shift_demo.py` | A small shift-value base matrix and its lifted parity-check matrix |
| `wavefront_qfc_lab.py` | Wavefront steering, cyclic phase quantisation, and QFC transport/noise figures |
| `cyclic_codes_pqc_lab.py` | CRC detection and the circulant-block structure of a toy QC-MDPC code |
| `residual_index_lab.py` | Predicted densities and finite prime counts for `(p-1)/ord_p(10)` |
| `ct_im_lab.py` | CT-IM constant-weight codebook, distance structure, and rate comparison |

Run the test suite with:

```bash
python -m pytest
```

The Hamming example encodes `1010` as `1010011`, introduces a one-bit error,
and demonstrates correction back to the original data.

## How it works

`Link` composes a modulation, channel, and optional FEC encoder. `run_link`
generates random data, pads it to complete codeword and symbol boundaries,
encodes it, modulates it, applies the channel, and then performs hard-decision
demodulation and FEC decoding. Its `LinkResult` provides both transmitted and
received bits and symbols, plus a BER property.

`ber_sweep` runs independent Monte Carlo trials at each Eb/N0 point. Pass it a
mapping of display names to modulation instances, or a sequence of configured
`Link` objects:

```python
import numpy as np

from modcode.metrics.ber import ber_sweep
from modcode.modulation import BPSK, QPSK, QAM

ebn0_db = np.arange(0.0, 11.0, 2.0)
curves = ber_sweep(
    {"BPSK": BPSK(), "QPSK": QPSK(), "16-QAM": QAM(16)},
    ebn0_db,
    n_bits=100_000,
    seed=2026,
)
```

The AWGN channel treats `Eb/N0` as energy per information bit. With unit
average symbol energy, the noise variance per real component is
`1 / (2 · bits_per_symbol · code_rate · Eb/N0)`. For an uncoded link, the code
rate is one.

## The math behind it

For square M-QAM, each axis has `sqrt(M)` equally spaced levels, and each
axis label is Gray-coded. The constellation is scaled so that its average
symbol energy is one. QPSK is the corresponding two-bit, unit-energy mapping;
BPSK maps bits 0 and 1 to +1 and -1, respectively.

For Hamming(7,4), the generator polynomial is `g(x) = x³ + x + 1`. Four data
bits followed by three zero bits are divided by the generator polynomial;
the three-bit remainder is appended to the original data. The decoder uses
the syndrome to correct one flipped bit per codeword.

### Cyclic groups, full-reptend primes, and Diffie–Hellman

`modcode.number_theory.full_reptend` computes multiplicative orders, identifies
full-reptend primes up to a finite limit, and constructs repetends such as
`1/7 = 0.(142857)`. The companion tests verify exact modular arithmetic and
the cyclic multiples of 142857. These computations illustrate the finite
group `(Z/pZ)^×`; they do not establish that infinitely many full-reptend
primes exist. That infinitude for base 10 is an instance of Artin's primitive
root conjecture. Hooley proved the corresponding result assuming an appropriate
Generalized Riemann Hypothesis (GRH); RH alone is not the claim being tested
here, and neither conjecture is proved by a finite search.

The Hamming tests also enumerate all 16 codewords and check closure under all
seven cyclic shifts. A separate test demonstrates that one particular
shortening loses six-bit cyclic closure. “Semi-cyclic” is not a single
universal property; shortened Reed–Solomon constructions require a precise
definition and separate algebraic tests, so this observation should not be
generalized to every shortened code.

A longer, status-aware treatment of cyclic numbers, Artin's conjecture, GRH,
cyclic codes, and discrete-logarithm cryptography is in
[`docs/cyclic-numbers-rh-fec-cryptography.md`](docs/cyclic-numbers-rh-fec-cryptography.md).
A companion note on radio modulation — PSK as a cyclic group, space–time group
codes, and Zadoff–Chu sequences — is in
[`docs/cyclic-groups-radio-modulation.md`](docs/cyclic-groups-radio-modulation.md).

### CT-IM: constant-weight index modulation

CT-IM replaces the combinatorial index lookup of OFDM-IM with a
constant-weight cyclic codebook. For `n=7`, the seven weight-3 codewords of
Hamming(7,4) form one cyclic orbit, have minimum distance 4, and correct one
index-bit error. They carry only `log2(7) ≈ 2.807` bits, so a 3-bit label
reserves one word; with QPSK this is 9 bits per 7-subcarrier block
(`1.286` bits/subcarrier) versus `1.571` for unconstrained OFDM-IM. The cyclic
shift can be chosen to reduce PAPR, but the receiver must be told the shift:
every shift is a codeword, so the syndrome cannot identify it. The full
theorem, proof, rate accounting, and limitations are in
[`docs/ct-im.md`](docs/ct-im.md).

### Residual-index densities

`modcode.number_theory.residual_index` computes the GRH density for primes with
residual index `m = (p-1)/ord_p(10)`. The Hooley–Lenstra field-degree sum is
`sum_n mu(n) / [Q(zeta_{mn}, 10^(1/(mn))):Q]`; for base 10 the degree halves
when `40 | mn` because `sqrt(10)` lies in the cyclotomic field. The small cases
are `m=1`: `A`; `m=2`: `3A/4`; `m=4`: `29A/152`. The often-quoted `29A/52` does
not match the formula or the direct counts. `examples/residual_index_lab.py`
plots the predicted densities, the finite-count convergence, and the observed
histogram.

### PSK and Zadoff–Chu waveforms

`modcode.waveforms` exposes `psk_constellation(order, phase_offset=0)`,
`zadoff_chu(length, root=1)`, `periodic_autocorrelation`, and
`cyclic_cross_correlation`. With zero phase offset, the `M` PSK points are the
`M`-th roots of unity and form the cyclic group `C_M`; a global offset gives a
coset, not a subgroup. ZC sequences use the quadratic phase
`exp(-j*pi*root*n*(n+1)/length)` and, for odd length with a coprime root, have
zero periodic autocorrelation at every nonzero shift. The ZC root index is a
phase coefficient, not a primitive root. The checks live in
`tests/test_waveforms.py`.

### CRC and quasi-cyclic post-quantum structures

`modcode.fec.crc` provides systematic CRC-8/16/32 encoding and detection. A
CRC appends the remainder of `x**r * m(x)` modulo a generator polynomial
`g(x)`; it is a cyclic code used for detection rather than correction.
`division_trace` exposes the shift-register states for plotting.

`modcode.pqc.QCMDPC` is a deliberately tiny QC-MDPC McEliece-style
demonstration. Its parity check is `H = [H0 | H1]` with circulant blocks, so
the public `Q` block is determined by one row (`n` bits) instead of a dense
`n x n` matrix (`n^2` bits). The same circulant structure used for QC-LDPC
error correction is therefore reused for key-size reduction in code-based
post-quantum cryptography. The security assumption is the hardness of decoding
a random quasi-cyclic code; the module’s parameters and toy decoder are **not
secure** and are included only to make the structure visible.

BIKE is the canonical QC-MDPC example; HQC is another quasi-cyclic code-based
KEM, but its construction is not simply “QC-LDPC,” so the shorthand should not
be read as a precise classification. `modcode.pqc.analysis` adds the
security experiments plotted in Figure 13: bit-flipping success versus error
weight, the exponential syndrome search space, and the low-weight dual-code
spectrum that key-recovery attacks target. The model is a diagnostic toy, not
a security proof or a parameter recommendation.

### Optical IM/DD experiments

`modcode.optical` provides OOK, conventional nonnegative PAM, a small
noise-margin search for a non-uniform 4-PAM candidate, and a symbol-rate
IM/DD impairment model. The Matplotlib helpers display intensity levels and
received symbols, and compare BER for configured cases:

```python
from modcode.optical import (
    IMDDChannel,
    IntensityPAM,
    OOK,
    ber_sweep_imdd,
    optimize_pam_levels,
)
from modcode.viz.optical import plot_intensity_constellations, plot_robustness

candidate = optimize_pam_levels(noise_std=0.05, amplitude_noise_std=0.1)
formats = {"OOK": OOK(), "PAM4": IntensityPAM(4), "candidate": candidate}
channels = {
    "noise": IMDDChannel(noise_std=0.05),
    "combined": IMDDChannel(
        noise_std=0.05,
        amplitude_noise_std=0.1,
        nonlinearity=0.05,
        bandwidth=0.05,
    ),
}
results = ber_sweep_imdd(formats, channels, n_bits=100_000, seed=2026)
plot_intensity_constellations(formats)
plot_robustness(results)
```

The candidate search maximizes a worst-neighbor margin under its specified
additive and multiplicative noise assumptions while fixing peak and average
intensity. It is an experiment, not a claim that the scheme is novel or
universally better; assess it against realistic optical link models and
measured BER, power efficiency, bandwidth, and implementation complexity.
The bandwidth impairment is currently only a three-tap symbol-rate ISI proxy.
This real-valued IM/DD model has no optical phase, so it does not model phase
noise or coherent optical transmission.

### APSK and cyclic structure in radio systems

An APSK constellation is arranged on concentric rings. The equally spaced
phase positions on one ring can be indexed by the cyclic group `C_n`, but a
multi-ring APSK constellation is not itself generally a group under complex
multiplication: ring radii and the union of ring points need not be closed
under that operation. `modcode.modulation.APSK` exposes this precisely as a
rotation action on each ring and tests 16-APSK with four inner and twelve
outer points. Ring geometry alone also does not guarantee better noise or
nonlinear-distortion performance than QAM; performance depends on ring radii,
point labeling, coding, channel, and amplifier constraints.

### Phased arrays, cyclic phase resolution, and quantum frequency conversion

`modcode.phasedarray` synthesises the scalar field of a uniform linear array
from the cylindrical-wave sum, not just the array-factor approximation, so the
near-field wavefront curvature is real. A phase ramp can be applied exactly, or
spent on a cyclic group `C_n` in two different ways:

- **Per-element rounding** keeps the average ramp slope, so the beam stays near
  the commanded angle, but the aperture acquires a staircase aberration. The
  Marechal estimate `exp(-rms^2)` predicts the resulting scan loss.
- **Step rounding** quantises the inter-element step once, keeping a perfectly
  flat wavefront but aiming the beam at the quantised angle. At 17 degrees a
  `C4` array keeps full gain and points at 30 degrees, 13 degrees off.

Those are separate hardware choices, not one error, and the module exposes
both. `array_factor` also accepts a frequency ratio, so wideband squint
(`sin(theta) = sin(theta0) / (f/f0)`) and the grating-lobe spacing limit are
numerical rather than asserted.

`modcode.qfc` carries the transduction side: photon energy and Bose-Einstein
occupation, a phase-insensitive converter with calibrated versus uncalibrated
efficiency, Gaussian-state covariance propagation through that converter, and
a transport comparison between coax, waveguide, free-space RF/optical, and
fibre. The converter is usually not what destroys a microwave state. A 4 K bath
at 10 GHz already carries 7.8 photons per mode, so the state is thermal before
conversion begins; quantum frequency conversion is an argument about transport
and detectors, not a way around cryogenics. Calibrated loss is likewise not
excess noise, which is why fibre CV-QKD survives tens of dB of it.

```python
import numpy as np

from modcode.phasedarray import LinearArray
from modcode.qfc import Transducer, thermal_occupation, transport_budget

array = LinearArray(elements=16, frequency=10e9)
phases = array.element_quantized_phases(17.0, order=8)
print(array.peak_angle(phases), array.scan_loss_db(17.0, phases))

device = Transducer(efficiency=0.5, added_noise=1e-3)
print(device.excess_noise(calibrated=True))
print(device.surviving_squeezing_db(10.0, bath_occupation=thermal_occupation(10e9, 4.0)))
for path in transport_budget(100e3):
    print(path.label, path.loss_db)
```

Run `python examples/wavefront_qfc_lab.py` to regenerate five figures in
`examples/figures/`: wavefront steering, phase quantisation, wideband squint,
the QFC transport/noise comparison, and the Gaussian-state survival plot. The
script also prints the numeric report that the test suite pins down.

### Tail-biting trellis coding and cyclic shifts

`TailBitingTCM8PSK` is a compact rate-2/3 trellis-coded example mapped to
8-PSK. Its four-state trellis is initialized with its final input pair
(tail-biting), and its transition rule uses the preceding pair cyclically.
Tests check both Viterbi recovery and equivariance: rotating the input-pair
sequence rotates its encoded symbol sequence. This is a specific constructed
example—not a claim that all trellis-coded modulation is cyclically closed,
that any arbitrary symbol rotation remains a valid path, or that cyclicity
automatically reduces decoder complexity.

### Probabilistic shaping and sequence symmetry

`modcode.shaping` computes Maxwell-Boltzmann symbol probabilities, samples a
shaped sequence, and computes periodic autocorrelation. Shaping changes
symbol probabilities; it does not by itself impose invariance on codewords,
FEC, or sequence correlations. A cyclic shift preserves periodic
autocorrelation, while an arbitrary permutation generally does not, as the
tests illustrate. The linked 2026 preprint is titled
[“How Does Probabilistic Constellation Shaping Control Ranging Sidelobe Levels
in ISAC?”](https://arxiv.org/abs/2610.21653); its topic is shaping and
ambiguity/ranging sidelobes, not a general theorem that FEC requires
permutation invariance. PAS systems need compatible distribution matching,
coding, and bit-metric decoding, but there is no universal cyclic-balance
condition that guarantees error-free correction.

Run `python examples/explore_cyclic_structures.py` to see the 16-APSK
constellation and a sample of its shaped symbol distribution. APSK's rings
are not a standard 5G NR modulation format; the project treats this as an
independent research experiment.

5G NR uses CP-OFDM with QPSK and specified QAM orders for the modulation
formats in its radio interface; APSK is not a standard NR modulation format.
The cyclic prefix is a copy of the end of an OFDM block. Under the usual
assumptions that the channel memory fits within the prefix and synchronization
is maintained, removing the prefix makes the useful received block behave as
circular convolution, which the DFT diagonalizes into per-subcarrier channel
gains. The prefix does not turn an arbitrary physical channel into a cyclic
one without those conditions. NR uses LDPC for shared-channel data and Polar
codes for specified control channels; quasi-cyclic parity-check structures
are an implementation technique, distinct from every codeword being closed
under arbitrary cyclic shifts.

### Quasi-cyclic LDPC matrix structure

`modcode.fec.qc_ldpc.expand_qc_base_matrix` expands a compact shift-value
matrix into a binary parity-check matrix. Each `-1` becomes a zero block;
each shift `s` becomes `I` rotated right by `s` columns modulo lifting factor
`Z`. For example, `[[0, 1], [2, -1]]` at `Z=3` expands into 6-by-6 blocks
containing `I`, `P`, `P²`, and a zero block. The tests verify the exact block
entries and the one-`1`-per-row-and-column permutation property. This is a
generic QC-LDPC construction for studying the matrix pattern, not a complete
5G NR LDPC encoder or the standardized NR base graphs and shift tables.
Run `python examples/qc_ldpc_shift_demo.py` for a 5-by-5 shift-by-2 block and
an illustrative lifted base matrix. The example is not a standardized 5G
base-graph fragment.

### Polar codes and the Shannon limit

`modcode.fec.polar.PolarBEC` implements an Arikan polar transform, selects
information positions using exact binary-erasure-channel (BEC) synthetic
channel probabilities, and decodes with successive cancellation (SC). For
BEC erasure probability `epsilon`, one polarization step produces
`epsilon_minus = 2*epsilon - epsilon**2` and `epsilon_plus = epsilon**2`.
Their average capacity is preserved:

```text
((1 - epsilon_minus) + (1 - epsilon_plus)) / 2 = 1 - epsilon
```

As the block length grows, the synthetic BEC channels polarize toward
nearly-perfect and nearly-useless channels. For the BEC, the fraction of
nearly-perfect channels approaches capacity `1 - epsilon`; selecting those
channels is the asymptotic reason Polar codes can approach Shannon capacity
with SC decoding. "Frozen" means the encoder sets selected source-bit
positions to known values, usually zero; it does not mean a physical channel
has become noiseless. SC estimates information bits sequentially using
previous decisions, and those decisions can be wrong. No finite blocklength
reaches capacity with zero error: the code has a finite-length rate/reliability
tradeoff. The included code is
an educational BEC implementation, not a 5G NR Polar encoder: it omits
rate-matching, CRC-aided list decoding, and construction for general channels.

Run `python examples/explore_qc_polar.py` to print an expanded lifted matrix,
the exact finite-length BEC synthetic erasure probabilities, the selected
information indices, and the capacity identity. For example, at erasure
probability `0.4`, the Shannon capacity is `0.6` bits/use; the example selects
9 of 16 bit-channels (rate `0.5625`) using their BEC reliability estimates,
then runs a seeded example frame through the encoder, BEC, and SC decoder.

`modcode.cryptography.toy_dh` demonstrates Diffie–Hellman in the cyclic
multiplicative group modulo a small prime and tests that both participants
derive the same value. It is intentionally insecure and for mathematical
experimentation only. Real cryptography depends on careful parameter selection,
protocol design, side-channel resistance, and reviewed implementations—not
merely on choosing a group that is non-cyclic. Many deployed systems use
cyclic subgroups; proposing a “stronger” group requires a precisely defined
hard problem and cryptanalysis, not just a different group structure.

## Figures

All figures are stored in the repository. The maintained figures are
regenerated by the example scripts: wavefront/QFC by
`python examples/wavefront_qfc_lab.py`, CRC and QC-MDPC by
`python examples/cyclic_codes_pqc_lab.py`, residual-index densities by
`python examples/residual_index_lab.py`, and CT-IM by
`python examples/ct_im_lab.py`. The legacy figures are retained from the
original scratch scripts under `src/modcode/cyclic-5g/`. Click a figure to
open the full PNG.

### Maintained wavefront and QFC lab

| Wavefront steering | Phase quantisation |
|---|---|
| <a href="examples/figures/wavefront_steering.png"><img src="examples/figures/wavefront_steering.png" alt="Wavefronts for exact phase ramps at 0, 15, 30, and 45 degrees" width="440"></a> | <a href="examples/figures/phase_quantization.png"><img src="examples/figures/phase_quantization.png" alt="Per-element and uniform-step C4, C8, and C16 phase quantisation at 17 degrees" width="440"></a> |
| <sub><b>Figure 1.</b> Near-field wavefronts and phase maps for exact ramps at 0°, 15°, 30°, and 45° (16 elements, 10 GHz, d = λ/2).</sub> | <sub><b>Figure 2.</b> C4/C8/C16 quantisation at 17°: per-element rounding keeps the pointing but adds wavefront ripple; uniform-step rounding keeps a flat wavefront but moves the beam.</sub> |

| Wideband squint | QFC transport and noise |
|---|---|
| <a href="examples/figures/beam_squint.png"><img src="examples/figures/beam_squint.png" alt="Realised pointing versus drive frequency and patterns at three frequencies" width="440"></a> | <a href="examples/figures/qfc_transport.png"><img src="examples/figures/qfc_transport.png" alt="Transport loss, thermal occupation, and transduction excess-noise budget" width="440"></a> |
| <sub><b>Figure 3.</b> A fixed phase ramp is pointed only at its design frequency; a 10% shift moves the beam about 3° at 30°.</sub> | <sub><b>Figure 4.</b> 100 km transport loss, thermal photons per mode, and the calibrated versus uncalibrated excess-noise budget.</sub> |

<p align="center">
  <a href="examples/figures/qfc_state.png"><img src="examples/figures/qfc_state.png" alt="Squeezed-state covariance through the converter and surviving squeezing versus efficiency and bath temperature" width="700"></a>
</p>
<p align="center"><sub><b>Figure 5.</b> A squeezed state through the converter: Wigner contours, surviving squeezing versus efficiency, and the 4 K bath that erases 10 dB of microwave squeezing.</sub></p>

### Legacy cyclic-5g outputs

These figures are retained from the original scratch lab. The maintained
figures above are authoritative; the scratch scripts are kept for provenance
and the `cyclic-5g` directory is not an importable package because its name
contains a hyphen.

| Legacy wavefront | Legacy quantisation |
|---|---|
| <a href="src/modcode/cyclic-5g/figures/wavefront_steering.png"><img src="src/modcode/cyclic-5g/figures/wavefront_steering.png" alt="Legacy wavefront figure for 0, 15, 30, and 45 degrees" width="440"></a> | <a href="src/modcode/cyclic-5g/figures/quantization_error.png"><img src="src/modcode/cyclic-5g/figures/quantization_error.png" alt="Legacy quantisation comparison at 17 degrees" width="440"></a> |
| <sub><b>Figure 6.</b> Original scratch-lab wavefront figure for 0°, 15°, 30°, and 45°.</sub> | <sub><b>Figure 7.</b> Original scratch-lab C4/C8 quantisation comparison at 17°.</sub> |

| Legacy QFC transduction | Legacy beam patterns |
|---|---|
| <a href="src/modcode/cyclic-5g/figures/qfc_transduction.png"><img src="src/modcode/cyclic-5g/figures/qfc_transduction.png" alt="Legacy QFC transport and noise figure" width="440"></a> | <a href="src/modcode/cyclic-5g/figures/beam_patterns.png"><img src="src/modcode/cyclic-5g/figures/beam_patterns.png" alt="Legacy polar beam patterns, quantisation damage, and phase ramps" width="440"></a> |
| <sub><b>Figure 8.</b> Original scratch-lab QFC transport and transduction-noise figure.</sub> | <sub><b>Figure 9.</b> Original scratch-lab polar patterns, quantisation damage, and phase ramps.</sub> |

### Cyclic coding and post-quantum structure

| CRC detection | QC-MDPC matrix |
|---|---|
| <a href="examples/figures/crc_encoding.png"><img src="examples/figures/crc_encoding.png" alt="CRC-8 codeword with an error, division register trace, and single-bit error syndromes" width="440"></a> | <a href="examples/figures/qc_mdpc_matrix.png"><img src="examples/figures/qc_mdpc_matrix.png" alt="QC-MDPC parity check with two circulant blocks, one block, and key-size scaling" width="440"></a> |
| <sub><b>Figure 10.</b> CRC-8 frame with one error, the remainder-register trace, and the nonzero remainder left by every single-bit error.</sub> | <sub><b>Figure 11.</b> The toy QC-MDPC check `H = [H0 | H1]`, the circulant block `H0`, and the linear key-size scaling from `n` bits instead of `n²`.</sub> |

<p align="center">
  <a href="examples/figures/qc_mdpc_decoding.png"><img src="examples/figures/qc_mdpc_decoding.png" alt="Injected QC-MDPC error, syndrome before decoding, and bit-flipping syndrome weight" width="900"></a>
</p>
<p align="center"><sub><b>Figure 12.</b> A weight-two QC-MDPC error, its syndrome, and the syndrome weight falling to zero during bit-flipping.</sub></p>

<p align="center">
  <a href="examples/figures/qc_mdpc_security.png"><img src="examples/figures/qc_mdpc_security.png" alt="QC-MDPC decoder success versus error weight, log2 syndrome search space, and dual-code low-weight spectrum" width="900"></a>
</p>
<p align="center"><sub><b>Figure 13.</b> The security experiments: bit-flipping success falls as the error weight grows, the syndrome search space grows exponentially, and the dual code contains the sparse secret rows at weight 6 (toy n=15).</sub></p>

<p align="center">
  <a href="examples/figures/residual_index_density.png"><img src="examples/figures/residual_index_density.png" alt="Predicted residual-index densities, finite-count convergence, and observed histogram" width="900"></a>
</p>
<p align="center"><sub><b>Figure 14.</b> Residual-index densities for base 10: predicted GRH values, finite-count convergence for m=1,2,4, and the observed histogram at X=1,000,000.</sub></p>

<p align="center">
  <a href="examples/figures/ct_im_codebook.png"><img src="examples/figures/ct_im_codebook.png" alt="CT-IM constant-weight codebook, pairwise Hamming distance matrix, and rate comparison" width="900"></a>
</p>
<p align="center"><sub><b>Figure 15.</b> The CT-IM codebook: the seven weight-3 Hamming(7,4) words, their pairwise distances (4 off the diagonal), and the rate cost of index error correction versus classical OFDM-IM.</sub></p>

## Module reference

| Module | Purpose | Main entry points |
|---|---|---|
| `modcode` | Core convenience exports | `BPSK`, `QPSK`, `QAM`, `Hamming74`, `AWGNChannel`, `Link`, `LinkResult`, `run_link` |
| `modcode.modulation` | Modulation families and the shared interface, plus CT-IM | `BPSK`, `QPSK`, `QAM`, `APSK`, `TailBitingTCM8PSK`, `CTIM`, `ConstantWeightCyclicCode`, `Modulation` |
| `modcode.waveforms` | PSK and Zadoff–Chu generation, periodic correlation | `psk_constellation`, `zadoff_chu`, `periodic_autocorrelation`, `cyclic_cross_correlation` |
| `modcode.fec` | Block-code interface, CRC, Hamming, and Polar | `FEC`, `CRC`, `Hamming74`, `PolarBEC`; `modcode.fec.cyclic` (`cyclic_shift`), `modcode.fec.qc_ldpc` (`expand_qc_base_matrix`) |
| `modcode.pqc` | Toy QC-MDPC McEliece-style structure and security experiments | `QCMDPC`, `DecodeResult`, `circulant`, `is_circulant`, `decoder_success_rates`, `syndrome_search_space_bits`, `dual_weight_distribution` |
| `modcode.channel` | Channel models | `AWGNChannel` |
| `modcode.pipeline` | End-to-end link simulation | `Link`, `LinkResult`, `run_link` |
| `modcode.metrics` | BER measurement | `ber_sweep` |
| `modcode.viz` | Constellation, BER, array, CRC, PQC, CT-IM, and QFC plotting | `plot_constellation`, `plot_ber`, `plot_wavefront`, `plot_pattern`, `plot_phase_ramp`, `plot_squint`, `plot_crc_frame`, `plot_crc_division`, `plot_crc_syndromes`, `plot_parity_check`, `plot_circulant_block`, `plot_key_scaling`, `plot_decoder_success`, `plot_search_space`, `plot_dual_weight_histogram`, `plot_codebook`, `plot_distance_matrix`, `plot_rate_budget`, `plot_transport_budget`, `plot_noise_budget`, `plot_wigner` |
| `modcode.phasedarray` | Array geometry and cyclic phase quantisation | `LinearArray`, `quantize_phases`, `squint_angle`, `grating_lobe_limit` |
| `modcode.qfc` | QFC noise, quantum states, and link budgets | `Transducer`, `thermal_occupation`, `transport_budget`, `wigner_function`, `covariance_squeezing_db` |
| `modcode.optical` | IM/DD experiments | `OOK`, `IntensityPAM`, `NonUniformPAM`, `IMDDChannel`, `optimize_pam_levels`, `ber_sweep_imdd` |
| `modcode.number_theory` | Full-reptend arithmetic and residual-index densities | `multiplicative_order`, `is_full_reptend_prime`, `full_reptend_cycle`, `cyclic_multiples`, `full_reptend_primes`, `artin_constant`, `residual_index_density`, `count_residual_indices` |
| `modcode.cryptography` | Deliberately insecure Diffie–Hellman | `ToyDiffieHellman` |
| `modcode.shaping` | Probabilistic shaping and autocorrelation | `maxwell_boltzmann_probabilities`, `sample_shaped_symbols`, `periodic_autocorrelation` |
| `src/modcode/cyclic-5g` | Standalone scratch lab; not importable | `qkd_phased_array_lab.py`, `wavefront_qfc_lab.py` |

## Documentation

The README is the map; the long-form notes are in `docs/`:

- [Cyclic numbers, the Riemann Hypothesis, cyclic codes, and cryptography](docs/cyclic-numbers-rh-fec-cryptography.md)
  — full-reptend primes, Artin/GRH, quasi-RH caveats, residual-index densities,
  cyclic codes, DLP, quasi-cyclic post-quantum codes, quantum cyclic topology,
  and a claims audit.
- [Cyclic groups in radio: PSK, space–time codes, and Zadoff–Chu](docs/cyclic-groups-radio-modulation.md)
  — PSK as `C_M`, unitary group codes, ZC autocorrelation, the KHAN proposal,
  and runnable checks backed by `modcode.waveforms`.
- [CT-IM: constant-weight cyclic index modulation](docs/ct-im.md)
  — the constant-weight theorem and proof, the Hamming(7,4) orbit, rate
  accounting, cyclic-shift detection limits, and literature positioning.

## Project structure

```text
modcode-lab/
├── .github/workflows/ci.yml
├── docs/
│   ├── ct-im.md
│   ├── cyclic-groups-radio-modulation.md
│   └── cyclic-numbers-rh-fec-cryptography.md
├── examples/
│   ├── figures/
│   │   ├── beam_squint.png
│   │   ├── crc_encoding.png
│   │   ├── ct_im_codebook.png
│   │   ├── phase_quantization.png
│   │   ├── qc_mdpc_decoding.png
│   │   ├── qc_mdpc_matrix.png
│   │   ├── qc_mdpc_security.png
│   │   ├── qfc_state.png
│   │   ├── qfc_transport.png
│   │   ├── residual_index_density.png
│   │   └── wavefront_steering.png
│   ├── cyclic_codes_pqc_lab.py
│   ├── ct_im_lab.py
│   ├── explore_cyclic_structures.py
│   ├── explore_qc_polar.py
│   ├── hamming_demo.py
│   ├── qc_ldpc_shift_demo.py
│   ├── residual_index_lab.py
│   ├── run_ber_sweep.py
│   └── wavefront_qfc_lab.py
├── src/modcode/
│   ├── __init__.py
│   ├── shaping.py
│   ├── channel/
│   │   ├── __init__.py
│   │   ├── awgn.py
│   │   └── base.py
│   ├── cryptography/
│   │   ├── __init__.py
│   │   └── toy_dh.py
│   ├── cyclic-5g/                  # scratch lab; not an importable package
│   │   ├── figures/
│   │   │   ├── beam_patterns.png
│   │   │   ├── qfc_transduction.png
│   │   │   ├── quantization_error.png
│   │   │   └── wavefront_steering.png
│   │   ├── qkd_phased_array_lab.py
│   │   ├── run_output.txt
│   │   └── wavefront_qfc_lab.py
│   ├── fec/
│   │   ├── __init__.py
│   │   ├── base.py
│   │   ├── crc.py
│   │   ├── cyclic.py
│   │   ├── hamming.py
│   │   ├── polar.py
│   │   └── qc_ldpc.py
│   ├── metrics/
│   │   ├── __init__.py
│   │   └── ber.py
│   ├── modulation/
│   │   ├── __init__.py
│   │   ├── apsk.py
│   │   ├── base.py
│   │   ├── bpsk.py
│   │   ├── ct_im.py
│   │   ├── cyclic_tcm.py
│   │   ├── qam.py
│   │   └── qpsk.py
│   ├── number_theory/
│   │   ├── __init__.py
│   │   ├── full_reptend.py
│   │   └── residual_index.py
│   ├── optical/
│   │   ├── __init__.py
│   │   ├── channel.py
│   │   ├── modulation.py
│   │   └── robustness.py
│   ├── phasedarray/
│   │   ├── __init__.py
│   │   └── array.py
│   ├── pipeline/
│   │   ├── __init__.py
│   │   └── chain.py
│   ├── pqc/
│   │   ├── __init__.py
│   │   ├── analysis.py
│   │   └── qc_mdpc.py
│   ├── qfc/
│   │   ├── __init__.py
│   │   ├── link.py
│   │   └── photons.py
│   ├── viz/
│   │   ├── __init__.py
│   │   ├── ber_plot.py
│   │   ├── constellation.py
│   │   ├── crc.py
│   │   ├── ct_im.py
│   │   ├── number_theory.py
│   │   ├── optical.py
│   │   ├── phasedarray.py
│   │   ├── pqc.py
│   │   └── qfc.py
│   └── waveforms/
│       ├── __init__.py
│       └── sequences.py
├── tests/
│   ├── test_apsk.py
│   ├── test_ber_awgn.py
│   ├── test_crc.py
│   ├── test_ct_im.py
│   ├── test_cyclic_fec.py
│   ├── test_cyclic_tcm.py
│   ├── test_full_reptend.py
│   ├── test_hamming.py
│   ├── test_optical_imdd.py
│   ├── test_phasedarray.py
│   ├── test_polar_bec.py
│   ├── test_pqc_analysis.py
│   ├── test_qam_gray.py
│   ├── test_qc_ldpc.py
│   ├── test_qc_mdpc.py
│   ├── test_qfc.py
│   ├── test_residual_index.py
│   ├── test_shaping.py
│   ├── test_toy_dh.py
│   └── test_waveforms.py
├── .gitignore
├── LICENSE
├── README.md
└── pyproject.toml
```

`src/modcode/cyclic-5g/` is kept as provenance for the original QKD and
phased-array scratch scripts. It is not importable because the directory name
contains a hyphen; the maintained implementations are `modcode.phasedarray`
and `modcode.qfc`, exercised by `examples/wavefront_qfc_lab.py`.

## Testing and reproducibility

```bash
python -m pytest          # full test suite
ruff check .              # lint
python examples/wavefront_qfc_lab.py     # wavefront/QFC figures
python examples/cyclic_codes_pqc_lab.py # CRC and QC-MDPC figures
```

The suite covers modulation geometry, Gray labeling, AWGN BER, CRC detection,
Hamming and QC-LDPC structure, Polar BEC construction, optical IM/DD,
full-reptend arithmetic, residual-index densities and counts, toy
Diffie–Hellman, phased-array synthesis, QFC noise, PSK/ZC waveforms, the
toy QC-MDPC key structure, decoder threshold, syndrome search space, and
dual-code weight spectrum, plus the CT-IM constant-weight codebook, distance
structure, index-error correction, and rate accounting. The figure-producing
scripts are deterministic given their fixed seeds and constants; the numeric
reports they print are the same ones the tests pin down.

## Extending modcode-lab

Implement `Modulation` to provide `bits_per_symbol`, `constellation`,
`modulate`, and `demodulate`. Implement `FEC` with its rate, input and encoded
block sizes, `encode`, and `decode` methods to use a new block code in `Link`.
Channels implement the `Channel.transmit` interface. The shared interfaces
make components composable without changing the simulation pipeline.

## Roadmap

| Milestone | Status | Content |
|---|---|---|
| **M0** | ✅ | Skeleton, BPSK, AWGN, BER sweep, CI |
| **M1** | ✅ | QPSK/QAM, Hamming(7,4), examples |
| **M2** | ⬜ | APSK, Rayleigh fading, multipath echo channel |
| **M3** | ⬜ | BCH, Reed-Solomon, convolutional codes |
| **M4** | ⬜ | TCM / C-TCM over 8-PSK |
| **M5** | ⬜ | Documentation site (MkDocs), Jupyter notebooks |

## Contributing

Create a virtual environment, install `.[dev]`, then run `ruff check .` and
`pytest` from the repository root. Contributions should include focused tests
for new modulation, coding, and channel behavior.

## References

- Proakis and Salehi, *Digital Communications*, 5th edition.
- Sklar, *Digital Communications: Fundamentals and Applications*, 2nd edition.
- Lin and Costello, *Error Control Coding*, 2nd edition.
- W. W. Peterson and D. T. Brown, “Cyclic Codes for Error Detection,”
  *Proceedings of the IRE* **49** (1961), 228–235.
- R. Misoczki, J.-P. Tillich, N. Sendrier, and P. S. L. M. Barreto,
  “MDPC-McEliece: New McEliece variants from moderate density parity-check
  codes,” *Proc. IEEE ISIT* (2013).
- BIKE team, *BIKE: Bit Flipping Key Encapsulation*, specification,
  <https://bikesuite.org/>.
- P. Moree, “Artin's primitive root conjecture — a survey,” *Integers* **12**
  (2012), no. 6.
- E. Başar, Ü. Aygölü, E. Panayırcı, and H. V. Poor, “Orthogonal
  frequency-division multiplexing with index modulation,” *IEEE Transactions
  on Signal Processing* **61** (2013), no. 22, 5536–5549.
- E. Başar, “Index modulation techniques for 5G wireless networks,” *IEEE
  Communications Magazine* **54** (2016), no. 9, 168–175.
- I. A. Simakov and I. S. Besedin, “Low-overhead quantum error correction
  codes with a cyclic topology,” *Physical Review A* **111**, 012444 (2025);
  [arXiv:2211.03094](https://arxiv.org/abs/2211.03094).

## License

This project is licensed under the GNU General Public License, version 3 or
(at your option) any later version. See [LICENSE](LICENSE).