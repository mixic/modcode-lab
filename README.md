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
- [Project structure](#project-structure)
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
| **Visualization** | Constellation diagrams with bit labels, BER-vs-Eb/N0 plots |
| **Tests** | pytest suite verifying Gray adjacency, Hamming syndrome table, BER monotonicity |
| **CI** | GitHub Actions on Python 3.10 / 3.11 / 3.12 |

### Experimental research prototypes

| Area | Prototype |
|---|---|
| **Cyclic geometry** | APSK rings with explicit cyclic phase actions |
| **Trellis coding** | Tail-biting rate-2/3 8-PSK example with pair-rotation equivariance |
| **Probabilistic shaping** | Maxwell-Boltzmann probabilities and periodic autocorrelation |
| **Structured FEC** | QC-LDPC base-matrix lifting and Polar coding over the BEC |

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

## Project structure

```text
src/modcode/
├── modulation/       # BPSK, QPSK, QAM, APSK, tail-biting 8-PSK TCM
├── fec/              # Hamming, QC-LDPC matrix, Polar BEC, and FEC interface
├── channel/          # Channel interface and Eb/N0-normalized AWGN
├── pipeline/         # Link and run_link
├── metrics/          # BER sweeps
├── viz/              # Constellation and BER plotting
├── number_theory/    # Multiplicative orders and full-reptend exploration
├── optical/          # IM/DD modulation and impairment experiments
├── cryptography/     # Insecure educational Diffie-Hellman demonstration
└── shaping.py        # Probabilistic symbol shaping and periodic autocorrelation
examples/              # BER, coding, cyclicity, QC-LDPC, and Polar demos
tests/                 # Signal, coding, matrix, polarization, and cyclicity tests
```

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

## License

This project is licensed under the GNU General Public License, version 3 or
(at your option) any later version. See [LICENSE](LICENSE).