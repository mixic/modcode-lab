<div align="center">

# modcode-lab

**A reproducible Python lab for digital modulation, forward error correction, and link-level BER.**

*Bits → FEC → modulation → channel → LLR or hard decision → FEC decode → BER.*

[![CI](https://github.com/mixic/modcode-lab/actions/workflows/ci.yml/badge.svg)](https://github.com/mixic/modcode-lab/actions/workflows/ci.yml)
[![Python](https://img.shields.io/badge/python-3.10%20%7C%203.11%20%7C%203.12-blue.svg)](https://www.python.org/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Code style: ruff](https://img.shields.io/badge/code%20style-ruff-000000.svg)](https://github.com/astral-sh/ruff)

</div>

---

## What this is

Want to see why soft decoding beats hard decoding by 2 dB, or how Hamming(7,4)
turns a 3-bit syndrome into an error location? This lab reproduces those
results end-to-end, with theory overlays that let you check the answer.

`modcode-lab` is a reproducible lab for the **core** of a digital
communication link:

- map bits to BPSK, QPSK, or square QAM symbols;
- protect them with one of seven families — Hamming(7,4), a K=3 convolutional
  code, BCH(15,7), Reed–Solomon(15,11), LDPC(1008,506), polar(1024,512), or
  four-state 8-PSK TCM;
- pass them through AWGN, block Rayleigh fading, or bursts;
- decode with hard decisions or **soft LLRs**, from syndrome lookup through
  Viterbi to min-sum belief propagation and successive-cancellation list
  decoding;
- measure BER and compare it with the closed-form references.

The core now spans the classical block codes, the convolutional baseline, and
both modern families (LDPC and polar), so the trade-offs are measured in one
pipeline rather than quoted: the LDPC-versus-polar BER and decoding-cost
comparison is Figure 24 in [`docs/figures.md`](docs/figures.md).  The boundary
is still a boundary: these are small pedagogical codes with explicit
[non-claims](docs/polar.md#6-what-is-claimed-measured-and-not-claimed), not
3GPP implementations.

The repository also contains **exploratory appendices** — optical IM/DD,
phased arrays, quantum frequency conversion, post-quantum toy codes, number
theory, CT-IM, and an I/Q phase-involution note.  They are worked examples, not
part of the core simulator and not required for the core BER claims.  See
[`docs/exploratory-labs.md`](docs/exploratory-labs.md).

If you want the theory behind the core, read
[`docs/core-simulator.md`](docs/core-simulator.md).  If you want the quickest
path to a BER curve, read the quickstart below.  If you are new to the
repository, pick one of the [learning paths](docs/learning-paths.md).

---

## Features

### Core simulator

| Category | What is implemented |
|---|---|
| **Modulation** | BPSK, Gray-mapped QPSK, square Gray-mapped QAM (4, 16, 64, 256, ...); unit average energy |
| **Soft demodulation** | Max-log bit LLRs for every finite constellation; positive means bit 0 |
| **FEC** | Hamming(7,4) with syndrome/soft ML; rate-1/2 convolutional code with hard/soft Viterbi; BCH(15,7) hard Berlekamp–Massey + soft ML; RS(15,11) hard Berlekamp–Massey/Chien; LDPC(1008,506) with normalized min-sum belief propagation; Polar(1024,512) with SC, SCL, and CRC-aided list decoding |
| **TCM** | Four-state rate-2/3 8-PSK trellis-coded modulation with soft Viterbi |
| **Interleaving** | Deterministic bit permutation; hard/soft deinterleaving in the link |
| **Channels** | AWGN; block Rayleigh fading with known CSI; synthetic burst channel; OFDM multipath/Doppler with pilot-based LS estimation and per-subcarrier LLRs |
| **Pipeline** | `Link` with `decoding="auto"`, `"hard"`, or `"soft"`; `LinkResult.decoding` records the path |
| **Metrics** | Simulated BER with adaptive stopping and Wilson confidence intervals, plus closed-form BPSK/QPSK/QAM references |
| **Plots** | Constellations and BER curves with dashed theory overlays |
| **Tests** | Full pytest suite, including LLR sign/scale, Hamming soft ML, Viterbi hard/soft, and theory values |
| **CI** | GitHub Actions on Python 3.10 / 3.11 / 3.12 |

### Exploratory appendices

| Area | Example | Note |
|---|---|---|
| PSK and Zadoff–Chu waveforms | `tests/test_waveforms.py` | [`cyclic-groups-radio-modulation.md`](docs/cyclic-groups-radio-modulation.md) |
| APSK, tail-biting 8-PSK, shaping | `examples/explore_cyclic_structures.py` | cyclic-group note |
| Optical IM/DD | `tests/test_optical_imdd.py` | no separate note |
| CRC, QC-LDPC, Polar | `examples/cyclic_codes_pqc_lab.py` | [`cyclic-numbers-rh-fec-cryptography.md`](docs/cyclic-numbers-rh-fec-cryptography.md) |
| Toy QC-MDPC | `examples/cyclic_codes_pqc_lab.py` | cyclic-numbers note §3.7 |
| Number theory / residual indices | `examples/residual_index_lab.py` | cyclic-numbers note |
| Phased arrays and QFC | `examples/wavefront_qfc_lab.py` | no separate note |
| CT-IM | `examples/ct_im_lab.py`, `examples/ct_im_experiments.py` | [`ct-im.md`](docs/ct-im.md) |
| I/Q phase involution | `examples/phase_imbalance_lab.py` | [`phase-shift-involution.md`](docs/phase-shift-involution.md) |
| Lattice modulation (D4/E8) | `examples/lattice_lab.py` | [`lattice-modulation.md`](docs/lattice-modulation.md) |
| Tau-invariant RM(1,3) with FHT | `examples/tau_invariant_lab.py` | [`tau-invariant.md`](docs/tau-invariant.md) |

The appendices are **not** claims that the project is about all of these
topics.  They are linked from one place so the core README stays readable.

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

---

## Quickstart

### 1. Uncoded BER against theory

```python
import numpy as np

from modcode import AWGNChannel, BPSK, Link, run_link
from modcode.metrics import theoretical_ber

ebn0_db = np.array([0.0, 2.0, 4.0, 6.0])
link = Link(BPSK(), AWGNChannel())

measured = np.array(
    [run_link(link, n_bits=200_000, ebn0_db=point).ber for point in ebn0_db]
)
print("measured:", measured)
print("theory:  ", theoretical_ber(BPSK(), ebn0_db))
```

The two arrays should agree within Monte Carlo error.  This is the sanity
check that catches a wrong noise variance or a wrong Eb/N0 convention.

### 2. Coded link: hard versus soft decoding

```python
from modcode import AWGNChannel, BPSK, ConvolutionalCode, Link, run_link

code = ConvolutionalCode(block_bits=64)

hard = run_link(
    Link(BPSK(), AWGNChannel(), code, decoding="hard"),
    n_bits=100_000,
    ebn0_db=3.0,
)
soft = run_link(
    Link(BPSK(), AWGNChannel(), code, decoding="soft"),
    n_bits=100_000,
    ebn0_db=3.0,
)

print(hard.ber, hard.decoding)   # hard
print(soft.ber, soft.decoding)   # soft
```

`decoding="auto"` is the default: it uses LLRs when the channel exposes a
noise variance and the FEC implements `decode_soft`, otherwise it falls back
to hard decisions.  `"soft"` is strict; it raises if soft decoding is not
available.

Hamming(7,4) also supports soft ML decoding:

```python
from modcode import Hamming74

coded = Link(BPSK(), AWGNChannel(), Hamming74(), decoding="soft")
```

BCH and Reed–Solomon use the same interface:

```python
from modcode import BCH15, ReedSolomon

bch = Link(BPSK(), AWGNChannel(), BCH15(), decoding="soft")
rs = Link(BPSK(), AWGNChannel(), ReedSolomon(), decoding="hard")
```

LDPC uses normalized min-sum and the soft path:

```python
from modcode import LDPC36

ldpc = Link(BPSK(), AWGNChannel(), LDPC36(), decoding="soft")
```

Polar uses successive cancellation, list decoding, and an optional CRC:

```python
from modcode import PolarAWGN
from modcode.fec.crc import CRC8

# list_size=1 is plain SC; list_size>1 is SCL
scl = Link(BPSK(), AWGNChannel(), PolarAWGN(1024, 512, list_size=4), decoding="soft")

# CRC-aided SCL: the CRC selects among the surviving list paths
ca_scl = Link(
    BPSK(), AWGNChannel(),
    PolarAWGN(1024, 512, crc=CRC8, list_size=4), decoding="soft",
)
```

An interleaver is selected on the `Link`:

```python
from modcode.pipeline import PermutationInterleaver

coded = Link(
    BPSK(), AWGNChannel(), Hamming74(),
    decoding="soft", interleaver=PermutationInterleaver(seed=7),
)
```

### 3. Run the examples

```bash
# core
python examples/run_ber_sweep.py       # BER + theory overlays
python examples/soft_decision_lab.py   # Hamming, RM(1,3) FHT, convolutional
python examples/interleaver_lab.py     # bursts with/without interleaving
python examples/fading_lab.py          # block Rayleigh + FEC
python examples/bch_rs_lab.py          # BCH, RS, Hamming, convolutional
python examples/ldpc_lab.py            # LDPC min-sum vs convolutional
python examples/polar_lab.py           # polar SC/SCL/CA-SCL vs LDPC
python examples/multipath_lab.py       # OFDM multipath with estimated CSI
python examples/tcm_lab.py             # TCM 8-PSK vs QPSK
python examples/hamming_demo.py

# exploratory appendices
python examples/explore_cyclic_structures.py
python examples/explore_qc_polar.py
python examples/wavefront_qfc_lab.py
python examples/cyclic_codes_pqc_lab.py
python examples/residual_index_lab.py
python examples/ct_im_lab.py
python examples/ct_im_experiments.py
python examples/phase_imbalance_lab.py
python examples/lattice_lab.py
python examples/tau_invariant_lab.py
```

---

## Core module reference

| Module | Purpose | Main entry points |
|---|---|---|
| `modcode` | Convenience exports | `BPSK`, `QPSK`, `QAM`, `BCH15`, `ReedSolomon`, `LDPC36`, `PolarAWGN`, `ConvolutionalCode`, `TCM8PSK`, `Hamming74`, `AWGNChannel`, `BlockRayleighChannel`, `Link`, `run_link` |
| `modcode.modulation` | Constellations, LLRs, TCM | `BPSK`, `QPSK`, `QAM`, `TCM8PSK`, `Modulation.soft_demodulate` |
| `modcode.fec` | FEC codes | `Hamming74`, `ConvolutionalCode`, `BCH15`, `ReedSolomon`, `LDPC36`, `PolarAWGN`, `CRC`, `PolarBEC` |
| `modcode.channel` | Channels | `AWGNChannel`, `BlockRayleighChannel`, `BurstChannel`, `MultipathChannel`/`TapProcess`; `PhaseImbalanceChannel` is an exploratory appendix |
| `modcode.pipeline` | Link composition | `Link(..., decoding="auto"/"hard"/"soft", interleaver=...)`, `PermutationInterleaver`, `LinkResult`, `run_link`; `OFDMLink`/`run_ofdm_link` for multipath with estimated CSI |
| `modcode.metrics` | BER and theory | `ber_sweep`, `ber_sweep_detailed`, `BERCurve`, `theoretical_ber`, `theoretical_rayleigh_ber`, `q_function` |
| `modcode.viz` | Plots | `plot_constellation`, `plot_ber(..., theory=...)` |

The full API and derivations are in [`docs/core-simulator.md`](docs/core-simulator.md).

---

## Documentation

| Document | Contents |
|---|---|
| [`learning-paths.md`](docs/learning-paths.md) | Four guided routes: soft decoding, modulation, modern codes (LDPC and polar), and the appendices |
| [`core-simulator.md`](docs/core-simulator.md) | Signal model, LLR convention, Hamming hard/soft, convolutional trellis/Viterbi, pipeline modes, theoretical BER |
| [`advanced-coding.md`](docs/advanced-coding.md) | Interleaving, block Rayleigh, BCH/RS, TCM, and their limits |
| [`ldpc.md`](docs/ldpc.md) | LDPC(1008,506) construction, min-sum decoding, and finite-length comparison |
| [`polar.md`](docs/polar.md) | Polar(1024,512) construction, SC/SCL/CA-SCL decoding, and the list-size trade-off |
| [`multipath.md`](docs/multipath.md) | Multipath/Doppler channel, pilot-based channel estimation, zero-forcing equalisation, and measured BER |
| [`exploratory-labs.md`](docs/exploratory-labs.md) | One place for the appendices: optical, arrays/QFC, PQC, number theory, CT-IM, phase involution |
| [`figures.md`](docs/figures.md) | Core BER figures and the complete exploratory gallery |
| [`ct-im.md`](docs/ct-im.md) | CT-IM constant-weight note and milestone status |
| [`phase-shift-involution.md`](docs/phase-shift-involution.md) | The DSIN τ / I/Q correspondence and its limits |
| [`lattice-modulation.md`](docs/lattice-modulation.md) | D4/E8 spherical cuts, dimension-matched BER, PAPR diagnostic, and complexity |
| [`tau-invariant.md`](docs/tau-invariant.md) | Convolution versus involution, RM(1,3), FHT decoding, and the non-cyclic coordinate flip |
| [`cyclic-groups-radio-modulation.md`](docs/cyclic-groups-radio-modulation.md) | PSK as a cyclic group, space–time codes, ZC sequences |
| [`cyclic-numbers-rh-fec-cryptography.md`](docs/cyclic-numbers-rh-fec-cryptography.md) | Full-reptend primes, Artin/GRH, cyclic codes, DLP, QC-MDPC, QEC |

---

## Project structure

```text
modcode-lab/
├── docs/
│   ├── learning-paths.md
│   ├── core-simulator.md
│   ├── advanced-coding.md
│   ├── ldpc.md
│   ├── polar.md
│   ├── multipath.md
│   ├── exploratory-labs.md
│   ├── figures.md
│   ├── ct-im.md
│   ├── phase-shift-involution.md
│   ├── lattice-modulation.md
│   ├── tau-invariant.md
│   ├── cyclic-groups-radio-modulation.md
│   └── cyclic-numbers-rh-fec-cryptography.md
├── examples/
│   ├── figures/
│   ├── run_ber_sweep.py
│   ├── soft_decision_lab.py
│   ├── hamming_demo.py
│   └── exploratory examples ...
├── src/modcode/
│   ├── channel/
│   ├── fec/
│   ├── modulation/
│   ├── pipeline/
│   ├── metrics/
│   ├── viz/
│   └── exploratory packages ...
├── tests/
├── pyproject.toml
├── CONTRIBUTING.md
└── README.md
```

The exact file tree is in git; the abbreviated tree above separates the core
from the appendices.

---

## Testing and reproducibility

```bash
python -m pytest          # full suite
ruff check .              # lint
```

Core tests include:

- `tests/test_ber_awgn.py` — noise normalization, Q-function comparison, adaptive stopping, and confidence intervals;
- `tests/test_soft_decision.py` — LLR sign/scale, Hamming soft ML, theory values;
- `tests/test_convolutional.py` — hard/soft Viterbi and the soft-decoding gain;
- `tests/test_interleaver.py`, `tests/test_burst_channel.py` — permutation inverse and burst spreading;
- `tests/test_rayleigh.py` — block fading against the known-CSI reference;
- `tests/test_bch_rs.py`, `tests/test_tcm.py` — BCH/RS correction and TCM gain;
- `tests/test_ldpc.py` — (3,6) regularity, rank, generator orthogonality, and min-sum decoding;
- `tests/test_polar_awgn.py` — GA construction, SC agreement with an independent reference, SCL work scaling, and list-decoding error reduction;
- `tests/test_multipath.py` — tap power/correlation, delay structure, subcarrier gains, and the OFDM estimation path;
- `tests/test_tau_invariant.py` — RM(1,3) invariance, non-cyclicity, ML, and FHT decoding;
- `tests/test_ct_im_experiments.py` — shortened/BCH index codes, ML detection, Rayleigh, and PAPR invariance;
- `tests/test_hamming.py`, `tests/test_qam_gray.py` — code and constellation structure.

The core figures are generated by:

```bash
python examples/run_ber_sweep.py
python examples/soft_decision_lab.py
```

Monte Carlo points are not pinned by tests; the codebooks, LLR conventions,
rates, and closed-form values are.  The 6 dB and 8 dB BPSK cases are pinned
explicitly: 500,000 bits at 6 dB must land within 15% of the published
$2.388\times10^{-3}$, and at 8 dB within 20% of
`0.5 * erfc(sqrt(Eb/N0))`.

---

## Roadmap

The core roadmap is intentionally short.  The appendices are not on it.

| Priority | Status | Item |
|---|---|---|
| P0 | ✅ | Soft LLR demodulation and soft Hamming ML |
| P0 | ✅ | Rate-1/2 convolutional code with hard/soft Viterbi |
| P0 | ✅ | Theoretical BER overlays for BPSK/QPSK/QAM |
| P0 | ✅ | O(N) QAM slicing, vectorized Hamming LUT, adaptive BER with Wilson intervals |
| P0 | ✅ | LDPC(1008,506) with normalized min-sum belief propagation |
| P0 | ✅ | Polar(1024,512) with SC, SCL, and CRC-aided list decoding |
| P1 | ✅ | Block Rayleigh fading with known CSI (hard decoding) |
| P1 | ✅ | Bit interleaving across FEC blocks and QAM symbols |
| P1 | ✅ | BCH(15,7) and Reed–Solomon RS(15,11) codes |
| P2 | ✅ | Four-state 8-PSK TCM example with soft Viterbi |
| P2 | ✅ | Soft LLRs for equalized fading (per-subcarrier variance) in the OFDM path |
| P2 | ✅ | Multipath/Doppler channels with pilot-based channel estimation |
| P2 | ⬜ | Documentation site and notebooks |

Not planned as core: QKD, optical hardware, quantum transduction, post-quantum
security proofs, or a proof of RH.

### Exploratory track

CT-IM and the DSIN/I-Q note are **appendix material**, not core commitments.
They are kept because they are tested and documented, not because they are
results.  The CT-IM milestone table is in
[`ct-im.md` §10](docs/ct-im.md#10-milestones); CT-IM-0, CT-IM-5, and the
initial write-up are complete, the experimental items are open.  The I/Q note
explicitly states that the correspondence is basis-dependent and does not
identify the two models.

No paper-level novelty claim is made for either appendix.

---

## Scope boundaries

- **Pedagogical, not production.** The codes are small; the channel models are
  simple.  This is a lab, not a modem.
- **Toy post-quantum code.** `QCMDPC` is a structure demonstration, not
  cryptography.
- **No security proof.** The QC-MDPC and phase-involution notes describe toy
  models and limits; they do not establish a protocol.
- **No RH claim.** The number-theory note records finite computations and
  conjectures; it does not prove RH or Artin's conjecture.
- **Exploratory is exploratory.** Optical, arrays, QFC, CT-IM, and the I/Q
  involution are appendices; see [`docs/exploratory-labs.md`](docs/exploratory-labs.md).

---

## Contributing

The contribution contract is in [`CONTRIBUTING.md`](CONTRIBUTING.md).  In short:

- a core addition must implement `Modulation`, `FEC`, or `Channel`, work with
  `Link`, and have a focused test that pins a number;
- an appendix must carry a scope box, a tested claim, a bounded non-claim, and a
  runnable example;
- Monte Carlo tests use fixed seeds and a stated tolerance;
- documentation states what is claimed and what is not.

Run `ruff check .` and `pytest` before opening a pull request.  If a change does
not fit one of those two layers, open an issue first.

---

## References

- J. G. Proakis and M. Salehi, *Digital Communications*, 5th edition.
- S. Lin and D. J. Costello, *Error Control Coding*, 2nd edition.
- E. Başar, Ü. Aygölü, E. Panayırcı, and H. V. Poor, “Orthogonal
  frequency-division multiplexing with index modulation,” *IEEE Transactions
  on Signal Processing* **61** (2013), no. 22, 5536–5549.
- E. Başar, “Index modulation techniques for 5G wireless networks,” *IEEE
  Communications Magazine* **54** (2016), no. 9, 168–175.

---

## License

GNU General Public License, version 3 or (at your option) any later version.
See [LICENSE](LICENSE).
