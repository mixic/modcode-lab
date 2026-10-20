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

`modcode-lab` is a small simulator for the **core** of a digital communication
link:

- map bits to BPSK, QPSK, or square QAM symbols;
- optionally protect them with a block or convolutional code;
- pass them through AWGN;
- decode with hard decisions or **soft LLRs**;
- measure BER and compare it with the closed-form uncoded reference.

The core is deliberately narrow: modulation, coding, channel, LLRs, decoding,
and BER.  The repository also contains **exploratory appendices** — optical
IM/DD, phased arrays, quantum frequency conversion, post-quantum toy codes,
number theory, CT-IM, and an I/Q phase-involution note.  They are worked
examples, not part of the core simulator and not required for the core BER
claims.  See [`docs/exploratory-labs.md`](docs/exploratory-labs.md).

If you want the theory behind the core, read
[`docs/core-simulator.md`](docs/core-simulator.md).  If you want the quickest
path to a BER curve, read the quickstart below.

---

## Features

### Core simulator

| Category | What is implemented |
|---|---|
| **Modulation** | BPSK, Gray-mapped QPSK, square Gray-mapped QAM (4, 16, 64, 256, ...); unit average energy |
| **Soft demodulation** | Max-log bit LLRs for every finite constellation; positive means bit 0 |
| **FEC** | Hamming(7,4) with syndrome decoding and soft ML over 16 codewords; rate-1/2 constraint-length-3 convolutional code with hard and soft Viterbi |
| **Channels** | AWGN with information-bit Eb/N0 normalization |
| **Pipeline** | `Link` with `decoding="auto"`, `"hard"`, or `"soft"`; `LinkResult.decoding` records the path |
| **Metrics** | Simulated BER and closed-form BPSK/QPSK/QAM references |
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
| CT-IM | `examples/ct_im_lab.py` | [`ct-im.md`](docs/ct-im.md) |
| I/Q phase involution | `examples/phase_imbalance_lab.py` | [`phase-shift-involution.md`](docs/phase-shift-involution.md) |

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

### 3. Run the examples

```bash
# core
python examples/run_ber_sweep.py       # BER + theory overlays
python examples/soft_decision_lab.py   # Hamming/convolutional hard vs soft
python examples/hamming_demo.py

# exploratory appendices
python examples/explore_cyclic_structures.py
python examples/explore_qc_polar.py
python examples/wavefront_qfc_lab.py
python examples/cyclic_codes_pqc_lab.py
python examples/residual_index_lab.py
python examples/ct_im_lab.py
python examples/phase_imbalance_lab.py
```

---

## Core module reference

| Module | Purpose | Main entry points |
|---|---|---|
| `modcode` | Convenience exports | `BPSK`, `QPSK`, `QAM`, `Hamming74`, `ConvolutionalCode`, `AWGNChannel`, `Link`, `run_link` |
| `modcode.modulation` | Constellations and LLRs | `BPSK`, `QPSK`, `QAM`, `Modulation.soft_demodulate`, `demodulate_phase_imbalance` |
| `modcode.fec` | FEC codes | `Hamming74` (hard + `decode_soft`), `ConvolutionalCode` (hard/soft Viterbi), `CRC`, `PolarBEC` |
| `modcode.channel` | Core channel | `AWGNChannel`; `PhaseImbalanceChannel` is an exploratory appendix |
| `modcode.pipeline` | Link composition | `Link(..., decoding="auto"/"hard"/"soft")`, `LinkResult`, `run_link` |
| `modcode.metrics` | BER and theory | `ber_sweep`, `theoretical_ber`, `q_function` |
| `modcode.viz` | Plots | `plot_constellation`, `plot_ber(..., theory=...)` |

The full API and derivations are in [`docs/core-simulator.md`](docs/core-simulator.md).

---

## Documentation

| Document | Contents |
|---|---|
| [`core-simulator.md`](docs/core-simulator.md) | Signal model, LLR convention, Hamming hard/soft, convolutional trellis/Viterbi, pipeline modes, theoretical BER |
| [`exploratory-labs.md`](docs/exploratory-labs.md) | One place for the appendices: optical, arrays/QFC, PQC, number theory, CT-IM, phase involution |
| [`figures.md`](docs/figures.md) | Core BER figures and the complete exploratory gallery |
| [`ct-im.md`](docs/ct-im.md) | CT-IM constant-weight note and milestone status |
| [`phase-shift-involution.md`](docs/phase-shift-involution.md) | The DSIN τ / I/Q correspondence and its limits |
| [`cyclic-groups-radio-modulation.md`](docs/cyclic-groups-radio-modulation.md) | PSK as a cyclic group, space–time codes, ZC sequences |
| [`cyclic-numbers-rh-fec-cryptography.md`](docs/cyclic-numbers-rh-fec-cryptography.md) | Full-reptend primes, Artin/GRH, cyclic codes, DLP, QC-MDPC, QEC |

---

## Project structure

```text
modcode-lab/
├── docs/
│   ├── core-simulator.md
│   ├── exploratory-labs.md
│   ├── figures.md
│   ├── ct-im.md
│   ├── phase-shift-involution.md
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

- `tests/test_ber_awgn.py` — noise normalization and link BER;
- `tests/test_soft_decision.py` — LLR sign/scale, Hamming soft ML, theory values;
- `tests/test_convolutional.py` — hard/soft Viterbi and the soft-decoding gain;
- `tests/test_hamming.py`, `tests/test_qam_gray.py` — code and constellation structure.

The core figures are generated by:

```bash
python examples/run_ber_sweep.py
python examples/soft_decision_lab.py
```

Monte Carlo points are not pinned by tests; the codebooks, LLR conventions,
rates, and closed-form values are.

---

## Roadmap

The core roadmap is intentionally short.  The appendices are not on it.

| Priority | Status | Item |
|---|---|---|
| P0 | ✅ | Soft LLR demodulation and soft Hamming ML |
| P0 | ✅ | Rate-1/2 convolutional code with hard/soft Viterbi |
| P0 | ✅ | Theoretical BER overlays for BPSK/QPSK/QAM |
| P1 | ⬜ | Rayleigh/block-fading channel and CSI assumptions |
| P1 | ⬜ | BCH and Reed–Solomon codes |
| P2 | ⬜ | TCM / C-TCM over 8-PSK |
| P2 | ⬜ | Documentation site, notebooks, and more channel models |

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
