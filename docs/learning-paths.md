# Learning paths

`modcode-lab` has a stable core and a set of exploratory appendices.  The
three routes below are the shortest paths through the material.  Each route is
a sequence, not a syllabus; run the commands from the repository root.

If you only have ten minutes, read the route headings and pick one.  If you
are new to the project, start with the core route that matches your question,
not with the appendices.

---

## Path 1 — "I want to understand FEC"

**Goal:** see a code correct an error, then see why soft information changes
the decoder, then see a trellis code with real coding gain.

| Step | Do this | What to notice |
|---:|---|---|
| 1 | `python examples/hamming_demo.py` | A 4-bit message becomes a 7-bit codeword, one bit is flipped, and the syndrome decoder recovers the message. |
| 2 | Read [`tests/test_hamming.py`](../tests/test_hamming.py) | The test pins the systematic codeword and every single-bit correction. The codebook is small enough to inspect. |
| 3 | Read [the Hamming section](core-simulator.md#3-hamming74) | Why the cyclic code has distance 3, how the syndrome maps to a position, and how soft ML uses all 16 codewords. |
| 4 | `python examples/soft_decision_lab.py` | Hamming hard versus soft, then convolutional hard versus soft Viterbi. The soft Viterbi gain is the main event. |
| 5 | Read [`tests/test_convolutional.py`](../tests/test_convolutional.py) | The noiseless round trip, a single-error correction, and the fixed-seed comparison that requires soft Viterbi to beat hard Viterbi. |

**You should be able to explain:** why a code with distance 3 corrects one
error; why a hard decision throws away reliability; why Viterbi on LLRs beats
Viterbi on sliced bits; and why the code rate makes low-SNR coded performance
worse than uncoded even when the decoder is correct.

**Checkpoint:** `python -m pytest tests/test_hamming.py
tests/test_convolutional.py -q` should pass.

---

## Path 2 — "I want to understand modulation"

**Goal:** connect bits to symbols, see the theoretical BER reference, and then
see a trellis code that spends redundancy on the modulation itself.

| Step | Do this | What to notice |
|---:|---|---|
| 1 | `python examples/run_ber_sweep.py` | BPSK, QPSK, and 16-QAM simulated curves sit on the closed-form references. The script prints counts and Wilson intervals. |
| 2 | Read [`tests/test_qam_gray.py`](../tests/test_qam_gray.py) | Gray adjacency, unit average power, round trip, and the large-batch slicer check. |
| 3 | Read [the LLR section](core-simulator.md#21-llr-convention) | The sign convention, the max-log formula, and the per-axis QAM decomposition. |
| 4 | `python examples/tcm_lab.py` | Four-state 8-PSK TCM versus uncoded QPSK at the same 2 bit/symbol. The TCM curve crosses QPSK around 3–4 dB. |
| 5 | Read [the TCM section](advanced-coding.md#5-tcm-8-psk) | The parity equation, the four-state trellis, soft Viterbi, and the honest low-SNR caveat. |

**You should be able to explain:** why BPSK and Gray QPSK have the same BER
per bit; what `soft_demodulate` returns and why the sign matters; why 16-QAM
needs more Eb/N0 for the same BER; and why TCM is a joint coding/modulation
design rather than a code followed by a mapper.

**Checkpoint:** the 4 dB BPSK simulated and theoretical values in
`run_ber_sweep.py` should agree within Monte Carlo error.

---

## Path 3 — "I want to understand the appendices"

**Goal:** enter one exploratory area through its own note and runnable
example, without mistaking it for the core.

Start with [`exploratory-labs.md`](exploratory-labs.md).  It is the index and
the scope statement.  Then pick one row:

| Area | Start with | Runnable |
|---|---|---|
| CT-IM constant-weight index modulation | [`ct-im.md`](ct-im.md) | `python examples/ct_im_lab.py` |
| I/Q phase involution | [`phase-shift-involution.md`](phase-shift-involution.md) | `python examples/phase_imbalance_lab.py` |
| Phased arrays and QFC | `exploratory-labs.md` §6 | `python examples/wavefront_qfc_lab.py` |
| CRC and toy QC-MDPC | [`cyclic-numbers-rh-fec-cryptography.md`](cyclic-numbers-rh-fec-cryptography.md) | `python examples/cyclic_codes_pqc_lab.py` |
| Number theory / residual indices | [`cyclic-numbers-rh-fec-cryptography.md`](cyclic-numbers-rh-fec-cryptography.md) | `python examples/residual_index_lab.py` |
| PSK, space–time codes, Zadoff–Chu | [`cyclic-groups-radio-modulation.md`](cyclic-groups-radio-modulation.md) | `tests/test_waveforms.py` |
| Optical IM/DD | `exploratory-labs.md` §3 | `tests/test_optical_imdd.py` |

**You should be able to explain:** what the appendix demonstrates, what it
does not demonstrate, and why it is separated from the core.  Each note has a
scope box; read that before the equations.

**Checkpoint:** pick one appendix, run its example, and state one limitation
from its note.  If you cannot state a limitation, the appendix has not been
read as an appendix.

---

## How the paths relate to the repository

- **Core:** modulation, FEC, AWGN/known-CSI fading, LLRs, decoding, BER, theory
  overlays.  The core path is Paths 1 and 2.
- **Coding extensions:** interleaving, Rayleigh, BCH/RS, TCM.  These are
  documented in [`advanced-coding.md`](advanced-coding.md) and are still
  pedagogical, but they are part of the core-facing material.
- **Appendices:** optical, arrays/QFC, PQC, number theory, CT-IM, phase
  involution.  These are Path 3 and are not required for the core claims.

The core theory is in [`core-simulator.md`](core-simulator.md); the core
quickstart is in the [README](../README.md#quickstart).  The figure index is
in [`figures.md`](figures.md).
