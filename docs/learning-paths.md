# Learning paths

`modcode-lab` has a stable core and a set of exploratory appendices.  The
four routes below are the shortest paths through the material.  Each route is a
sequence, not a syllabus; run the commands from the repository root.

If you only have ten minutes, read the route headings and pick one.  If you
are new to the project, start with the core route that matches your question,
not with the appendices.

---

## Path 1 — "I want to understand soft decoding"

**Goal:** see why reliability matters, then see a soft decoder turn that
reliability into a real coding gain.

| Step | Do this | What to notice |
|---:|---|---|
| 1 | `python examples/soft_decision_lab.py` | Hamming hard versus soft, then convolutional hard versus soft Viterbi. The soft Viterbi gain is the main event. |
| 2 | Read [`tests/test_soft_decision.py`](../tests/test_soft_decision.py) | The LLR sign convention, Hamming soft ML, and the closed-form BER values. |
| 3 | Read [the LLR section](core-simulator.md#21-llr-convention) | What `soft_demodulate` returns, why positive means bit 0, and why max-log is used. |
| 4 | Read [`tests/test_convolutional.py`](../tests/test_convolutional.py) | The noiseless round trip, a single-error correction, and the fixed-seed comparison that requires soft Viterbi to beat hard Viterbi. |

**You should be able to explain:** why a hard decision throws away reliability;
why Viterbi on LLRs beats Viterbi on sliced bits; what the sign of an LLR
means; and why the code rate can make low-SNR coded performance worse than
uncoded even when the decoder is correct.

**Checkpoint:** `python -m pytest tests/test_soft_decision.py
tests/test_convolutional.py -q` should pass.

**Next, the modern codes:** Path 3 below is the guided route through LDPC and
polar: `python examples/ldpc_lab.py` → [`ldpc.md`](ldpc.md) →
`python examples/polar_lab.py` → [`polar.md`](polar.md).

**Optional FEC detour (10 minutes):** `python examples/hamming_demo.py` →
[`tests/test_hamming.py`](../tests/test_hamming.py) → [the Hamming
section](core-simulator.md#3-hamming74). This shows the block-code side before
the soft-decision path.

---

## Path 2 — "I want to understand modulation"

**Goal:** connect bits to symbols, see the theoretical BER reference, and then
see a trellis code that spends redundancy on the modulation itself.

| Step | Do this | What to notice |
|---:|---|---|
| 1 | `python examples/run_ber_sweep.py` | BPSK, QPSK, and 16-QAM simulated curves sit on the closed-form references. The script prints counts and Wilson intervals. |
| 2 | Read [the LLR section](core-simulator.md#21-llr-convention) | The sign convention, the max-log formula, and the per-axis QAM decomposition. |
| 3 | `python examples/tcm_lab.py` | Four-state 8-PSK TCM versus uncoded QPSK at the same 2 bit/symbol. The TCM curve crosses QPSK around 3–4 dB. |
| 4 | Read [the TCM section](advanced-coding.md#5-tcm-8-psk) | The parity equation, the four-state trellis, soft Viterbi, and the honest low-SNR caveat. |

**Optional constellation detour:** [`tests/test_qam_gray.py`](../tests/test_qam_gray.py)
pins Gray adjacency, unit average power, round trip, and the large-batch
slicer.

**You should be able to explain:** why BPSK and Gray QPSK have the same BER
per bit; what `soft_demodulate` returns and why the sign matters; why 16-QAM
needs more Eb/N0 for the same BER; and why TCM is a joint coding/modulation
design rather than a code followed by a mapper.

**Checkpoint:** `test_bpsk_ber_at_6db_matches_published_value` in
[`tests/test_ber_awgn.py`](../tests/test_ber_awgn.py) pins BPSK at 6 dB to the
published value $2.388\times10^{-3}$, and the 8 dB test pins the high-SNR
tail.

---

## Path 3 — "I want to understand the modern codes (LDPC and polar)"

**Goal:** see the two code families used in 5G, measure which one wins at a
given block length, and see what list decoding and CRC-aided selection buy.

| Step | Do this | What to notice |
|---:|---|---|
| 1 | `python examples/ldpc_lab.py` | LDPC(1008,506) min-sum versus the K=3 convolutional code: a steep waterfall and a crossover between 1 and 2 dB. |
| 2 | Read [`ldpc.md`](ldpc.md) | The quasi-cyclic construction, rank 502, the min-sum update, and the measured gain of roughly 1.5–2 dB at $10^{-3}$. |
| 3 | `python examples/polar_lab.py` | SC → SCL → CA-SCL at $N=1024$: the list gain, its saturation at 2 dB, and the measured decoder cost per block. |
| 4 | Read [`polar.md`](polar.md) | Gaussian-approximation construction, the exact check update, the LDPC/polar crossover near 1.2 dB, and block-length/design-point sensitivity. |
| 5 | Read [`tests/test_ldpc.py`](../tests/test_ldpc.py) | (3,6) regularity, rank 502, and $G H^T = 0$. |
| 6 | Read [`tests/test_polar_awgn.py`](../tests/test_polar_awgn.py) | GA numerics, agreement with an independent SC decoder, and the $N\log_2 N$ work pin. |

**You should be able to explain:** why LDPC and polar are the two 5G families;
which one wins at $N \approx 1000$ and rate $1/2$ in this lab (the measured
crossover is near 1.2 dB); why SCL is roughly a dB better than SC at
$10^{-3}$; why the list saturates at high SNR; and why none of the codes here
is 3GPP conformant.

**Checkpoint:** `python -m pytest tests/test_ldpc.py tests/test_polar_awgn.py
-q` should pass, and `examples/figures/polar.png` (Figure 24 in
[`figures.md`](figures.md)) should show the crossover, not just one family.

---

## Path 4 — "I want to understand the appendices"

**Goal:** enter one exploratory area through its own note and runnable
example, without mistaking it for the core.

Start with [`exploratory-labs.md`](exploratory-labs.md).  It is the index and
the scope statement.  Then pick one row:

| Area | Start with | Runnable |
|---|---|---|
| CT-IM constant-weight index modulation | [`ct-im.md`](ct-im.md) | `python examples/ct_im_lab.py` |
| I/Q phase involution | [`phase-shift-involution.md`](phase-shift-involution.md) | `python examples/phase_imbalance_lab.py` |
| Lattice modulation (D4/E8) | [`lattice-modulation.md`](lattice-modulation.md) | `python examples/lattice_lab.py` |
| Sphere packing and the 1.53 dB ceiling | [`sphere-packing-modulation.md`](sphere-packing-modulation.md) | `python examples/sphere_packing_lab.py` |
| Tau-invariant RM(1,3) | [`tau-invariant.md`](tau-invariant.md) | `python examples/tau_invariant_lab.py` |
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
  overlays.  The core path is Paths 1, 2, and 3.
- **Coding extensions:** interleaving, Rayleigh, BCH/RS, TCM, LDPC, polar, and
  the OFDM multipath/Doppler path with pilot-based channel estimation.  These
  are documented in [`advanced-coding.md`](advanced-coding.md),
  [`ldpc.md`](ldpc.md), [`polar.md`](polar.md), and
  [`multipath.md`](multipath.md); they are still pedagogical, but they are part
  of the core-facing material.
- **Appendices:** optical, arrays/QFC, PQC, number theory, CT-IM, phase
  involution.  These are Path 4 and are not required for the core claims.

The core theory is in [`core-simulator.md`](core-simulator.md); the core
quickstart is in the [README](../README.md#quickstart).  The figure index is
in [`figures.md`](figures.md).
