# Exploratory labs and appendices

The core of `modcode-lab` is the modulation/FEC/BER simulator described in
[`core-simulator.md`](core-simulator.md).  The modules below are **appendices**:
they are useful, tested demonstrations, but they do not define the project and
they are not required for the core BER claims.

They are kept because they connect the core to adjacent questions: cyclic
structure in radio, optical IM/DD, post-quantum code structure, number theory,
phased arrays, quantum transduction, and index modulation.  Treat them as
worked examples, not as a claim that all of these topics belong in one
research programme.

---

## What is where

| Area | Code | Runnable example | Note |
|---|---|---|---|
| PSK and Zadoff–Chu waveforms | `modcode.waveforms` | `tests/test_waveforms.py` | [`cyclic-groups-radio-modulation.md`](cyclic-groups-radio-modulation.md) |
| APSK, tail-biting 8-PSK, shaping | `modcode.modulation`, `modcode.shaping` | `examples/explore_cyclic_structures.py` | cyclic-group note |
| Optical IM/DD | `modcode.optical` | `tests/test_optical_imdd.py` | no separate note |
| CRC and cyclic codes | `modcode.fec.crc` | `examples/cyclic_codes_pqc_lab.py` | [`cyclic-numbers-rh-fec-cryptography.md`](cyclic-numbers-rh-fec-cryptography.md) |
| QC-LDPC and Polar | `modcode.fec.qc_ldpc`, `modcode.fec.polar` | `examples/explore_qc_polar.py`, `examples/qc_ldpc_shift_demo.py` | cyclic-numbers note |
| Toy QC-MDPC post-quantum structure | `modcode.pqc` | `examples/cyclic_codes_pqc_lab.py` | cyclic-numbers note §3.7 |
| Full-reptend primes and residual indices | `modcode.number_theory` | `examples/residual_index_lab.py` | cyclic-numbers note |
| Phased arrays and QFC | `modcode.phasedarray`, `modcode.qfc` | `examples/wavefront_qfc_lab.py` | no separate note |
| CT-IM constant-weight index modulation | `modcode.modulation.ct_im` | `examples/ct_im_lab.py` | [`ct-im.md`](ct-im.md) |
| I/Q phase involution | `modcode.channel.phase_imbalance` | `examples/phase_imbalance_lab.py` | [`phase-shift-involution.md`](phase-shift-involution.md) |
| Lattice modulation (D4/E8) | `modcode.modulation.lattice` | `examples/lattice_lab.py` | [`lattice-modulation.md`](lattice-modulation.md) |
| Tau-invariant RM(1,3) | `modcode.fec.tau_invariant` | `examples/tau_invariant_lab.py` | [`tau-invariant.md`](tau-invariant.md) |
| Legacy cyclic-5g scratch lab | `src/modcode/cyclic-5g` | not importable | provenance only |

---

## 1. Waveforms and cyclic radio structure

`modcode.waveforms` generates PSK constellations and Zadoff–Chu sequences and
computes periodic autocorrelation and cyclic cross-correlation.  With zero
phase offset, the $M$ PSK points are the $M$-th roots of unity and form the
cyclic group $C_M$; a global offset gives a coset, not a subgroup.

The companion note develops PSK as a cyclic group, unitary space–time group
codes, ZC autocorrelation, and the KHAN proposal.  It is a radio-modulation
note, not a claim about the core BER simulator.

See [`cyclic-groups-radio-modulation.md`](cyclic-groups-radio-modulation.md)
and `tests/test_waveforms.py`.

---

## 2. APSK, tail-biting 8-PSK, and shaping

`APSK` exposes each ring as a cyclic phase orbit.  The multi-ring
constellation is **not** asserted to be a group under complex multiplication;
the code documents the ring action separately.

`TailBitingTCM8PSK` is a small rate-2/3 8-PSK example with pair-rotation
equivariance.  It is a demonstration of tail-biting trellis structure, not a
competitive TCM design.  The core four-state TCM example is `TCM8PSK`,
documented in [`advanced-coding.md`](advanced-coding.md); the tail-biting
variant remains an appendix.

`modcode.shaping` implements Maxwell–Boltzmann probabilities and periodic
autocorrelation for shaped symbols.

See `examples/explore_cyclic_structures.py` and `tests/test_apsk.py`,
`tests/test_cyclic_tcm.py`, `tests/test_shaping.py`.

---

## 3. Optical IM/DD

`modcode.optical` compares OOK, intensity PAM, and nonnegative PAM over an
IM/DD channel, including clipping and level optimization.  This is a separate
intensity-modulation model; it does not reuse the complex-baseband BER
normalization without saying so.

See `tests/test_optical_imdd.py`.

---

## 4. CRC, QC-LDPC, Polar, and post-quantum structure

`modcode.fec.crc` implements systematic CRC-8/16/32 and exposes the division
trace.  CRC is a detection code, not a correction code.

`modcode.fec.qc_ldpc` demonstrates quasi-cyclic lifting from a base matrix;
`modcode.fec.polar` also contains `PolarBEC`, which constructs a Polar code for
the binary erasure channel. Both are structural demonstrations, not decoded BER
comparisons against the core codes. The **core** LDPC decoder is `LDPC36`,
documented in [`ldpc.md`](ldpc.md); the **core** BI-AWGN polar decoder is
`PolarAWGN`, documented in [`polar.md`](polar.md). Neither is part of this
appendix.

`modcode.pqc.QCMDPC` is a deliberately tiny QC-MDPC McEliece-style toy.  It
shows circulant blocks and linear key-size scaling; it is **not cryptography**
and must not be used for security.  The analysis functions measure the decoder
threshold, the syndrome search space, and the dual-code weight spectrum at toy
parameters.

See [`cyclic-numbers-rh-fec-cryptography.md`](cyclic-numbers-rh-fec-cryptography.md)
and `examples/cyclic_codes_pqc_lab.py`.

---

## 5. Number theory and residual indices

`modcode.number_theory` computes multiplicative orders, full-reptend primes,
cyclic multiples, and the GRH residual-index densities.  The finite searches
do not prove Artin's conjecture or RH; the note records the distinction.

The often-quoted $m=4$ density $29A/52$ is not supported by the formula or
the direct counts; the implemented value is $29A/152$.

See [`cyclic-numbers-rh-fec-cryptography.md`](cyclic-numbers-rh-fec-cryptography.md)
and `examples/residual_index_lab.py`.

---

## 6. Phased arrays and quantum frequency conversion

`modcode.phasedarray` synthesizes linear-array wavefronts, quantizes phases,
and computes squint and grating-lobe limits.  `modcode.qfc` computes transport
loss, thermal occupation, transduction noise, squeezing survival, and Wigner
functions.

These are self-contained physics demonstrations.  They are not part of the
digital link pipeline and do not share its BER normalization.

See `examples/wavefront_qfc_lab.py`.

---

## 7. CT-IM: constant-weight index modulation

CT-IM replaces the combinatorial OFDM-IM index lookup with the seven weight-3
codewords of Hamming(7,4).  The codebook is one cyclic orbit, has minimum
distance 4, and corrects one index-bit error.  The cost is a smaller index
alphabet: $\log_2 7\approx2.807$ bits, or 9 bits per 7-subcarrier block with
QPSK versus 11 for classical OFDM-IM.

This is an **exploratory appendix**, not a core result.  The implemented
object is the unpunctured weight-3 orbit; the punctured variant, baselines, BER
and PAPR experiments, Rayleigh extension, and BCH(15,7) variant are not done.
The blind cyclic-shift theorem is an impossibility result: every shift is a
codeword, so the syndrome cannot identify it.

See [`ct-im.md`](ct-im.md) and `examples/ct_im_lab.py`.

---

## 8. The I/Q phase-shift involution

The DSIN relative-phase flip $\tau=\operatorname{diag}(I,-I)$ and the I/Q
phase error at $\phi=\pi$ share the real matrix $\operatorname{diag}(1,-1)$.
The correspondence is basis-dependent and is not an identity between the two
models.

The engineering result is an **identifiability** statement: for a symmetric
constellation and prior, the blind hypotheses $\phi=0$ and $\phi=\pi$ have
identical likelihoods.  A pilot, preamble, asymmetric constellation, or code
constraint can break the tie.  This is an interpretive appendix, not a new
modulation scheme and not evidence for RH.

See [`phase-shift-involution.md`](phase-shift-involution.md) and
`examples/phase_imbalance_lab.py`.

---

## 9. Lattice modulation (D4 and E8)

`LatticeModulation` builds spherical-cut D4 and E8 constellations and compares
them with QAM at matched bits and energy per **real dimension**.  The
dimension matching matters: D4-16 is comparable to QPSK, not 16-QAM.

The tested claim is a minimum-distance claim: D4-256 beats 16-QAM and E8-256
beats QPSK at matched $\eta$.  The BER measurement is more nuanced: D4-16
does not beat QPSK, E8-256 wins above roughly 5–6 dB, and D4-256 crosses
16-QAM only near 11–12 dB.  The component-stream PAPR diagnostic and the naive
software demodulation cost are reported as costs, not as a system-level
verdict.

This is a standalone appendix; multi-dimensional symbols are not carried by
the scalar-complex `Link` pipeline.

See [`lattice-modulation.md`](lattice-modulation.md) and
`examples/lattice_lab.py`.

---

## 10. Tau-invariant codes (RM(1,3))

`TauInvariantCode` implements a code invariant under a coordinate involution.
RM(1,3) is the canonical example: length 8, dimension 4, distance 4, invariant
under the non-cyclic coordinate flip $i\mapsto i\oplus1$, and not cyclic.
It uses a hard nearest-codeword decoder and a soft ML decoder over the 16
codewords.  The comparison with cyclic Hamming(7,4) is structural, not a gain
claim.

See [`tau-invariant.md`](tau-invariant.md) and `examples/tau_invariant_lab.py`.

---

## 11. Provenance: the cyclic-5g scratch lab

`src/modcode/cyclic-5g/` is retained as provenance for the original QKD and
phased-array scratch scripts.  It is not importable because the directory name
contains a hyphen.  The maintained implementations are `modcode.phasedarray`
and `modcode.qfc`; use `examples/wavefront_qfc_lab.py`, not the scratch files.

---

## What the appendices are not

- They are not a claim that the project's scope is “all of cyclic mathematics”.
- They are not a security proof, an RH proof, or a QKD implementation.
- They are not on the core roadmap and do not block core releases.
- Their figures are kept in [`figures.md`](figures.md), separate from the
  quickstart.
