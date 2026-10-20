# Figures

All figures are stored in the repository.  The core BER figures are generated
by the two core examples; the exploratory figures are kept for the appendices
described in [`exploratory-labs.md`](exploratory-labs.md).

Click a figure to open the full PNG.

---

## Core simulator

| BER with theoretical references | Hard versus soft decoding |
|---|---|
| <a href="../examples/figures/ber_sweep.png"><img src="../examples/figures/ber_sweep.png" alt="BPSK, QPSK, and 16-QAM BER versus Eb/N0 with Q-function theory overlays" width="440"></a> | <a href="../examples/figures/soft_decision.png"><img src="../examples/figures/soft_decision.png" alt="Uncoded BPSK theory, Hamming hard/soft, and convolutional hard/soft BER" width="440"></a> |
| <sub><b>Figure 1.</b> Uncoded BPSK, QPSK, and 16-QAM against the closed-form references. The dashed curves are theory, not fits.</sub> | <sub><b>Figure 2.</b> Hamming(7,4) and rate-1/2 convolutional decoding: hard decisions versus LLRs. The convolutional soft curve is the main coding gain.</sub> |

---

## Radio, waveforms, and trellis examples

| Wideband squint | Phase quantisation |
|---|---|
| <a href="../examples/figures/beam_squint.png"><img src="../examples/figures/beam_squint.png" alt="Beam squint versus frequency" width="440"></a> | <a href="../examples/figures/phase_quantization.png"><img src="../examples/figures/phase_quantization.png" alt="Per-element and uniform-step phase quantisation" width="440"></a> |
| <sub><b>Figure 3.</b> Wideband squint: a fixed phase ramp points at its design frequency only.</sub> | <sub><b>Figure 4.</b> C4/C8/C16 phase quantisation at 17 degrees: per-element rounding versus a uniform step.</sub> |

The APSK/TCM/shaping examples are code-first and do not currently emit a
dedicated figure; see `examples/explore_cyclic_structures.py`.

---

## Phased arrays and quantum frequency conversion

| Wavefront steering | QFC transport and noise |
|---|---|
| <a href="../examples/figures/wavefront_steering.png"><img src="../examples/figures/wavefront_steering.png" alt="Wavefronts for exact phase ramps" width="440"></a> | <a href="../examples/figures/qfc_transport.png"><img src="../examples/figures/qfc_transport.png" alt="Transport loss and thermal occupation" width="440"></a> |
| <sub><b>Figure 5.</b> Near-field wavefronts for exact phase ramps on a linear array.</sub> | <sub><b>Figure 6.</b> 100 km transport loss, thermal photons, and the transduction excess-noise budget.</sub> |

<p align="center">
  <a href="../examples/figures/qfc_state.png"><img src="../examples/figures/qfc_state.png" alt="Squeezed-state covariance and surviving squeezing" width="700"></a>
</p>
<p align="center"><sub><b>Figure 7.</b> A squeezed state through the converter: Wigner contours, surviving squeezing, and the thermal floor.</sub></p>

---

## Coding, number theory, and exploratory index modulation

| CRC detection | QC-MDPC matrix |
|---|---|
| <a href="../examples/figures/crc_encoding.png"><img src="../examples/figures/crc_encoding.png" alt="CRC-8 frame, division trace, and single-bit syndromes" width="440"></a> | <a href="../examples/figures/qc_mdpc_matrix.png"><img src="../examples/figures/qc_mdpc_matrix.png" alt="QC-MDPC circulant blocks and key-size scaling" width="440"></a> |
| <sub><b>Figure 8.</b> CRC-8 frame with one error, the remainder-register trace, and the nonzero remainder left by each single-bit error.</sub> | <sub><b>Figure 9.</b> The toy QC-MDPC check `H = [H0 | H1]` and the linear key-size scaling of a circulant block.</sub> |

<p align="center">
  <a href="../examples/figures/qc_mdpc_decoding.png"><img src="../examples/figures/qc_mdpc_decoding.png" alt="QC-MDPC error, syndrome, and bit-flipping weight" width="900"></a>
</p>
<p align="center"><sub><b>Figure 10.</b> A weight-two QC-MDPC error and the syndrome weight falling to zero during bit-flipping.</sub></p>

<p align="center">
  <a href="../examples/figures/qc_mdpc_security.png"><img src="../examples/figures/qc_mdpc_security.png" alt="QC-MDPC decoder success, search space, and dual weights" width="900"></a>
</p>
<p align="center"><sub><b>Figure 11.</b> Toy QC-MDPC security experiments: decoder success, syndrome search space, and dual-code low-weight spectrum.</sub></p>

<p align="center">
  <a href="../examples/figures/residual_index_density.png"><img src="../examples/figures/residual_index_density.png" alt="Residual-index densities and convergence" width="900"></a>
</p>
<p align="center"><sub><b>Figure 12.</b> Residual-index densities for base 10: GRH predictions, finite-count convergence, and the observed histogram.</sub></p>

<p align="center">
  <a href="../examples/figures/ct_im_codebook.png"><img src="../examples/figures/ct_im_codebook.png" alt="CT-IM codebook, distance matrix, and rate comparison" width="900"></a>
</p>
<p align="center"><sub><b>Figure 13.</b> The CT-IM weight-3 Hamming orbit, pairwise distances, and the rate cost against classical OFDM-IM.</sub></p>

<p align="center">
  <a href="../examples/figures/phase_imbalance.png"><img src="../examples/figures/phase_imbalance.png" alt="QPSK reflection at phi=pi, BER, channel determinant, and blind LLR" width="900"></a>
</p>
<p align="center"><sub><b>Figure 14.</b> The I/Q reflection at phi = pi, hard/oracle/pilot/blind BER, the channel determinant, and the identically zero blind LLR.</sub></p>

---

## Legacy scratch-lab outputs

These are retained for provenance from `src/modcode/cyclic-5g/`; the
maintained implementations are `modcode.phasedarray` and `modcode.qfc`.

| Legacy wavefront | Legacy quantisation |
|---|---|
| <a href="../src/modcode/cyclic-5g/figures/wavefront_steering.png"><img src="../src/modcode/cyclic-5g/figures/wavefront_steering.png" alt="Legacy wavefront figure" width="440"></a> | <a href="../src/modcode/cyclic-5g/figures/quantization_error.png"><img src="../src/modcode/cyclic-5g/figures/quantization_error.png" alt="Legacy quantisation figure" width="440"></a> |
| <sub><b>Figure 15.</b> Original scratch-lab wavefront figure.</sub> | <sub><b>Figure 16.</b> Original scratch-lab C4/C8 quantisation comparison.</sub> |

| Legacy QFC | Legacy beam patterns |
|---|---|
| <a href="../src/modcode/cyclic-5g/figures/qfc_transduction.png"><img src="../src/modcode/cyclic-5g/figures/qfc_transduction.png" alt="Legacy QFC figure" width="440"></a> | <a href="../src/modcode/cyclic-5g/figures/beam_patterns.png"><img src="../src/modcode/cyclic-5g/figures/beam_patterns.png" alt="Legacy beam patterns" width="440"></a> |
| <sub><b>Figure 17.</b> Original QFC transport and noise figure.</sub> | <sub><b>Figure 18.</b> Original polar beam patterns and phase ramps.</sub> |

---

## Regenerate

```bash
python examples/run_ber_sweep.py
python examples/soft_decision_lab.py
python examples/wavefront_qfc_lab.py
python examples/cyclic_codes_pqc_lab.py
python examples/residual_index_lab.py
python examples/ct_im_lab.py
python examples/phase_imbalance_lab.py
```
