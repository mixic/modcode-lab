# Lattice modulation: D4 and E8 spherical cuts

*Exploratory appendix. This is not part of the core `Link` pipeline and it is
not a claim that lattice modulation replaces QAM.*

Multi-dimensional constellations are a real research direction, but the
comparison has to be dimensionally fair.  This note implements spherical-cut
D4 and E8 constellations, checks their minimum distances, and reports the BER,
PAPR, and complexity costs measured in this repository.

---

## 1. What is implemented

`modcode.modulation.LatticeModulation` supports:

- **D4**: integer points in \(\mathbb R^4\) with even coordinate sum.
- **E8**: integer points in \(\mathbb R^8\) with even sum, plus half-integer
  points with even sum.
- Spherical cut: enumerate lattice points by squared norm, keep the `M`
  lowest-norm points, and normalize to unit average **symbol** energy.
- `lex` bit labels, or a nearest-neighbour-walk `gray` heuristic.
- Vector `modulate`/`demodulate`: a D4 symbol is two complex components, an
  E8 symbol is four.

The class is **not** a `Modulation` subclass.  The core `Link` carries one
complex value per symbol; a D4/E8 symbol is a vector.  The lattice appendix
therefore uses its own vector AWGN simulation in
`modcode.metrics.lattice`, not `run_link`.

### 1.1 Packing facts

For D4 the minimum norm is \(\sqrt2\), so the packing radius is
\(\sqrt2/2\); the fundamental cell volume is 2 and the center density is
\(1/8\).

For E8 the minimum norm is also \(\sqrt2\), so the packing radius is
\(\sqrt2/2\); the fundamental cell volume is 1 and the center density is
\(1/16\).  E8 is the densest lattice known in 8 dimensions (Viazovska 2016),
but that statement is about asymptotic packing, not about a particular
small spherical cut.

---

## 2. The dimension-matching correction

A D4-16 symbol carries 4 bits over 4 real dimensions, so it carries
\(4/4=1\) bit per real dimension.  16-QAM carries 4 bits over 2 real
dimensions, or \(4/2=2\) bits per real dimension.  Comparing D4-16 with
16-QAM at equal bits per **symbol** is therefore not a matched comparison: D4
uses twice as many real dimensions for the same bits.

The fair comparisons used here are:

| Comparison | \(\eta=b/K\) | Why |
|---|---:|---|
| D4-16 vs QPSK | 1 | both carry 1 bit per real dimension |
| D4-256 vs 16-QAM | 2 | both carry 2 bits per real dimension |
| E8-256 vs QPSK | 1 | both carry 1 bit per real dimension |

The simulation normalizes every constellation to unit average energy per real
dimension and uses

\[
\sigma^2=\frac{1}{2\eta\,(E_b/N_0)}
\]

per real dimension.  This is the only way the `d_min` comparison and the BER
simulation answer the same question.

---

## 3. Minimum distances

At unit energy per real dimension:

| Constellation | \(\eta\) | \(d_{\min}\) | Nominal gain vs 2D reference |
|---|---:|---:|---:|
| QPSK | 1 | 2.000 | — |
| D4-16 | 1 | 2.000 | 0 dB |
| E8-256 | 1 | 2.744 | 2.75 dB |
| 16-QAM | 2 | 0.894 | — |
| D4-256 | 2 | 1.069 | 1.55 dB |

The D4-16 result is the important correction: its **raw** \(d_{\min}\) looks
larger than 16-QAM only because the dimensions are mismatched.  At equal
\(\eta\) and equal energy per dimension it ties QPSK, and its BER is worse
because the labeling is not Gray-optimal and the constellation has more
nearest neighbours.

The D4-256 and E8-256 rows are genuine \(d_{\min}\) gains at matched
dimensions, but \(d_{\min}\) is an asymptotic proxy, not a BER result.

---

## 4. BER at matched dimensions

`examples/lattice_lab.py` simulates hard-decision BER with the same
per-dimension noise for every scheme.  Selected values:

| Eb/N0 | QPSK | D4-16 | E8-256 | 16-QAM | D4-256 |
|---:|---:|---:|---:|---:|---:|
| 4 dB | 1.26e-2 | 3.09e-2 | 1.71e-2 | 5.94e-2 | 1.44e-1 |
| 6 dB | 2.37e-3 | 6.78e-3 | 1.11e-3 | 2.79e-2 | 6.74e-2 |
| 8 dB | 1.40e-4 | 6.20e-4 | 0 | 9.35e-3 | 2.00e-2 |
| 12 dB | 0 | 0 | 0 | 1.10e-4 | 8.00e-5 |

What this shows:

- **D4-16 is not a win.**  It ties QPSK in \(d_{\min}\) but loses in BER at
  low/moderate SNR because of its labelling and neighbour structure.
- **E8-256 beats QPSK above roughly 5–6 dB**, but loses at 4 dB.  The measured
  gain at \(10^{-3}\) is modest, well below the nominal 2.75 dB \(d_{\min}\)
  advantage.
- **D4-256 beats 16-QAM only near 11–12 dB**, not across the waterfall.  The
  nominal 1.55 dB \(d_{\min}\) advantage is asymptotic and the finite
  constellation/labelling erodes it.

These are measured results for this construction, not literature claims.

---

## 5. PAPR diagnostic

`examples/lattice_lab.py` treats each of the \(K/2\) complex component streams
as a separate OFDM block of 64 components and reports the maximum PAPR across
streams.  This is a **diagnostic**, not a full system comparison: the streams
could be interleaved, mapped to different time slots, or assigned different
power in a real system.

Measured 99.9%-tile values:

| Modulation | 99.9% PAPR |
|---|---:|
| 256-QAM | 10.01 dB |
| D4-256 | 10.33 dB |
| E8-256 | 10.85 dB |

So this diagnostic shows a **0.3–0.9 dB** penalty, not the 3–7 dB figures that
appear when the number of effective subcarriers is multiplied without keeping
the comparison fixed.  A real PAPR claim needs a defined frame mapping and a
defined power normalization; this note does not make one.

---

## 6. Complexity

The implementation uses a naive nearest-neighbour search: `O(M K)` per
symbol.  Measured on this machine with 2,000-symbol batches:

| Modulation | Naive demodulation |
|---|---:|
| 256-QAM | ~1–2 × 10⁷ symbols/s |
| D4-256 | ~1.4 × 10⁵ symbols/s |
| E8-256 | ~8 × 10⁴ symbols/s |

The per-axis QAM slicer is vectorized and \(O(N)\), so the contrast is large
(roughly two orders of magnitude in this Python implementation).  Lattice-aware
algorithms (Conway–Sloane for D4, sphere decoding for E8) can reduce this to
\(O(K)\)-ish per symbol in hardware or optimized code, but those algorithms are
**not implemented here**.  The table is the honest software cost of this
appendix, and the exact throughput is machine-dependent.

---

## 7. What is claimed, measured, and not claimed

**Tested claims**

- D4-16 and QPSK have the same \(d_{\min}\) at matched \(\eta\) and energy per
  dimension.
- D4-256 has larger \(d_{\min}\) than 16-QAM at matched \(\eta\).
- E8-256 has larger \(d_{\min}\) than QPSK at matched \(\eta\).
- Round trips, unit average energy, and label uniqueness.

**Measured, with caveats**

- BER curves at matched \(\eta\) and per-dimension energy.
- A component-stream PAPR diagnostic.
- Naive software demodulation throughput.

**Not claimed**

- That lattice modulation replaces QAM in radio.
- The 2/3/3.5 dB literature figures quoted in the proposal.
- That D4-16 beats 16-QAM; that comparison is dimensionally unfair.
- A full OFDM PAPR result; the diagnostic is not a frame-level study.
- Link/`ber_sweep` integration, soft LLRs, or coded lattice links.
- Optical-channel gains; the proposal's optical extension is future work.

---

## 8. Reproduce

```bash
python examples/lattice_lab.py
python -m pytest tests/test_lattice.py -q
```

The script writes `examples/figures/lattice.png` and prints the distance,
BER, PAPR, and throughput tables.

---

## 9. References

- J. H. Conway and N. J. A. Sloane, *Sphere Packings, Lattices and Groups*,
  3rd edition, Springer, 1999.
- M. Viazovska, “The sphere packing problem in dimension 8,” *Annals of
  Mathematics* **185** (2017), 991–1015.
- G. D. Forney and G. Ungerboeck, “Modulation and coding for linear Gaussian
  channels,” *IEEE Transactions on Information Theory* **44** (1998),
  2384–2415.
