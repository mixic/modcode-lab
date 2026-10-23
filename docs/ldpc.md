# LDPC(1008, 506): a short quasi-cyclic code with min-sum decoding

*Core FEC extension. This is a pedagogical finite-length LDPC code, not a 5G
NR LDPC implementation.*

The convolutional code in the core is a low-latency, low-complexity baseline.
The natural next comparison is an iteratively decoded code.  This note adds a
short regular (3,6) LDPC code and a normalized min-sum decoder, then compares
it with the convolutional code in the existing `Link` pipeline.  The other
family used in 5G, polar coding with successive-cancellation list decoding, is
documented in [`polar.md`](polar.md).

---

## 1. Construction

The parity-check matrix is quasi-cyclic:

- block size $z=84$;
- 6 block rows and 12 block columns;
- $N=1008$, $M=504$;
- every check has degree 6, every variable has degree 3.

Each nonzero block is a circulant permutation with a deterministic pseudo-random
shift.  Two of the 504 checks are linearly dependent, so the rank is 502 and
the code carries $k=506$ information bits:

$$
R=\frac{506}{1008}\approx0.502.
$$

The generator is obtained by GF(2) row reduction of $H$.  The code is
systematic in the **free columns** of the reduced matrix, not necessarily in
the first $k$ positions.  `encode` maps $k$ data bits to the corresponding
free positions; `decode` reads those positions back.

The parity-check matrix and generator satisfy

$$
G H^T = 0 \pmod 2,
$$

which is tested directly.

### 1.1 Scale against the standards

| Code | $N$ | $k$ | Reference |
|---|---:|---:|---|
| this repository | 1008 | 506 | — |
| 5G NR LDPC, base graph 1 | up to 25344 | up to 8448 | 3GPP TS 38.212 |
| DVB-S2X normal frame | 64800 | up to 58320 | ETSI EN 302 307-2 |

The code is $8.4$ times shorter than a 5G NR information block and $64$ times
shorter than a DVB-S2X frame, and its degree distribution is a regular $(3,6)$
construction, not a density-evolved one.  Nothing measured at $N=1008$ transfers
to those lengths: the waterfall slope, the error floor, and the iteration count
all depend on $N$.  The extrapolation rule and the block-length evidence for it
are in [`scope-and-limitations.md`](scope-and-limitations.md) §2.

---

## 2. Decoder

`LDPC36.decode_soft` runs normalized min-sum belief propagation:

- variable-to-check messages start at the channel LLR;
- each check sends the minimum magnitude excluding the destination edge, with
  the sign of the other edges and a normalization factor `alpha=0.8`;
- each variable sums the channel LLR and all incoming check messages;
- hard decisions are tested against the syndrome after every iteration;
- decoding stops early on a valid codeword, or after 50 iterations.

This is a practical approximation to sum-product, not exact log-domain BP.  It
uses no layered schedule, offset correction, or early-stopping heuristic beyond
the syndrome check.

`decode` converts hard bits to $\pm4$ LLRs and runs the same decoder.

---

## 3. Results

`examples/ldpc_lab.py` compares:

- convolutional K=3, rate 0.485, soft Viterbi;
- LDPC(1008,506), rate 0.502, min-sum;
- uncoded BPSK theory.

Measured at 100,000 bits per point with fixed seeds:

| Eb/N0 | convolutional | LDPC |
|---:|---:|---:|
| 0 dB | 9.20e-2 | 1.58e-1 |
| 1 dB | 4.67e-2 | 1.02e-1 |
| 2 dB | 1.65e-2 | 1.35e-3 |
| 3 dB | 3.91e-3 | 0 |
| 4 dB | 5.00e-4 | 0 |
| 5 dB | 1.00e-4 | 0 |

The finite-length LDPC curve crosses the convolutional curve between 1 and
2 dB and then falls steeply.  The measured gain at $10^{-3}$ is roughly
1.5–2 dB for this pair, not the 3–4 dB quoted for long standardized LDPC codes.
That is expected: this code is short, the decoder is min-sum, and no girth
optimization or density evolution was performed.  The crossing point is a
property of this pair at this length, not a ranking of the two families; §1.1
and [`scope-and-limitations.md`](scope-and-limitations.md) §2 give the reason.

The zeros at 3–5 dB are zero-error runs at this sample size, not an error-floor
claim.  A proper floor study needs much larger blocks and many more bits.

The companion note [`polar.md`](polar.md) adds the other 5G family on the same
`Link` pipeline; the two-family BER and decoding-cost comparison is Figure 24 in
[`figures.md`](figures.md).

---

## 4. Complexity and latency

The decoder cost is $O(I \cdot E)$, where $I$ is the iteration count and
$E=3024$ is the number of edges.  The block latency is $N=1008$ bits; the
code is block-oriented, unlike the streaming convolutional code.

No operation-count benchmark against 5G LDPC is claimed.  The comparison here
is a reproducible Python implementation of two decoders, not a hardware study.

---

## 5. What is claimed, measured, and not claimed

**Tested**

- H is (3,6) regular.
- rank 502, k 506, rate 506/1008.
- $G H^T=0$.
- noiseless hard and soft round trips.
- BP link decoding at 3 dB with BER below $10^{-3}$.

**Measured**

- finite-length BER comparison with the K=3 convolutional code.
- the crossover and the steep waterfall.

**Not claimed**

- a 3–4 dB gain over convolutional coding in general;
- the crossover point as evidence about the two families at standards-scale
  block lengths;
- 5G/Wi-Fi LDPC performance;
- an error floor;
- optimal girth, degree distribution, or decoding schedule;
- a complexity or latency advantage over convolutional codes.

---

## 6. Reproduce

```bash
python examples/ldpc_lab.py
python -m pytest tests/test_ldpc.py -q
```

Figure: `examples/figures/ldpc.png`.

---

## 7. References

- R. G. Gallager, *Low-Density Parity-Check Codes*, MIT Press, 1963.
- D. J. C. MacKay, “Good error-correcting codes based on very sparse
  matrices,” *IEEE Transactions on Information Theory* **45** (1999),
  399–431.
- T. J. Richardson and R. L. Urbanke, *Modern Coding Theory*, Cambridge
  University Press, 2008.
