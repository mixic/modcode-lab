# Cyclic Trellis-Constrained Index Modulation (CT-IM)

### A constant-weight index code, its cyclic invariance, and the limits of the proposal

*Research note for [modcode-lab](../README.md)*

---

## Abstract

Classical OFDM index modulation (OFDM-IM) conveys bits by choosing which
subcarriers are active, but the index pattern itself has no error protection.
This note develops a concrete constant-weight version of the proposed CT-IM:

- the index pattern is one of the seven weight-3 codewords of Hamming(7,4);
- those seven words form a single cyclic orbit of size 7;
- their minimum Hamming distance is 4, so one index-bit error is correctable;
- the price is a smaller index alphabet: \(\log_2 7\approx2.807\) bits, not
  the 4 bits of the full Hamming(7,4) code, and a measurable rate loss against
  unconstrained OFDM-IM.

Two claims in the original proposal do not survive scrutiny and are corrected
here:

1. **Constant weight is not free.** Hamming(7,4) has codeword weights 0, 3, 4,
   and 7. Fixing the number of active subcarriers requires selecting the
   weight-3 orbit; it is a nonlinear constant-weight codebook, not a linear
   cyclic code, and it carries fewer than 4 index bits.
2. **Blind cyclic-shift detection from the syndrome is impossible.** Every
   cyclic shift of a codeword is again a codeword, so every shift has zero
   syndrome. The syndrome cannot identify the shift. Shift selection needs
   side information (signalling, a pilot, or a key), and that information is
   not free.

The note gives the algebraic theorems, the rate accounting, the implementation
map, and a literature-positioning checklist. It does **not** claim novelty or
publishability; those require a systematic literature search and the
BER/PAPR experiments listed at the end.

---

## 1. Classical OFDM-IM in one paragraph

Split \(N\) subcarriers into groups of \(n\). In each group, activate
\(k<n\) subcarriers and leave the rest silent. The active-set choice carries
\(\lfloor\log_2\binom{n}{k}\rfloor\) index bits; the \(k\) QAM symbols carry
\(k\log_2 M\) symbol bits. For \(n=7,k=3\), there are \(\binom{7}{3}=35\)
patterns, so a practical power-of-two table uses 32 patterns and 5 index bits.

The active-set choice is a combinatorial lookup. A single activity error moves
the received pattern to a different valid pattern; the index bits are then
wrong, and there is no redundancy in the index space.

---

## 2. The constant-weight index code

### 2.1 Hamming(7,4) weight spectrum

The binary \([7,4,3]\) Hamming code has weight distribution

| weight | number of words |
|---:|---:|
| 0 | 1 |
| 3 | 7 |
| 4 | 7 |
| 7 | 1 |

The weight-3 words are the object of interest:

```text
0001011
0010110
0101100
0110001
1000101
1011000
1100010
```

They all have weight 3, so the number of active subcarriers is fixed at
\(k=3\).

### 2.2 One cyclic orbit

**Theorem 1 (cyclic closure).** The set \(W\) of seven weight-3 Hamming
codewords is invariant under every cyclic shift, and the action of \(C_7\) on
\(W\) is transitive.

*Proof.* Hamming(7,4) is a cyclic code, so a cyclic shift of a codeword is a
codeword. A cyclic shift preserves Hamming weight, so it maps \(W\) into
itself. The shift is invertible, hence a permutation of \(W\).

It remains to see that \(W\) is one orbit. A word fixed by a nonzero cyclic
shift is constant on each cycle of that shift. Since 7 is prime, a nonzero
shift is a 7-cycle, so a fixed binary word is either all zero or all one;
neither has weight 3. Thus no nonzero shift fixes a word in \(W\). The orbit
sizes therefore divide 7 and are not 1, so every orbit has size 7. Since
\(|W|=7\), there is exactly one orbit. \(\square\)

The set is not a linear subspace: it does not contain the zero word, and the
sum of two distinct weight-3 words is a weight-4 word. It is a
**constant-weight cyclic codebook** (a cyclic orbit), not a linear cyclic code
in the sense of an ideal in \(\mathbb{F}_2[C_7]\).

### 2.3 Minimum distance and correction

**Theorem 2 (distance).** Distinct words in \(W\) have Hamming distance 4.

*Proof.* The difference of two distinct binary words is a nonzero Hamming
codeword. The nonzero Hamming weights are 3, 4, and 7. The sum of two words of
odd weight 3 has even weight; the only even possibility is 4. Hence the
distance is 4. \(\square\)

A code with minimum distance \(d=4\) corrects \(t=\lfloor(d-1)/2\rfloor=1\)
error. The decoder is a nearest-word search over the seven patterns; it is
exact for one index-bit error and may tie for larger errors.

### 2.4 Index rate

There are only seven patterns. A fixed-length binary label needs 3 bits for
eight labels, so one label is reserved:

| quantity | value |
|---|---:|
| index patterns | 7 |
| index entropy | \(\log_2 7\approx2.807\) bits |
| fixed-length label | 3 bits (label 7 reserved) |
| active subcarriers | 3 |
| QPSK symbol bits | 6 |
| bits per 7-subcarrier block | 9 |
| bits per subcarrier | \(9/7\approx1.286\) |

The effective index entropy is 2.807 bits, not 3; the 3-bit label is a
convenient fixed-length interface with one forbidden word. A 4-bit index, as
in the full Hamming code, would require variable weight (3 or 4) and therefore
a variable symbol budget.

---

## 3. The CT-IM mapper

For QPSK (\(M=4\)):

```text
input: 9 bits = 3 index bits + 6 symbol bits
index label in 0..6 -> one weight-3 pattern
3 QPSK symbols -> the three active positions
silent positions -> 0
```

The decoder:

```text
hard activity: a_j = 1 if |y_j| > threshold
nearest weight-3 pattern -> index label
demodulate the three active positions -> symbol bits
```

The implementation is `modcode.modulation.CTIM`; the codebook is
`ConstantWeightCyclicCode`. The threshold detector is deliberately simple and
is not a channel-aware activity detector.

---

## 4. Rate comparison

Classical OFDM-IM with \(n=7,k=3\):

- 35 patterns, 5 index bits (32 used);
- 3 QPSK symbols, 6 symbol bits;
- 11 bits per block, \(11/7\approx1.571\) bits/subcarrier;
- no index error correction.

CT-IM:

- 7 patterns, \(\log_2 7\approx2.807\) bits of index entropy;
- 3 QPSK symbols, 6 symbol bits;
- 9 bits per block, \(9/7\approx1.286\) bits/subcarrier;
- one index-bit error corrected.

The rate cost is about \(0.285\) bits/subcarrier (18%) at QPSK. Whether the
index correction is worth that cost is exactly the experiment: under a channel
where index errors dominate, it can win; under a clean high-SNR channel, the
unconstrained alphabet wins on rate. The note does not assert a crossover
without BER data.

---

## 5. The cyclic-shift problem

### 5.1 Every shift is a codeword

**Theorem 3.** For every \(i\in W\) and every \(s\in\mathbb{Z}/7\mathbb{Z}\),
the shifted word \(\sigma_s(i)\) is also in \(W\), and therefore
\(H\sigma_s(i)^T=0\).

This is immediate from Theorem 1 and the definition of the parity check.

### 5.2 The syndrome cannot identify the shift

**Corollary.** A zero syndrome is consistent with all seven shifts. The
syndrome therefore cannot tell the receiver which shift the transmitter chose.

The proposal's “blind detection by zero syndrome” is not valid: the correct
shift and all six wrong shifts have the same syndrome. Adding noise does not
help uniquely, because the nearest-codeword set still contains the whole
orbit.

### 5.3 What shift selection actually costs

If the transmitter chooses a shift to minimize PAPR, the receiver must be told
the shift (or must know it as a key). The side information costs
\(\lceil\log_2 7\rceil=3\) bits per block. If those bits come from the same
payload budget, the CT-IM payload falls from 9 to 6 bits per block
(\(6/7\approx0.857\) bits/subcarrier), before any symbol-rate accounting.
A smaller shift set reduces the side information but also reduces the PAPR
gain.

Alternatives:

- **Pilot/control channel:** send the shift explicitly; clean but consumes
  resources.
- **Keyed shift:** use the shift as a covert-channel key; not a standard
  communication scheme, and it must be synchronized.
- **Orbit data:** make the data label the orbit rather than the word. Since
  the seven words form one orbit, this collapses the index alphabet to one
  orbit and destroys index modulation. This is the mathematical reason the
  “free shift” idea cannot work.

PAPR reduction by shift selection is therefore an **upper bound with perfect
side information**, not a free gain. The current implementation supports
`best_shift` and `encode(..., return_shift=True)` for measuring that upper
bound; it does not claim blind detection.

---

## 6. What is actually novel, and what is not

The implemented object is:

- a constant-weight cyclic codebook of size 7 and distance 4;
- an OFDM-IM mapper using it for the index pattern;
- a nearest-word index decoder;
- optional PAPR-oriented shift selection with explicit side information.

The **tail-biting trellis** in the proposal's name is not implemented. The
codebook is cyclic; a minimal trellis realization could be built, but the
current code is a codebook lookup, not a convolutional trellis. The phrase
“trellis-constrained” is retained as the proposal's name, not as a claim that
a trellis is running.

The literature positioning in the proposal is a hypothesis, not a result. A
proper search should cover:

- OFDM-IM and generalized/index-modulation variants;
- trellis-coded and memory-based index selection (the EGIM family);
- cyclic codes and constant-weight codes used as index alphabets;
- tail-biting trellises and their use in modulation;
- PAPR reduction by cyclic shift or phase selection in OFDM.

The claim that “index modulation + cyclic code” or “index modulation +
tail-biting” is essentially absent must be verified against IEEE, arXiv, and
patent literature before any novelty claim. This note makes no such claim.

---

## 7. Limitations of the current implementation

- Only \(n=7\), weight 3, and the Hamming(7,4) orbit are implemented.
- The index code is nonlinear; the decoder is a nearest-word search.
- The activity detector is a fixed magnitude threshold; no channel estimation,
  fading, or noise-adaptive threshold.
- Symbol bits are unprotected; a correct index pattern does not guarantee
  correct QAM symbols.
- No BER, PAPR, or throughput baseline is implemented yet.
- No trellis, no BCH/RS extension, and no Rayleigh channel.
- The shift-selection path is an upper-bound tool; it requires side
  information and is not a blind receiver.

---

## 8. Next experiments

1. **BER, AWGN:** CT-IM versus classical OFDM-IM with a matched activity
   detector; separate index-bit and symbol-bit errors.
2. **Rayleigh:** the case where index correction is supposed to matter.
3. **PAPR CCDF:** no shift, best shift with perfect side information, and
   best shift with a realistic side-information budget.
4. **Constant-weight families:** the \([7,3,4]\) simplex code (weight 4),
   shortened BCH codes, and other cyclic orbits.
5. **Activity detection:** energy threshold versus ML over the codebook.
6. **Trellis realization:** a minimal trellis for the seven-word orbit, if it
   buys decoding complexity over lookup.
7. **Literature audit:** the novelty claims in §6.

---

## 9. Reproducibility

```bash
python -m pytest tests/test_ct_im.py
python examples/ct_im_lab.py
```

The demo writes `examples/figures/ct_im_codebook.png`: the seven patterns, the
pairwise distance matrix, and the rate comparison. It also verifies all 49
single index-bit errors.

Files:

- `src/modcode/modulation/ct_im.py`
- `src/modcode/viz/ct_im.py`
- `examples/ct_im_lab.py`
- `tests/test_ct_im.py`

---

## References

1. E. Basar, U. Aygolu, E. Panayirci, and H. V. Poor, “Orthogonal frequency
   division multiplexing with index modulation,” *IEEE Transactions on Signal
   Processing* **61** (2013), 5536–5549.
2. E. Basar, “Index modulation techniques for 5G wireless networks,” *IEEE
   Communications Magazine* **54** (2016), 168–175.
3. F. J. MacWilliams and N. J. A. Sloane, *The Theory of Error-Correcting
   Codes*, North-Holland, 1977.
4. S. Lin and D. J. Costello, *Error Control Coding*, 2nd edition, Pearson,
   2004.
5. The `TailBitingTCM8PSK` example in this repository, for the tail-biting
   trellis background.
