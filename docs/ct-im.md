# Cyclic Trellis-Constrained Index Modulation (CT-IM)

### A constant-weight index code, its cyclic invariance, and the limits of the proposal

*Research note for [modcode-lab](../README.md)*

> **Scope:** this is an exploratory appendix, not part of the core
> modulation/FEC simulator and not a paper-level result. It is kept because
> the codebook and its limits are tested; the open milestones are marked in
> §10. See [`exploratory-labs.md`](exploratory-labs.md) for the project-level
> separation between the core and the appendices, and
> [`scope-and-limitations.md`](scope-and-limitations.md) for the modelling
> limits.

---

## Abstract

Classical OFDM index modulation (OFDM-IM) conveys bits by choosing which
subcarriers are active, but the index pattern itself has no error protection.
This note develops a concrete constant-weight version of the proposed CT-IM:

- the index pattern is one of the seven weight-3 codewords of Hamming(7,4);
- those seven words form a single cyclic orbit of size 7;
- their minimum Hamming distance is 4, so one index-bit error is correctable;
- the price is a smaller index alphabet: $\log_2 7\approx2.807$ bits, not
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
3. **A cyclic shift cannot reduce PAPR.** A shift of the whole frequency-domain
   block is a phase ramp in time, so the PAPR is unchanged to machine
   precision (§5.3). The PAPR argument for shift selection is void, not merely
   expensive; the measured deviation is below $10^{-14}$ dB.

The experiments answer the index-protection question, not the per-bandwidth
question:

- **The index code does what it is designed to do.**  Its index-bit BER is below
  the unprotected OFDM-IM alphabet at every measured point, in AWGN and in
  known-CSI Rayleigh: at 4 dB, $4.2\times10^{-4}$ against $1.5\times10^{-2}$
  in AWGN, and $2.5\times10^{-2}$ against $1.3\times10^{-1}$ in Rayleigh.
- **But the comparison is not rate-matched.**  CT-IM carries $1.286$
  bit/subcarrier and classical OFDM-IM $1.571$, so CT-IM also has the lower
  total BER at equal Eb/N0.  Declaring a winner per unit bandwidth would need
  CT-IM at a higher modulation order or OFDM-IM at a lower one; that experiment
  is not run.  What the data support is narrower: the constant-weight index
  layer buys index-bit protection, and the price is rate, not error rate.


The note gives the algebraic theorems, the rate accounting, the implementation
map, and a literature-positioning checklist. It does **not** claim novelty or
publishability; those require a systematic literature search and the
BER/PAPR experiments listed at the end.

---

## 1. Classical OFDM-IM in one paragraph

Split $N$ subcarriers into groups of $n$. In each group, activate
$k<n$ subcarriers and leave the rest silent. The active-set choice carries
$\lfloor\log_2\binom{n}{k}\rfloor$ index bits; the $k$ QAM symbols carry
$k\log_2 M$ symbol bits. For $n=7,k=3$, there are $\binom{7}{3}=35$
patterns, so a practical power-of-two table uses 32 patterns and 5 index bits.

The active-set choice is a combinatorial lookup. A single activity error moves
the received pattern to a different valid pattern; the index bits are then
wrong, and there is no redundancy in the index space.

---

## 2. The constant-weight index code

### 2.1 Hamming(7,4) weight spectrum

The binary $[7,4,3]$ Hamming code has weight distribution

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
$k=3$.

### 2.2 One cyclic orbit

**Theorem 1 (cyclic closure).** The set $W$ of seven weight-3 Hamming
codewords is invariant under every cyclic shift, and the action of $C_7$ on
$W$ is transitive.

*Proof.* Hamming(7,4) is a cyclic code, so a cyclic shift of a codeword is a
codeword. A cyclic shift preserves Hamming weight, so it maps $W$ into
itself. The shift is invertible, hence a permutation of $W$.

It remains to see that $W$ is one orbit. A word fixed by a nonzero cyclic
shift is constant on each cycle of that shift. Since 7 is prime, a nonzero
shift is a 7-cycle, so a fixed binary word is either all zero or all one;
neither has weight 3. Thus no nonzero shift fixes a word in $W$. The orbit
sizes therefore divide 7 and are not 1, so every orbit has size 7. Since
$|W|=7$, there is exactly one orbit. $\square$

The set is not a linear subspace: it does not contain the zero word, and the
sum of two distinct weight-3 words is a weight-4 word. It is a
**constant-weight cyclic codebook** (a cyclic orbit), not a linear cyclic code
in the sense of an ideal in $\mathbb{F}_2[C_7]$.

### 2.3 Minimum distance and correction

**Theorem 2 (distance).** Distinct words in $W$ have Hamming distance 4.

*Proof.* The difference of two distinct binary words is a nonzero Hamming
codeword. The nonzero Hamming weights are 3, 4, and 7. The sum of two words of
odd weight 3 has even weight; the only even possibility is 4. Hence the
distance is 4. $\square$

A code with minimum distance $d=4$ corrects $t=\lfloor(d-1)/2\rfloor=1$
error. The decoder is a nearest-word search over the seven patterns; it is
exact for one index-bit error and may tie for larger errors.

### 2.4 Index rate

There are only seven patterns. A fixed-length binary label needs 3 bits for
eight labels, so one label is reserved:

| quantity | value |
|---|---:|
| index patterns | 7 |
| index entropy | $\log_2 7\approx2.807$ bits |
| fixed-length label | 3 bits (label 7 reserved) |
| active subcarriers | 3 |
| QPSK symbol bits | 6 |
| bits per 7-subcarrier block | 9 |
| bits per subcarrier | $9/7\approx1.286$ |

The effective index entropy is 2.807 bits, not 3; the 3-bit label is a
convenient fixed-length interface with one forbidden word. A 4-bit index, as
in the full Hamming code, would require variable weight (3 or 4) and therefore
a variable symbol budget.

---

## 3. The CT-IM mapper

For QPSK ($M=4$):

```text
input: 9 bits = 3 index bits + 6 symbol bits
index label in 0..6 -> one weight-3 pattern
3 QPSK symbols -> the three active positions
silent positions -> 0
```

The decoder has two index detectors:

```text
threshold: a_j = 1 if |y_j| > threshold, then nearest weight-3 pattern
ML:        argmin over patterns of the squared distance to the nearest
           active symbol on the pattern plus the energy elsewhere
```

The ML detector is exact for equal-energy symbols and known CSI, and the noise
variance cancels because it is common to all patterns.  It is the matched
baseline used in the experiments; the threshold detector is kept to show what a
non-channel-aware activity detector costs (about 6 dB at $10^{-2}$ index BER in
the AWGN sweep).

The implementation is `modcode.modulation.CTIM`; the codebook is
`ConstantWeightCyclicCode`.  `CTIM(M=4, detector="ml")` selects the ML path.

---

## 4. Rate comparison

Classical OFDM-IM with $n=7,k=3$:

- 35 patterns, 5 index bits (32 used);
- 3 QPSK symbols, 6 symbol bits;
- 11 bits per block, $11/7\approx1.571$ bits/subcarrier;
- no index error correction.

CT-IM:

- 7 patterns, $\log_2 7\approx2.807$ bits of index entropy;
- 3 QPSK symbols, 6 symbol bits;
- 9 bits per block, $9/7\approx1.286$ bits/subcarrier;
- one index-bit error corrected.

The rate cost is about $0.285$ bits/subcarrier (18%) at QPSK. Whether the
index correction is worth that cost is exactly the experiment: under a channel
where index errors dominate, it can win; under a clean high-SNR channel, the
unconstrained alphabet wins on rate.  §8 reports the measured answer.

Two further constant-weight families are implemented for the experiments:

| family | $n$ | active | words | $d$ | index bits | symbol bits | bits/block | bits/subcarrier |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Hamming(7,4) orbit | 7 | 3 | 7 | 4 | 3 | 6 | 9 | 1.286 |
| shortened orbit | 6 | 3 | 4 | 4 | 2 | 6 | 8 | 1.333 |
| BCH(15,7) weight-7 orbit | 15 | 7 | 15 | 8 | 4 | 14 | 18 | 1.200 |
| classical OFDM-IM | 7 | 3 | 32 | 2 | 5 | 6 | 11 | 1.571 |

The shortened family is CT-IM-1.  Puncturing the weight-3 orbit is not
constant-weight (deleting an active coordinate gives weight 2, deleting a silent
one gives weight 3); *shortening* — fixing coordinate 0 to zero and deleting
it — keeps weight 3 and distance 4.  The price is a four-word alphabet and two
index bits, and the budget is explicit: three active subcarriers in six.

The BCH family is CT-IM-6.  The linear BCH(15,7) code is not constant-weight;
its weight-7 codewords form one cyclic orbit of 15 words with pairwise distance
8, and that orbit is the index alphabet.  It is the most protected and the
slowest of the three index codes.

---

## 5. The cyclic-shift problem

### 5.1 Every shift is a codeword

**Theorem 3.** For every $i\in W$ and every $s\in\mathbb{Z}/7\mathbb{Z}$,
the shifted word $\sigma_s(i)$ is also in $W$, and therefore
$H\sigma_s(i)^T=0$.

This is immediate from Theorem 1 and the definition of the parity check.

### 5.2 The syndrome cannot identify the shift

**Corollary.** A zero syndrome is consistent with all seven shifts. The
syndrome therefore cannot tell the receiver which shift the transmitter chose.

The proposal's “blind detection by zero syndrome” is not valid: the correct
shift and all six wrong shifts have the same syndrome. Adding noise does not
help uniquely, because the nearest-codeword set still contains the whole
orbit.

### 5.3 What shift selection actually costs

If the transmitter chooses a shift, the receiver must be told the shift (or
must know it as a key). The side information costs
$\lceil\log_2 7\rceil=3$ bits per block. If those bits come from the same
payload budget, the CT-IM payload falls from 9 to 6 bits per block
($6/7\approx0.857$ bits/subcarrier), before any symbol-rate accounting.

There is a stronger problem: **the shift does not reduce PAPR at all.**  A
cyclic shift of the whole frequency-domain block multiplies the time-domain
waveform by $e^{j2\pi s n/N}$, so $|x[n]|$ is unchanged.  The experiment in
`examples/ct_im_experiments.py` computes the PAPR of the 7-subcarrier block
placed in a 64-subcarrier band at offsets $0,\dots,6$; the largest deviation
across 500 blocks is $5.3\times10^{-15}$ dB.  The "best shift" is therefore any
shift, and the three bits of side information buy nothing.

The implemented `best_shift` and `encode(..., return_shift=True)` are retained
so the invariance can be demonstrated and tested; they are not a PAPR tool.
A pattern-only rotation (rotating the active set without rotating the data
symbols) would change the waveform, but it also changes the symbol-to-subcarrier
mapping and therefore the scheme; it is not a free optimisation of the same
transmission.

Alternatives:

- **Pilot/control channel:** send the shift explicitly; clean but consumes
  resources.
- **Keyed shift:** use the shift as a covert-channel key; not a standard
  communication scheme, and it must be synchronized.
- **Orbit data:** make the data label the orbit rather than the word. Since
  the seven words form one orbit, this collapses the index alphabet to one
  orbit and destroys index modulation. This is the mathematical reason the
  “free shift” idea cannot work.

PAPR reduction by cyclic shift is therefore not an upper bound with perfect side
information — it is a null operation.  The implementation supports
`best_shift` and `encode(..., return_shift=True)` for measuring the invariance;
it does not claim a PAPR gain or blind detection.

---

## 6. What is actually novel, and what is not

The implemented object is:

- a constant-weight cyclic codebook of size 7 and distance 4;
- an OFDM-IM mapper using it for the index pattern;
- a nearest-word decoder and an exact ML index detector;
- a shortened constant-weight variant, a BCH(15,7) weight-7 orbit, and a
  classical OFDM-IM baseline, all with tests and BER curves;
- a measured PAPR result showing that shift selection is a null operation.

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
- PAPR reduction by pattern selection, symbol permutation, or phase rotation in
  OFDM (cyclic shift is ruled out by §5.3, but the literature on other levers
  still needs checking).

The claim that “index modulation + cyclic code” or “index modulation +
tail-biting” is essentially absent must be verified against IEEE, arXiv, and
patent literature before any novelty claim. This note makes no such claim.

---

## 7. Limitations of the current implementation

- The index codes are small by construction: length 6, 7, and 15 with two to
  four index bits.  Nothing here scales to large groups.
- The index code is nonlinear; the decoder is a nearest-word search.  The ML
  detector is exact only for equal-energy symbols and known CSI.
- Symbol bits are unprotected; a correct index pattern does not guarantee
  correct QAM symbols.
- The channel is AWGN or flat known-CSI Rayleigh per subcarrier.  There is no
  estimation error, correlation across subcarriers, or frequency selectivity
  inside a group; the multipath/Doppler path is a separate core note.
- The Rayleigh simulation gives every subcarrier an independent gain, which is
  optimistic for a real group of adjacent subcarriers.
- The PAPR is computed on the block waveform and, in the shift experiment, on a
  64-subcarrier placement; it is not a full standard-compliant OFDM symbol.
- No trellis realisation, no tail-biting, no soft index decoding, and no
  iterative index/symbol decoding.
- The shift-selection path is demonstrated to be PAPR-invariant; it is not a
  blind receiver and not a PAPR reduction tool.
- The literature audit (CT-IM-6 in the original list, now §10.4) is not done:
  no novelty claim is made.

---

## 8. Experiments and results

`examples/ct_im_experiments.py` runs every scheme through
`metrics.index_modulation.simulate_index_modulation`, which fixes the channel,
the Eb/N0 definition, and the detector.  Eb/N0 is per information bit and the
block energy is three unit-energy QPSK symbols, so the noise variance is
$3/(2\cdot\text{bits per block}\cdot\mathrm{Eb}/N_0)$ per real component.

**AWGN, 800 blocks per point.**  Index-bit / symbol-bit BER:

| Eb/N0 | CT-IM threshold | CT-IM ML | shortened ML | BCH orbit ML | OFDM-IM ML |
|---:|---:|---:|---:|---:|---:|
| 0 dB | 2.55e-1 / 2.08e-1 | 7.79e-2 / 8.87e-2 | 7.44e-2 / 8.56e-2 | 6.34e-2 / 9.37e-2 | 1.92e-1 / 1.26e-1 |
| 2 dB | 1.29e-1 / 1.07e-1 | 2.00e-2 / 2.90e-2 | 1.75e-2 / 2.81e-2 | 3.13e-3 / 2.34e-2 | 7.12e-2 / 4.65e-2 |
| 4 dB | 4.54e-2 / 3.79e-2 | 4.17e-4 / 4.38e-3 | 1.87e-3 / 7.50e-3 | 0 / 5.62e-3 | 1.50e-2 / 1.02e-2 |
| 6 dB | 3.75e-3 / 2.92e-3 | 0 / 4.17e-4 | 0 / 1.25e-3 | 0 / 8.93e-4 | 1.00e-3 / 0 |
| 8 dB | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 |

Three observations:

- **The ML detector is the matched baseline.**  It gains about 6 dB over the
  naive threshold detector at $10^{-2}$ index BER, and symbol-bit errors, not
  index-bit errors, dominate CT-IM above 2 dB.
- **The comparison is not rate-matched.**  CT-IM carries $1.286$
  bit/subcarrier and OFDM-IM $1.571$, so CT-IM's lower total BER is partly the
  lower rate.  Per index bit, where the comparison is like-for-like in
  protection, CT-IM is still better: $4.2\times10^{-4}$ against
  $1.5\times10^{-2}$ at 4 dB.  A per-bandwidth comparison needs a rate-matched
  pair and is not run.
- **The shortened family (CT-IM-1) trades alphabet size for rate.**  It carries
  $8/6=1.333$ bit/subcarrier against $9/7=1.286$ for the orbit with the same
  distance 4, but it has four words instead of seven; its measured BER is close
  to the orbit's.  The rate gain is real and small.

**Known-CSI Rayleigh, index-bit BER.**  Each subcarrier gets an independent
$\mathcal{CN}(0,1)$ gain, known at the receiver:

| Eb/N0 | CT-IM AWGN | CT-IM Rayleigh | OFDM-IM AWGN | OFDM-IM Rayleigh |
|---:|---:|---:|---:|---:|
| 0 dB | 7.79e-2 | 1.33e-1 | 1.92e-1 | 2.88e-1 |
| 2 dB | 2.00e-2 | 6.42e-2 | 7.12e-2 | 1.98e-1 |
| 4 dB | 4.17e-4 | 2.50e-2 | 1.50e-2 | 1.33e-1 |
| 6 dB | 0 | 7.50e-3 | 1.00e-3 | 7.70e-2 |
| 8 dB | 0 | 3.75e-3 | 0 | 4.05e-2 |
| 10 dB | 0 | 0 | 0 | 1.92e-2 |

CT-IM keeps the lower index-bit BER in Rayleigh as well, but the gap is not a
free win: the rate is lower, and the Rayleigh model is optimistic (independent
gains per subcarrier, known at the receiver).  The multipath note shows that
estimation error can cost more than that.  The honest reading is that the
constant-weight index layer protects index bits, and that its rate cost is the
price — a per-bandwidth verdict requires the rate-matched experiment.

**PAPR.**  Cyclic shift of the whole block is PAPR-invariant (measured maximum
deviation $5.3\times10^{-15}$ dB), so the three side-information bits buy
nothing.  The CT-IM and OFDM-IM CCDFs are close: mean PAPR 4.45 dB and 4.37 dB
respectively on the same 500 blocks, with $P(\mathrm{PAPR}>5\ \mathrm{dB})$
below the sample resolution.  The remaining PAPR lever is the active pattern
itself (CT-IM fixes three active subcarriers, OFDM-IM does too), not a shift.

---

## 9. Reproducibility

```bash
python -m pytest tests/test_ct_im.py tests/test_ct_im_experiments.py
python examples/ct_im_lab.py            # codebook figure
python examples/ct_im_experiments.py    # BER and PAPR figures
```

Figures:

- `examples/figures/ct_im_codebook.png` — patterns, distance matrix, rate
  comparison (the algebraic demo);
- `examples/figures/ct_im_ber.png` — AWGN and Rayleigh BER, index versus symbol
  bits;
- `examples/figures/ct_im_papr.png` — PAPR CCDF for CT-IM and OFDM-IM.

Files:

- `src/modcode/modulation/ct_im.py`
- `src/modcode/metrics/index_modulation.py`
- `src/modcode/viz/ct_im.py`
- `examples/ct_im_lab.py`
- `examples/ct_im_experiments.py`
- `tests/test_ct_im.py`
- `tests/test_ct_im_experiments.py`

---

## 10. Milestones

The project-level placement is in the
[README exploratory track](../README.md#exploratory-track).  This section records
the status and what each checkbox rests on.  Status: ✅ complete, ◐ partial,
⬜ open.

| Milestone | Status | Deliverable | Evidence |
|---|---|---|---|
| **CT-IM-0** | ✅ | Hamming(7,4) weight-3 index layer, mapper, decoder | `CTIM`, `ConstantWeightCyclicCode`, `tests/test_ct_im.py` |
| **CT-IM-1** | ✅ | Constant-weight variant with an explicit symbol budget | `ShortenedConstantWeightCode`: length 6, weight 3, 4 words, $d=4$, two index bits |
| **CT-IM-2** | ◐ | Baselines: classical OFDM-IM and EGIM | `ClassicalOFDMIM` (32 patterns, 5 index bits, 6 symbol bits); EGIM is **not** implemented because no publication-level definition is pinned in this note |
| **CT-IM-3** | ✅ | BER sweep and PAPR CCDF | `examples/ct_im_experiments.py`, `ct_im_ber.png`, `ct_im_papr.png`, §8 |
| **CT-IM-4** | ✅ | Rayleigh extension | `simulate_index_modulation(..., rayleigh=True)`, known-CSI table in §8 |
| **CT-IM-5** | ✅ | Shift-detection impossibility | §5, plus the PAPR-invariance measurement in §5.3 |
| **CT-IM-6** | ✅ | BCH(15,7) index layer | `BCHIndexCode`: weight-7 orbit, 15 words, $d=8$, four index bits |
| **CT-IM-7** | ✅ | Write-up with measured results | this note, §8 |

### 10.1 Status notes

- **CT-IM-1 is a shortening, not a puncturing.**  Puncturing the weight-3 orbit
  destroys constant weight (weights 2 and 3 appear).  The completed variant
  fixes coordinate 0 to zero and deletes it, keeping weight 3 and distance 4
  over four words.  The budget is explicit: three active subcarriers in six,
  two index bits, six QPSK symbol bits, eight bits per block.  The rate is
  $8/6=1.333$ bit/subcarrier against $9/7=1.286$ for the orbit; §8 shows the
  rate gain does not translate into a BER advantage at these rates.
- **CT-IM-2 is partial by design.**  The classical OFDM-IM baseline is fully
  specified and tested.  "EGIM" is a proposal label with no pinned reference in
  this note; inventing a construction and calling it EGIM would be a claim, not
  a baseline.  If a citation is added, the variant must use the same $n=7$,
  $k=3$, QPSK, channel, and detector as the other schemes.
- **CT-IM-3 answered a different question than planned.**  The PAPR experiment
  showed that cyclic shift is PAPR-invariant, so the three planned cases
  coincide.  The measured result is therefore a null result, recorded in §5.3
  and tested in `tests/test_ct_im_experiments.py`.
- **CT-IM-4 is known-CSI only.**  Each subcarrier gets an independent Rayleigh
  gain known at the receiver; there is no estimation error and no correlation
  across the group's subcarriers.  That is optimistic, and it is stated as
  such.
- **CT-IM-5 rules out a detector.**  It does not produce a blind receiver.
  §5.3 also rules out the PAPR rationale for the shift.
- **CT-IM-6 is an orbit, not the linear code.**  The full BCH(15,7) code is not
  constant-weight; the weight-7 orbit is constant-weight with $d=8$, and that
  orbit is the index alphabet.
- **CT-IM-7 is now a measured report.**  §8 replaces the earlier "future work"
  list with data.  The remaining open item is the literature audit (§6), which
  is outside the milestone table.

### 10.2 What is still open

- **EGIM baseline:** needs a pinned definition before coding.
- **Literature audit:** the novelty questions in §6 require IEEE/arXiv/patent
  search; no novelty claim is made.
- **Trellis realisation:** a minimal trellis for an orbit could replace the
  nearest-word lookup; not attempted.
- **Soft index decoding:** the ML detector is hard-output; soft index LLRs
  (for a coded index layer) are not implemented.
- **Frequency-selective groups:** the Rayleigh model gives independent gains per
  subcarrier; correlated gains would be the realistic next step.

### 10.3 Study gate

The original gate asked for tests and figures for CT-IM-1 through CT-IM-4 and
CT-IM-6, a revised write-up, and the stated limitations.  That gate is met for
the known-CSI, AWGN, and flat-Rayleigh scope above; the EGIM baseline and the
literature audit remain open and are listed in §10.2.  This track is still not a
core-project commitment.

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
