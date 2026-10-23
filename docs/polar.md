# Polar(1024, 512): SC, SCL, and CRC-aided list decoding

*Core FEC extension. This is a pedagogical finite-length Polar code for the
BI-AWGN channel, not a 5G NR implementation.*

[`ldpc.md`](ldpc.md) adds an iteratively decoded code and compares it with the
convolutional baseline.  This note adds the other code family used in 5G: a
polar code decoded by successive cancellation (SC), by successive-cancellation
list (SCL) decoding, and by CRC-aided SCL.  All three run through the same
`Link`, so the comparison with LDPC is a controlled one.

The block length is $N = 1024 = 2^{10}$ with $k = 512$ message bits at rate
$0.5$, chosen to sit next to `LDPC36` ($N = 1008$, rate $0.502$).

---

## 1. Construction

The encoder applies the polar transform

$$
x = u F^{\otimes 10}, \qquad F = \begin{pmatrix} 1 & 0 \\ 1 & 1 \end{pmatrix},
$$

in the natural index order: no bit reversal is applied, and the information
positions are therefore not the last $k$ positions of $u$.  The recursive
butterfly used by `encode` is

$$
x_L = \mathcal{T}(u_L) \oplus \mathcal{T}(u_R), \qquad x_R = \mathcal{T}(u_R),
$$

where $(u_L, u_R)$ splits $u$ in half and $\mathcal{T}$ is the same transform on
the half-size block.  The $N-K$ frozen positions are set to zero.

The $K$ information positions are the positions with the largest synthetic
reliability under the Gaussian approximation (GA):

- the channel LLR mean is $m_0 = 2/\sigma^2$, with
  $\sigma^2 = 1/(2 R E_b/N_0)$ at the design point;
- a check-like (left) child has mean $\phi^{-1}\bigl(1 - (1-\phi(m))^2\bigr)$;
- a variable-like (right) child has mean $2m$;
- $\phi(x) = 1 - \mathbb{E}[\tanh(U/2)]$ for $U \sim \mathcal{N}(x, 2x)$.

`gaussian_phi` evaluates the expectation with a 401-point quadrature and
`gaussian_phi_inverse` inverts a fixed table by log-log interpolation.  The
design point defaults to 2 dB; `design_ebn0_db` exposes it because the
information set depends on it (§4).  The construction is approximate: GA
replaces the exact LLR density with a Gaussian, and no density evolution,
weight enumeration, or minimum-distance optimization is performed.

**CRC variant.**  With `crc=CRC8`, 520 bits are placed in the non-frozen
positions: 512 message bits followed by 8 CRC parity bits.  The message rate is
still $512/1024 = 0.5$, and because the eight CRC bits are determined by the
message, the number of codewords is $2^{512}$: the code rate is 0.5 as well.
The polar construction now protects 520 non-frozen positions instead of 512,
and the CRC gives the list decoder a selection rule that plain SCL does not
have.  The CRC is the lab's existing `CRC8` ($x^8 + x^2 + x + 1$), placed in
the last eight information positions.

`PolarBEC` in the same module remains the appendix-level erasure-channel
construction; it is not used by this decoder.

---

## 2. Decoder

`decode_soft` walks the factor tree depth-first.  At every node it forms the
two child LLR vectors

$$
L_{\text{left}} = 2\tanh^{-1}\!\left(\tanh\frac{a}{2}\,\tanh\frac{b}{2}\right),
\qquad
L_{\text{right}} = b + (1-2c)\,a,
$$

where $a$ and $b$ are the halves of the parent LLR vector and $c$ is the
already-decoded left partial-sum bit.  The check update is evaluated in the
Jacobian-logarithm form

$$
L_{\text{left}} = \operatorname{sgn}(a)\operatorname{sgn}(b)\min(|a|,|b|) +
\log\!\left(1+e^{-|a+b|}\right) - \log\!\left(1+e^{-|a-b|}\right),
$$

which is exact, numerically stable, and invariant to a positive rescaling of
the LLRs.

**SC** (`list_size=1`) decides each information bit from the sign of its leaf
LLR.  **SCL** keeps $L$ paths and updates each path metric by

$$
\mathrm{PM} \mathrel{+}= \log\!\left(1 + e^{-(1-2b)L}\right),
$$

the negative log-likelihood of the decision; after every information bit the
$L$ lowest-metric paths survive.  **CA-SCL** keeps the minimum-metric path whose
CRC checks.  If no path checks, it falls back to the minimum-metric path; the
CRC false-accept probability is about $2^{-8}$, so a wrong path can
occasionally be selected.

`decode` (hard bits) maps $0/1$ to $\pm4$ and calls the same decoder.
`last_work` counts element-level check and variable updates per block, so one SC
pass is exactly $N\log_2 N = 10{,}240$ updates; an SCL pass is roughly $L$ times
that, discounted because the list grows from one path and often saturates.

---

## 3. Results

`examples/polar_lab.py` sweeps the codes with BPSK over AWGN, 80,000 bits per
point with one fixed seed per point.  Measured BER:

| Eb/N0 | polar SC | polar SCL $L=4$ | polar CA-SCL $L=4$ | LDPC min-sum | convolutional |
|---:|---:|---:|---:|---:|---:|
| 0.0 dB | 4.43e-1 | 4.07e-1 | 4.32e-1 | 1.56e-1 | 9.21e-2 |
| 0.5 dB | 3.73e-1 | 2.30e-1 | 2.84e-1 | 1.42e-1 | 6.83e-2 |
| 1.0 dB | 2.59e-1 | 1.03e-1 | 1.16e-1 | 9.88e-2 | 5.03e-2 |
| 1.5 dB | 7.93e-2 | 9.18e-3 | 1.04e-2 | 2.58e-2 | 2.57e-2 |
| 2.0 dB | 2.20e-2 | 5.75e-4 | 2.72e-3 | 1.55e-3 | 1.76e-2 |
| 2.5 dB | 3.06e-3 | 0 | 0 | 0 | 8.23e-3 |

What the table shows:

- **SC is the weakest decoder.**  At 2 dB it is still at $2.2\times10^{-2}$
  while the LDPC code is at $1.55\times10^{-3}$.
- **List decoding is the main gain.**  At 2 dB, $L=4$ lowers the polar BER from
  $2.2\times10^{-2}$ to $5.8\times10^{-4}$, a factor of about 38 at the same
  $E_b/N_0$.
- **The two 5G families cross.**  The LDPC code is better below about 1.2 dB;
  the polar code is better from 1.5 dB upward in this pair, and both are
  error-free at 2.5 dB in this sample.  At 1.5 dB the polar SCL code reaches
  $9.2\times10^{-3}$ against $2.6\times10^{-2}$ for the LDPC code.
- **The zeros are zero-error runs**, 80,000 bits at one seed.  They are not an
  error floor, and the difference between "0" and a few errors is not resolved
  at this sample size.

The LDPC curve comes from the same `LDPC36` min-sum decoder as
[`ldpc.md`](ldpc.md); the convolutional curve is the K=3 soft Viterbi code.
Both are pedagogical implementations, not 5G or 3GPP conformant.

---

## 4. List size, block length, and construction

**List size.**  Decoder work and measured decoding time at 2 dB, 16 blocks of
1024 coded bits each:

| decoder | list | work/block | ms/block | BER at 2 dB |
|---|---:|---:|---:|---:|
| polar SC | 1 | 10,240 | 15.6 | 2.20e-2 |
| polar SCL | 2 | 18,688 | 17.8 | 5.75e-4 |
| polar SCL | 4 | 34,432 | 19.0 | 5.75e-4 |
| polar SCL | 8 | 64,900 | 20.0 | 5.75e-4 |
| polar CA-SCL, CRC-8 | 2 | 18,688 | 18.0 | 5.49e-3 |
| polar CA-SCL, CRC-8 | 4 | 34,432 | 19.0 | 2.72e-3 |
| polar CA-SCL, CRC-8 | 8 | 64,900 | 20.2 | 0 |
| LDPC min-sum | - | - | 39.0 | 1.55e-3 |
| convolutional | - | - | 1.2 | 1.76e-2 |

The millisecond column varies by a few percent between runs; the work column
is deterministic.

Two things are worth noticing.  First, **the list saturates at high SNR**: at
2 dB the plain SCL decoder produces the same BER for $L=2,4,8$, because the
surviving paths already contain the maximum-likelihood word.  A spot check on
the failing frames confirms that the decoded codeword has a higher likelihood
than the transmitted one, so the residual errors are ML errors, not list
misses.  Below about 1.5 dB the list still helps: measured with 80,000 bits at
1.5 dB, $L=1,2,4,8,16$ gives $7.93\times10^{-2}$, $2.55\times10^{-2}$,
$9.18\times10^{-3}$, $7.53\times10^{-3}$, $7.53\times10^{-3}$.

Second, **CRC-aided selection is not monotone in list size**.  At 2 dB, CA-SCL
with $L=2$ or $L=4$ is worse than plain SCL: the list is too small for the CRC
to reliably pick the right path, and a $2^{-8}$ false accept can select a wrong
one.  At $L=8$ the same CRC changes the answer: the residual ML-error frames
are replaced by CRC-valid paths, and the block is error-free in this sample.
The CRC is a selection device, not extra coding gain.

**Block length.**  SCL $L=4$ at 1.5 dB, 50,000 bits:

| $N$ | 256 | 512 | 1024 | 2048 |
|---|---:|---:|---:|---:|
| BER | 3.67e-2 | 1.17e-2 | 6.12e-3 | 5.68e-3 |

The improvement is the usual finite-length effect: the same decoder and rate get
better as the block grows, with diminishing returns by $N=1024$.  A comparison
against `LDPC36` at a *different* block length (for example $N=256$) would
confound code family with block length, which is why the main table uses
$N=1024$.

**Construction.**  SCL $L=4$, $N=1024$, 50,000 bits at 1.5 dB, with the
information set designed at three different $E_b/N_0$:

| design point | 1.0 dB | 2.0 dB | 3.0 dB |
|---|---:|---:|---:|
| BER at 1.5 dB | 1.30e-2 | 5.92e-3 | 1.28e-3 |

The design point is a real parameter, not a formality: at 1.5 dB the 3 dB design
has about 4.6 times fewer errors than the 2 dB design, and the 1 dB design is
worse than both.  No attempt is made to optimize it, and no claim is made that
3 dB is optimal for other rates or block lengths.

---

## 5. Complexity and latency

The decoder is $O(N\log_2 N)$ per path and $O(LN\log_2 N)$ element operations
per block, with $L$ up to 16 in the experiments here.  The block latency is
$N$ bits and the list decoder is inherently serial in the bit index; the
convolutional code remains the low-latency option and the LDPC code spreads its
work over iterations.

The measured times in §4 are Python on one machine.  They are useful for the
*relative* picture — list decoding costs a few milliseconds more per block while
the LDPC min-sum decoder costs about twice as much in this implementation — and
they are not a hardware or throughput claim.

The decoder is a straightforward list decoder: no lazy-copy path merging, no
SIMD or lookup-table node processing, and no simplified SC (SSC) fast nodes.
A production decoder for the same code would be considerably cheaper; the
relative comparison here is between two plain Python implementations.

---

## 6. What is claimed, measured, and not claimed

**Tested**

- the GA construction reproduces the exact BEC reliability order for $N=8$;
- SC agrees with an independent recursive decoder on random LLR vectors;
- noiseless hard and soft round trips for SC, SCL, and CA-SCL;
- one SC pass is exactly $N\log_2 N$ element updates, and SCL work grows with
  the list size;
- decisions are invariant to a positive rescaling of the LLRs;
- SCL reduces errors on a fixed noisy sample, and a 3 dB link BER is below
  $10^{-3}$;
- invalid lengths, non-finite LLRs, complex LLRs, and non-power-of-two blocks
  raise.

**Measured**

- the finite-length BER table and the SC/SCL gain above;
- list-size saturation at 2 dB and the ML-consistency of the residual errors;
- CRC-aided selection improving at $L=8$ but not at $L=2,4$;
- block-length and construction-point sensitivity;
- relative decoder cost for these Python implementations.

**Not claimed**

- 5G NR conformance: no rate matching, no distributed CRC, no interleaving,
  no polar sequence table, no agreement with 3GPP TS 38.212;
- that polar decoding beats LDPC in general, or that this ordering holds at
  other rates, channels, or block lengths;
- an error floor or an asymptotic threshold;
- optimality of the construction point, list size, or CRC polynomial;
- a hardware, throughput, or latency claim beyond the Python measurements;
- that SCL or CA-SCL is "better" at every SNR: CA-SCL is worse than plain SCL
  at 2 dB with $L=2,4$.

---

## 7. Related families that are not in this lab

Three modern families are frequently compared with LDPC in the literature.
None of them is implemented here, and this note does not reproduce their
reported gains.

- **Spatially coupled LDPC (SC-LDPC).**  Under belief propagation, spatially
  coupled chains exhibit threshold saturation: the BP threshold approaches the
  MAP threshold of the underlying ensemble, a gain of fractions of a dB.  It
  needs long chains and high decoding latency, and demonstrating it requires
  density-evolution style analysis and long-block simulations.  A short
  block code like `LDPC36` cannot show it.
- **GRAND and soft-output GRAND.**  Guessing random additive noise decoding
  searches noise patterns in order of likelihood instead of exchanging messages
  on a graph.  Product-code variants with soft outputs have been reported to be
  competitive with 5G LDPC; that is a different code construction and a
  different decoder, and this lab does not evaluate it.
- **5G NR itself.**  The standard uses polar codes for control channels and
  LDPC for data, with rate matching and interleaving designed around those
  choices.  Nothing in this repository implements the 3GPP chains.

The honest summary is the one this note measures: at $N \approx 1000$ and rate
$1/2$, a short polar code with list decoding crosses a min-sum LDPC code near
1.2 dB, and the ordering at other lengths or rates is a separate experiment.

---

## 8. Reproduce

```bash
python examples/polar_lab.py
python -m pytest tests/test_polar_awgn.py -q
```

Figure: Figure 24 in [`figures.md`](figures.md); source
`examples/figures/polar.png`.  The LDPC side of the comparison is documented in
[`ldpc.md`](ldpc.md).

---

## 9. References

- E. Arıkan, “Channel polarization: A method for constructing capacity-achieving
  codes for symmetric binary-input memoryless channels,” *IEEE Transactions on
  Information Theory* **55** (2009), 3051–3073.
- I. Tal and A. Vardy, “List decoding of polar codes,” *IEEE Transactions on
  Information Theory* **61** (2015), 2213–2226.
- K. Niu and K. Chen, “CRC-aided decoding of polar codes,” *IEEE Communications
  Letters* **16** (2012), 1668–1671.
- A. Balatsoukas-Stimming, M. B. Parizi, and A. Burg, “LLR-based successive
  cancellation list decoding of polar codes,” *IEEE Transactions on Signal
  Processing* **63** (2015), 5165–5179.
- 3GPP TS 38.212, *Multiplexing and channel coding*.
