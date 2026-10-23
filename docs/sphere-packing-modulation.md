# Sphere Packing, Multidimensional Modulation, and the Poincaré Conjecture

### What the packing actually buys, what it cannot, and what only shares a sphere

*Scientific note for [modcode-lab](../README.md)*

*Measurement companion: [Lattice modulation](lattice-modulation.md), which
holds the BER, PAPR, and complexity tables quoted in §4.*

---

## Abstract

Multi-dimensional constellations are a real research direction, and the
sphere-packing problem in high dimensions is genuinely solved in dimensions $8$
and $24$. Both statements are theorems, and it is reasonable to ask how much of
that mathematics a modem can cash in. This note separates the three gains that
are usually added together by mistake — the nominal packing gain, the shaping
gain, and the coding gain — and puts a number on each. The conclusion is
narrow but firm: at a fixed rate and power in AWGN, the entire advantage
available from constellation *geometry* is the shaping gain, at most
$\pi e/6 = 1.53$ dB, and no constellation can be "much better" than square QAM.
The measured tables in this repository show the nominal gains of $D_4$ and $E_8$
eroding to a fraction of their asymptotic values at finite size, exactly as the
theory predicts.

The note also states, and justifies, a negative result. The sphere-packing
chain has no connection to the Poincaré conjecture, which Grigori Perelman
proved in 2002–2003. No theorem, conjecture, or heuristic links a
homeomorphism classification of three-manifolds to a packing density in
$\mathbb{R}^n$. The two subjects share a word, not a mechanism, and a paper
that asserts a correlation between them will be rejected by any competent
referee. §6 gives the invariance argument that rules the correlation out, and
the references that belong in its place.

---

## 0. Scope and status labels

| Label | Meaning |
|---|---|
| **Definition** | A convention or name; neither true nor false. |
| **Theorem** | Proved, unconditionally. |
| **Conjecture** | Open; supported by evidence and heuristics. |
| **Measured** | A number produced by the code in this repository. |
| **Analogy** | A structural resemblance; it implies nothing by itself. |
| **False** | Contradicted by an argument given here. |

Two conventions are used throughout:

1. **Dimensions are real.** A complex symbol is two real dimensions. A scheme
   with $M$ points in $n$ real dimensions carries $\eta = \log_2 M / n$ bits
   per real dimension, and every comparison below holds $\eta$ and the energy
   per real dimension fixed.
2. **Gains are in dB.** A "gain of $x$ dB" means the ratio in received
   signal-to-noise ratio needed for the same error rate, not a change in
   throughput.

### 0.1 The claim under test

> **Claim.** A sphere-like constellation in more dimensions is *much* better
> than QAM, its advantage is the sphere-packing problem in higher dimensions,
> and the advantage is connected to Perelman's proof of the Poincaré
> conjecture.

The verdict is mixed, and the three parts come apart cleanly:

| Part of the claim | Verdict | Where |
|---|---|---|
| Packing in high dimensions determines the best possible gain | **Theorem**, with a ceiling | §2, §3 |
| A sphere-like constellation is *much* better than QAM | **False** as stated; bounded by 1.53 dB | §3 |
| The gain correlates with the Poincaré conjecture | **False**; no implication exists | §6 |

The useful part of the intuition survives: shaping a constellation to a
sphere-like region is a real, measurable, standard gain, and it is exactly the
gain that dense lattices provide. The rest of the note is about how large that
gain is, and where it does and does not appear.

---

## 1. What a constellation actually optimises

### 1.1 The object

A modulation is a finite set $S \subset \mathbb{R}^n$ of $M$ points together
with a labelling $S \to \{0,1\}^{\log_2 M}$ and a decoder. Its power efficiency
at fixed $\eta$ depends on four things, not one:

- the **minimum distance** $d_{\min}$, which sets the high-SNR slope;
- the **number of nearest neighbours** $N_{\min}$, which multiplies the error
  probability at that slope;
- the **shape of the region** the points fill, which sets the shaping gain;
- the **labelling and decoder**, which decide how much of $d_{\min}$ is
  actually realised.

A note that reports only $d_{\min}$ is reporting one of four numbers.

### 1.2 Three gains, not one

**Definition (nominal coding gain).** For a lattice
$\Lambda \subset \mathbb{R}^n$ with minimum squared norm $d^2$ and fundamental
cell volume $\det \Lambda$, the Hermite parameter is

$$
\gamma(\Lambda) = \frac{d^2}{(\det \Lambda)^{2/n}},
$$

and the nominal coding gain over the integer lattice $\mathbb{Z}^n$ is
$10 \log_{10} \gamma(\Lambda)$.

**Definition (shaping gain).** For a constellation region $R$ with normalized
second moment $G(R)$ — the second moment about its centroid, divided by
$(\mathrm{vol}\,R)^{1+2/n}$ — the shaping gain over a cube of equal volume is

$$
\gamma_s(R) = \frac{1/12}{G(R)} .
$$

The cube has $G = 1/12$ by definition, so it is the reference.

**Definition (coding gain).** Everything beyond $d_{\min}$ and shaping that
comes from the channel code and its decoder. This is not a property of the
constellation at all.

The three are independent, and only the first two are geometry. That
distinction is the whole argument of §3.

---

## 2. The packing that does enter

### 2.1 Lattices and their packing figures

A lattice is a discrete additive subgroup of $\mathbb{R}^n$. Its Voronoi cell
is the set of points closer to the origin than to any other lattice point; it
tiles $\mathbb{R}^n$ by translation, and its volume is $\det \Lambda$. The
density of the associated sphere packing is

$$
\delta(\Lambda) = \frac{\rho^n}{\det \Lambda},
\qquad \rho = \frac{d}{2},
$$

which is the fraction of space covered by spheres of radius $\rho$.

These are the standard lattices, with the figures computed by
`examples/sphere_packing_lab.py`:

| Lattice | $n$ | $d^2$ | $\det$ | $\gamma$ | Nominal gain | $\delta$ |
|---|---:|---:|---:|---:|---:|---:|
| $\mathbb{Z}^2$ | 2 | 1.0 | 1.0 | 1.0000 | 0.000 dB | 0.25000 |
| $D_4$ | 4 | 2.0 | 2.0 | 1.4142 | 1.505 dB | 0.12500 |
| $E_8$ | 8 | 2.0 | 1.0 | 2.0000 | 3.010 dB | 0.06250 |
| Leech | 24 | 4.0 | 1.0 | 4.0000 | 6.021 dB | 1.00000 |

The pattern is the source of the intuition: the packing density *improves*
with dimension, and the nominal gain rises with it. It is also why $E_8$ and
the Leech lattice dominate discussions of high-dimensional modulation.

### 2.2 The theorem at the end of the chain

The sphere-packing problem is solved in exactly four dimensions, and each
solution is a theorem about **asymptotic density**:

| Dimension | Densest packing | Status |
|---:|---|---|
| 2 | hexagonal lattice | **Theorem** (Thue 1910) |
| 3 | face-centred cubic | **Theorem** (Hales 2005) |
| 8 | $E_8$ lattice | **Theorem** (Viazovska 2017) |
| 24 | Leech lattice | **Theorem** (Cohn–Kumar–Miller–Radchenko–Viazovska 2017) |

Every one of these proofs is analysis or number theory: the Cohn–Elkies linear
programming bound, and then a magic function built from modular forms that
attains the bound in dimensions $8$ and $24$. None of them is topology, and
none of them says anything about a constellation of $256$ points.

This is the first place the claim under test goes wrong in practice. A theorem
of the form "no packing in dimension $24$ is denser than Leech" constrains the
limit of an infinite sequence of packings. A 256-point spherical cut of $E_8$
is not $E_8$; it is a finite subset with its own $d_{\min}$, its own
neighbour count, and its own rate. §4 measures the difference.

### 2.3 The limit that the lattice does reach

There is a genuine and quantitative chain, and it is worth stating precisely
because it is the one the claim is reaching for.

**Theorem (Poltyrev 1994).** For the unconstrained-power AWGN channel, there is
a limit on the logarithm of the lattice volume per dimension below which
reliable decoding is possible, and lattices with vanishing second moment can
approach it.

**Theorem (Erez–Zamir 2004).** Nested lattice codes with dithering achieve the
capacity of the AWGN channel, $\tfrac12 \log_2(1+\mathrm{SNR})$ per real
dimension, with lattice decoding.

So dense lattices do connect to the theoretical limit — but through *coding*,
as the infinite-dimensional codebook, and with a decoder that does lattice
decoding. The connection is real; it is also not a statement that a
sphere-shaped constellation is much better than QAM at the same size. It says
that a lattice *code* can be as good as any code, which is what a capacity
proof needs.

---

## 3. The ceiling: 1.53 dB

### 3.1 Shaping

Take a uniform source over a region of fixed volume. The sphere minimises the
second moment, so it maximises the shaping gain over the cube. The bound is
reached in the limit of large dimension and is not exceedable:

$$
\gamma_{s,\max} = \frac{1/12}{1/(2\pi e)} = \frac{\pi e}{6} = 1.4233,
$$

or $1.53$ dB. Finite-dimensional regions get a fraction of it:

| Region | $G$ | Shaping gain |
|---|---:|---:|
| cube (QAM) | 0.0833333 | 0.000 dB |
| $D_4$ Voronoi cell | 0.0766032 | 0.366 dB |
| $E_8$ Voronoi cell | 0.0716821 | 0.654 dB |
| Leech Voronoi cell | 0.0657710 | 1.028 dB |
| sphere (limit) | 0.0585498 | 1.533 dB |

The $G$ values are Conway and Sloane's table; the gains are computed here.
This table is the honest version of "spheres are better": the best conceivable
shaping ever recovers $1.53$ dB, and the best lattice in dimension $24$ recovers
$1.03$ dB of it.

### 3.2 The gap that shaping is supposed to close

Uncoded square QAM at a bit error rate of $10^{-5}$ sits roughly $7.5$ dB above
the Shannon limit:

| $\eta$ | Scheme | Shannon limit | Uncoded at $10^{-5}$ | Gap |
|---:|---|---:|---:|---:|
| 1 | QPSK | 1.761 dB | 9.588 dB | 7.827 dB |
| 2 | 16-QAM | 5.740 dB | 13.435 dB | 7.694 dB |
| 3 | 64-QAM | 10.212 dB | 17.787 dB | 7.575 dB |
| 4 | 256-QAM | 15.035 dB | 22.503 dB | 7.469 dB |

The uncoded column is the standard closed-form Gray reference implemented in
`modcode.metrics.theory`; the limit column is
$(2^{2\eta}-1)/(2\eta)$, the real-input AWGN capacity inverted for $\eta$ bits
per real dimension. Both columns are reproduced by
`examples/sphere_packing_lab.py`.

**Measured (arithmetic).** Shaping can remove at most $1.53$ dB of a
$\approx 7.5$ dB gap. The remaining $\approx 6$ dB must come from coding — and
a binary code delivers it just as well on QAM as on a lattice. This is why
every serious comparison in the literature separates the two, and why no
constellation can be "much better" than QAM at fixed rate and power.

**Status: Theorem.** The $1.53$ dB ceiling is a bound over all regions of
fixed volume, not an engineering limitation.

### 3.3 What would have to be false for the claim to hold

For a sphere-like constellation to beat QAM by more than a couple of dB in
AWGN, one of these would have to fail:

- capacity as the upper bound on achievable rate at a given SNR;
- the isoperimetric-type bound $G \ge 1/(2\pi e)$;
- the fact that a channel code is blind to the shape of the constellation it
  is applied to.

All three are theorems. The claim, as stated, is therefore false, and the
interesting question is not whether it can be evaded but where the geometry
does something the bound does not capture — §5.

---

## 4. Why nominal gains erode at finite size

A finite spherical cut is not the lattice. Three effects consume the nominal
gain, and the repository measures all three.

### 4.1 The $d_{\min}$ accounting at matched $\eta$

From the measurement note, at unit energy per real dimension:

| Constellation | $\eta$ | $d_{\min}$ | Nominal gain |
|---|---:|---:|---:|
| QPSK | 1 | 2.000 | reference |
| $D_4$-16 | 1 | 2.000 | 0 dB |
| $E_8$-256 | 1 | 2.744 | 2.75 dB |
| 16-QAM | 2 | 0.894 | reference |
| $D_4$-256 | 2 | 1.069 | 1.55 dB |

The $D_4$-16 row is the useful warning. Its raw $d_{\min}$ looks better than
16-QAM's only because the comparison is dimensionally unfair: $D_4$-16 spends
four real dimensions on four bits, so it is a QPSK-class format, not a
16-QAM-class one. At matched $\eta$ it ties QPSK exactly.

### 4.2 The measured BER

From the same note, hard-decision BER at matched $\eta$ and equal energy per
dimension:

| $E_b/N_0$ | QPSK | $D_4$-16 | $E_8$-256 | 16-QAM | $D_4$-256 |
|---:|---:|---:|---:|---:|---:|
| 4 dB | 1.26e-2 | 3.09e-2 | 1.71e-2 | 5.94e-2 | 1.44e-1 |
| 6 dB | 2.37e-3 | 6.78e-3 | 1.11e-3 | 2.79e-2 | 6.74e-2 |
| 8 dB | 1.40e-4 | 6.20e-4 | 0 | 9.35e-3 | 2.00e-2 |
| 12 dB | 0 | 0 | 0 | 1.10e-4 | 8.00e-5 |

Three readings matter:

- The $D_4$-16 format **loses to QPSK** despite equal $d_{\min}$, because its
  nearest neighbours are more numerous and its labelling is not Gray-optimal.
  Equal $d_{\min}$ is not equal performance.
- The $E_8$-256 format **beats QPSK by about $0.8$ dB at $10^{-3}$**, read off
  the $6$ dB row, against a nominal $d_{\min}$ advantage of $2.75$ dB. Roughly
  $30\%$ of the nominal gain survives.
- The $D_4$-256 format **only catches 16-QAM between $11$ and $12$ dB**, and
  loses by $\approx 0.8$ dB at $8$ dB, against a nominal $1.55$ dB advantage.

### 4.3 The residual gap

At $10^{-3}$ the $E_8$-256 format needs about $6$ dB, while the Shannon limit at
$\eta = 1$ is $1.76$ dB. The residual gap is therefore about $4.2$ dB. Geometry
narrowed the $7.8$ dB gap of uncoded QPSK; it did not close it, and the
remaining $4.2$ dB is a coding problem.

### 4.4 Why the erosion is expected

The union bound for a constellation with minimum distance $d_{\min}$ and
$N_{\min}$ nearest neighbours has the form

$$
P_e \;\lesssim\; N_{\min}\, Q\!\left(\frac{d_{\min}}{2\sigma}\right),
$$

so two constellations with the same $d_{\min}$ differ by
$10\log_{10}(N_{\min,2}/N_{\min,1})$ dB. A finite spherical cut has many
nearest neighbours — that is what "dense packing" means locally — and the
labelling of those neighbours is not free. The nominal gain is an asymptotic
quantity that assumes an ML decoder and a labelling that realises the distance
structure.

**Status: Measured** for the tables, **Theorem** for the bound form.

---

## 5. Where multidimensional modulation genuinely wins

The geometry is not useless; it is simply bounded in AWGN. Four regimes where
it is not bounded by the same number:

### 5.1 Fading: signal-space diversity

In a Rayleigh channel the gain is not packing density but **diversity**. A
lattice rotation that spreads each real symbol across all coordinates gives
full diversity order, and the gain over an unrotated QAM is real, large, and
demonstrable. This is a channel effect: the noise is no longer isotropic in the
relevant coordinates, so region shape is not the limiting quantity.

**Status: Theorem** (Boutros–Viterbo), and testable in this repository's
`channel/rayleigh.py`.

### 5.2 Coherent optical: 4-D and 8-D formats

In coherent optical links the four real dimensions of the Jones space carry the
signal, and 4-D formats are standard research. The reported gains over
polarisation-multiplexed QPSK are of order $1$ dB, and a good part of that
comes from nonlinearity tolerance and set partitioning rather than from packing
density. This is the most defensible place for a "more-dimensional modulation"
claim, precisely because the channel is not AWGN.

**Status: Measured** in the literature; the modelling is not in this repository.

### 5.3 Peak-power-limited links

The repository measures a PAPR penalty of $0.3$ to $0.9$ dB for the lattice
formats in its component-stream diagnostic. Packing under a *crest* constraint
is not solved by the classical packing theory, which assumes an average-power
sphere, and the crest-constrained problem is genuinely open.

**Status: Open** for the design problem; **Measured** for the diagnostic.

### 5.4 Short blocks and non-AWGN channels

Asymptotic packing arguments assume large dimension and long blocks. At short
blocklengths the finite-blocklength expansion dominates, and at high phase
noise or Doppler the Euclidean isotropy that makes spheres optimal is the wrong
model. This repository already contains the channels to test the last point:
`channel/phase_imbalance.py`, `channel/multipath.py`, and the optical IM/DD
path.

**Status: Open.**

---

## 6. Why the Poincaré conjecture is not in this chain

This section exists because the claim is made often enough to be worth
answering properly, and because the answer is short and decisive.

### 6.1 What the theorem says

**Theorem (Poincaré conjecture; Perelman 2002–2003).** Every closed,
simply connected three-manifold is homeomorphic to the three-sphere $S^3$.

Perelman proved it as part of geometrization, using Ricci flow with surgery —
a controlled deformation of a Riemannian metric that pinches off
high-curvature regions. The statement, like the proof, is about
**homeomorphism classes of three-manifolds**. It is a classification result.

### 6.2 The invariance argument

Suppose, for contradiction, that the Poincaré conjecture had a consequence for
packing. Any such consequence would have to be stated in terms of a
homeomorphism-invariant quantity, because the theorem concludes only that two
spaces are homeomorphic.

But every quantity in a packing or modulation problem is **not**
homeomorphism-invariant. Packing density, minimum distance, second moment, and
capacity all depend on the metric structure of $\mathbb{R}^n$: a homeomorphism
can take a ball to a shape of arbitrarily large or small volume, and can move
points arbitrarily far apart. A theorem that says "these two spaces have the
same topology" therefore cannot imply anything about which of two
constellations is denser, closer, or more power-efficient.

This is not a gap in the literature that a paper could fill. It is the wrong
kind of statement to imply the right kind of conclusion.

### 6.3 The dimensional mismatch

The packing problems that modulate a constellation live in dimensions $2$, $4$,
$8$, and $24$ — the dimensions where good lattices exist, and the dimensions
where the packing problem has been solved. The Poincaré conjecture is about
dimension $3$ and only dimension $3$. The two subjects do not even meet in
dimension count, and the one place they might have met, the Kepler conjecture
in dimension $3$, was settled by Hales with a computer-assisted exhaustion of
cases, not by topology.

### 6.4 What the real geometric chain is

For completeness, the chain that does exist, and that a paper on this topic
should cite:

```
sphere packing (asymptotic density)
  ⇐ linear programming bounds (Cohn–Elkies 2003)
  ⇐ modular forms and magic functions (Viazovska 2017; d = 8, 24)

dense lattice  ⇒  large Hermite parameter γ
  ⇒  nominal coding gain 10 log₁₀ γ        (asymptotic, ML decoding)

constellation region  ⇒  small G  ⇒  shaping gain (1/12)/G ≤ πe/6 = 1.53 dB

lattice code + dither  ⇒  AWGN capacity    (Erez–Zamir 2004)
```

Every arrow is proved, and every arrow has a hypothesis attached: asymptotic
dimension, ML decoding, known CSI, or uniform power. None of them has a
topological hypothesis.

### 6.5 Status

**Status: False.** There is no theorem, no conjecture, and no heuristic
connecting the Poincaré conjecture to sphere packing or to modulation
performance, and §6.2 shows that no such implication can exist in principle.
A manuscript that asserts the correlation would be rejected at review, and a
manuscript that presents it as a new research direction would damage the
author's credibility on the parts that are correct.

The correct citations for the geometric content are Conway–Sloane for the
lattices, Viazovska and Cohn–Elkies for the packing theorems, and Forney and
Erez–Zamir for the communication theory. Perelman has no place in the
bibliography of a modulation paper.

---

## 7. What is actually open

Four directions where the geometry is not already bounded, and where this
repository can produce evidence:

- **Crest-constrained packing.** Optimise a constellation for minimum
  distance under a PAPR or cubic-metric constraint. The classical theory
  assumes an average-power sphere; the crest-constrained problem is open, and
  the diagnostic in §5.3 is the departure point.
- **Finite-blocklength joint design.** Constellation, shaping, and code chosen
  together for short blocks, where the normal approximation replaces the
  asymptotic gains. The repository's `polar.py` and `ldpc.py` give the code
  side.
- **Non-AWGN geometry.** Rotation, shaping, and packing under phase noise,
  Doppler, and IM/DD optical front ends, using the channels already
  implemented.
- **Decoder complexity at matched performance.** The naive lattice demodulator
  in this repository is two orders of magnitude slower than a QAM slicer;
  whether the $0.8$ dB of §4.2 survives a decoder that fits a real budget is a
  measurable question.

A negative result is a legitimate outcome here, and this repository is
unusually well equipped to produce one honestly.

---

## 8. Reproduce

```bash
python examples/sphere_packing_lab.py      # tables in §2.1, §3.1, §3.2
python examples/lattice_lab.py             # BER, PAPR, throughput in §4
python -m pytest tests/test_sphere_packing.py -q
```

The script prints the three tables and no figure; the claims in §2 and §3 are
arithmetic on published invariants, and the claims in §4 are the measured
tables of the companion note.

---

## 9. Claims audit

| Claim | Status | Comment |
|---|---|---|
| Dense lattices give a nominal gain over $\mathbb{Z}^n$ | **Theorem** | $1.51$ dB for $D_4$, $3.01$ dB for $E_8$, $6.02$ dB for Leech |
| The best constellation region is sphere-like | **Theorem** | Isoperimetric bound $G \ge 1/(2\pi e)$ |
| The shaping gain is at most $1.53$ dB | **Theorem** | $\pi e/6$, attained only as $n \to \infty$ |
| Sphere-like modulation is *much* better than QAM | **False** | Bounded by $1.53$ dB at fixed rate and power in AWGN |
| A finite $E_8$ cut realises the $E_8$ gain | **False** | $2.75$ dB nominal, about $0.8$ dB measured at $10^{-3}$ |
| Equal $d_{\min}$ implies equal performance | **False** | Nearest-neighbour count and labelling intervene |
| Lattice codes can achieve AWGN capacity | **Theorem** | Erez–Zamir 2004, with dithering and lattice decoding |
| The packing problem is open in all high dimensions | **False** | Solved in $d = 2, 3, 8, 24$; open elsewhere |
| Packing gains correlate with the Poincaré conjecture | **False** | §6.2: no homeomorphism-invariant consequence can exist |
| Perelman's proof is relevant to modulation | **False** | Classification of three-manifolds; no quantitative content |
| Multidimensional formats win in coherent optical | **Measured** | Of order $1$ dB, largely nonlinearity and set partitioning |
| Rotation gives diversity in Rayleigh fading | **Theorem** | Boutros–Viterbo; a channel effect, not packing |
| Lattice demodulation is cheap | **False** | Two orders of magnitude slower than a QAM slicer here |

---

## References

1. C. E. Shannon, “Probability of error for optimal codes in a Gaussian
   channel,” *Bell System Technical Journal* **38** (1959), 611–656. The
   sphere-packing bound: the one place where packing is a *bound* on coding.
2. J. H. Conway and N. J. A. Sloane, *Sphere Packings, Lattices and Groups*,
   3rd edition, Springer, 1999. Lattices, Voronoi cells, and Table 2.3 for the
   normalized second moments used in §3.1.
3. G. D. Forney Jr., “Trellis shaping,” *IEEE Transactions on Information
   Theory* **38** (1992), 281–300. Shaping gain and its $1.53$ dB limit.
4. G. D. Forney Jr. and G. Ungerboeck, “Modulation and coding for linear
   Gaussian channels,” *IEEE Transactions on Information Theory* **44** (1998),
   2384–2415.
5. G. Ungerboeck, “Channel coding with multilevel/phase signals,” *IEEE
   Transactions on Information Theory* **28** (1982), 55–67. Trellis-coded
   modulation, the practical route to coding gain on QAM.
6. G. Poltyrev, “On coding without restrictions for the AWGN channel,”
   *Problems of Information Transmission* **30** (1994), 239–247. The
   unconstrained-power limit that lattice codes approach.
7. U. Erez and R. Zamir, “Achieving $1/2\log(1+\mathrm{SNR})$ on the AWGN
   channel with lattice encoding and decoding,” *IEEE Transactions on
   Information Theory* **50** (2004), 2293–2314.
8. H.-A. Loeliger, “Signal sets matched to groups,” *IEEE Transactions on
   Information Theory* **37** (1991), 1675–1682. Group structure buys
   decoding, not packing density.
9. J. Boutros and E. Viterbo, “Signal space diversity: a power- and
   bandwidth-efficient diversity technique for the Rayleigh fading channel,”
   *IEEE Transactions on Information Theory* **44** (1998), 1453–1467.
10. H. Cohn and N. Elkies, “New bounds on sphere packings I,” *Annals of
    Mathematics* **157** (2003), 689–714; arXiv:math/0110009.
11. M. S. Viazovska, “The sphere packing problem in dimension 8,” *Annals of
    Mathematics* **185** (2017), 991–1015; arXiv:1603.04246.
12. H. Cohn, A. Kumar, S. D. Miller, D. Radchenko, and M. Viazovska, “The
    sphere packing problem in dimension 24,” *Annals of Mathematics* **185**
    (2017), 1017–1033; arXiv:1603.06515.
13. T. C. Hales, “A proof of the Kepler conjecture,” *Annals of Mathematics*
    **162** (2005), 1065–1185.
14. A. Thue, “Über die dichteste Zusammenstellung von kongruenten Kreisen in
    einer Ebene,” *Norske Vid. Selsk. Skr.* **1** (1910), 1–9.
15. P. J. Karlsson and E. Agrell, “Which is the most power-efficient modulation
    format in optical links?” *Optics Express* **17** (2009), 10814–10819.
16. G. Perelman, “The entropy formula for the Ricci flow and its geometric
    applications,” arXiv:math/0211159, 2002.
17. G. Perelman, “Ricci flow with surgery on three-manifolds,”
    arXiv:math/0303109, 2003.
18. J. Morgan and G. Tian, *Ricci Flow and the Poincaré Conjecture*, American
    Mathematical Society, 2007. An exposition of the proof used in §6.
