# The phase-shift involution: a relative sector phase and its I/Q counterpart

*Research note for [modcode-lab](../README.md), with a companion note in the
`riemann-framework` repository.*

> **Scope:** this is an interpretive appendix, not part of the core
> modulation/FEC simulator. It records a structural analogy and an
> identifiability result; it is not a new modulation scheme and not evidence
> for RH. See [`exploratory-labs.md`](exploratory-labs.md).

---

## Abstract

The DSIN analysis uses a relative phase flip between a bosonic and a
fermionic sector,

$$
\tau = \begin{pmatrix} I & 0 \\ 0 & -I \end{pmatrix},
\qquad \tau^2 = I.
$$

It is a unitary involution. It is **not** the sector-swap involution
$\sigma$, and its fixed locus is not the equal-weight subspace. The
equal-weight symmetric states are the fixed locus of $\sigma$.

At $\phi=\pi$, an I/Q phase imbalance produces the same real $2\times2$
matrix, $\mathrm{diag}(1,-1)$, after the two sectors are identified
with the I and Q axes. That is a genuine structural correspondence, but it is
not an identity of the two models. The logical action of the matrix depends on
the encoding basis. In DSIN's $\sigma$-eigenbasis, $\tau$ flips the encoded
bit; in the standard I/Q basis of QPSK, it flips only the Q bit.

The engineering counterpart of DSIN blindness is an **identifiability**
statement, not a claim that a coherent receiver cannot see a reflection. If a
constellation and its data prior are invariant under the reflection, the two
hypotheses $\phi=0$ and $\phi=\pi$ have identical marginal likelihoods for
every observation. No blind detector can distinguish them better than chance.
A pilot, a known preamble, an asymmetric constellation, or a code constraint
breaks the symmetry. That is the same structural lesson as DSIN: an observable
that is invariant under a transformation cannot detect that transformation.

A runnable experiment in this repository demonstrates all four cases:
standard coherent detection, an oracle that knows $\phi=\pi$, a
pilot-assisted detector, and a blind detector. The blind detector is pinned at
BER $1/2$; the oracle and pilot recover. The experiment does **not** claim a
new modulation scheme, a new blind I/Q estimator, or evidence for the Riemann
Hypothesis.

---

## 1. The DSIN object

DSIN works on a two-sector Hilbert space

$$
\mathcal H = \mathcal H_b \oplus \mathcal H_f,
$$

with two distinct involutions.

### 1.1 The sector swap

The sector-swap involution is

$$
\sigma = \begin{pmatrix} 0 & I \\ I & 0 \end{pmatrix},
\qquad
\sigma |b,k\rangle = |f,k\rangle,
\qquad
\sigma |f,k\rangle = |b,k\rangle.
$$

Its eigenvalues are $+1$ and $-1$, with eigenspaces

$$
|+\rangle = \frac{|b,k\rangle + |f,k\rangle}{\sqrt2},
\qquad
|-\rangle = \frac{|b,k\rangle - |f,k\rangle}{\sqrt2}.
$$

The fixed subspace is

$$
\mathrm{Fix}(\sigma)
= \{(u,u)/\sqrt2 : u \in \mathcal H_b\},
$$

the equal-weight symmetric subspace. The antisymmetric subspace is the
$-1$ eigenspace.

### 1.2 The relative phase flip

The object called $\tau$ in the DSIN analysis is the fermionic phase flip

$$
\tau = \begin{pmatrix} I & 0 \\ 0 & -I \end{pmatrix},
\qquad
\tau |b,k\rangle = |b,k\rangle,
\qquad
\tau |f,k\rangle = -|f,k\rangle.
$$

It satisfies

$$
\tau^2 = I,
\qquad
\tau^\dagger = \tau,
\qquad
\tau^\dagger \tau = I.
$$

It is a unitary involution, and it is Hermitian, so it is also an observable
with eigenvalues $\pm1$. Its fixed subspace is the bosonic sector,

$$
\mathrm{Fix}(\tau)
= \{(u,0) : u \in \mathcal H_b\},
$$

not the equal-weight subspace.

### 1.3 The two involutions anticommute

A direct multiplication gives

$$
\sigma\tau = -\tau\sigma.
$$

Consequently $\tau$ maps the $+1$ eigenspace of $\sigma$ onto the $-1$
eigenspace:

$$
\tau |+\rangle = |-\rangle,
\qquad
\tau |-\rangle = |+\rangle.
$$

If DSIN encodes a classical bit in the $\sigma$ eigenvalue, then $\tau$
inverts the bit. It also changes the sign of $\langle\sigma\rangle$, but
leaves its magnitude unchanged:

$$
\langle\sigma\rangle \longmapsto -\langle\sigma\rangle,
\qquad
|\langle\sigma\rangle| \longmapsto |\langle\sigma\rangle|.
$$

This is the DSIN blindness result. The parity **meter** thresholds the
deviation of $|\langle\sigma\rangle|$ from 1, so it sees nothing, while the
decoded data are inverted. The data decoder is not blind: it measures
$\sigma$ and returns the opposite bit. The blind object is the
symmetry-violation monitor.

### 1.4 Fixed loci are different

The two statements that are easy to conflate are:

| Involution | Action | Fixed locus |
|:---|:---|:---|
| $\sigma$ | swaps the sectors | equal-weight symmetric subspace |
| $\tau$ | multiplies the fermionic sector by $-1$ | bosonic sector |

The sentence “the fixed locus is the set of states with equal weight in both
sectors” describes $\sigma$, not $\tau$. For $\tau$, a state
$\alpha|b\rangle+\beta|f\rangle$ is fixed exactly when $\beta=0$.

---

## 2. What the phase-shift involution is not

### 2.1 It is not the sector swap

$\sigma$ and $\tau$ are different matrices. They anticommute rather than
commute, and their fixed loci are different. A phase flip is not a swap.

### 2.2 It is not the Riemann functional equation

The functional equation acts on the complex variable as
$s \mapsto 1-s$, and the completed zeta function is invariant under it. The
prototype dimension-shift involution is
$\sigma(s)=1-\overline{s}$, whose fixed locus is $\mathrm{Re}s=1/2$.
Neither map is the relative phase between two quantum sectors. They are
structural analogues, not the same object.

### 2.3 It is not a continuous phase rotation

A phase rotation $e^{i\theta}$ is a one-parameter family. $\tau$ is the
discrete element at $\theta=\pi$ acting on one sector only. It is an
involution, not a generator of a continuous flow.

### 2.4 It is not automatically a logical bit flip

$\tau$ flips the logical bit only in a basis in which the logical states are
the $\sigma$ eigenstates. In the sector basis, $\tau$ is diagonal and
$\sigma$ is the swap. In the $\sigma$ eigenbasis, the roles reverse:

$$
\sigma \longmapsto \begin{pmatrix} 1 & 0 \\ 0 & -1 \end{pmatrix},
\qquad
\tau \longmapsto \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}
$$

on the $\{|+\rangle,|-\rangle\}$ basis. This basis dependence is the source
of most of the confusion in the engineering analogy.

### 2.5 It is not a standard named object in the RH literature

Relative phase gates, fermionic parity $(-1)^F$, and $\mathbb Z_2$-graded
phase flips are standard objects in quantum information and supersymmetry.
The specific name “phase-shift involution” is not a standard term of art.
The literature search in the proposal did not establish a named RH tool with
this action; the honest conclusion is that $\tau$ is a natural
representation-theoretic object of the DSIN model, not an established
number-theoretic instrument. No novelty claim is made here.

---

## 3. The I/Q channel, correctly

### 3.1 The model

For a complex baseband symbol $x=I+jQ$, a phase error $\phi$ in the
quadrature local oscillator gives

$$
y = I + j Q e^{j\phi}
  = (I - Q\sin\phi) + j(Q\cos\phi)
$$

before additive noise. In real vector form, with column vector $(I,Q)^T$,

$$
R(\phi) = \begin{pmatrix}
1 & -\sin\phi \\
0 & \cos\phi
\end{pmatrix},
\qquad
\det R(\phi) = \cos\phi.
$$

The common shorthand

```text
y = I + exp(j phi) Q
```

is missing the $j$ that makes $Q$ the imaginary axis. It gives
$y=I-Q$ at $\phi=\pi$, a rank-one projection onto the real axis, not a
reflection and not an involution.

The noise variance in the experiment is normalized to the **input** symbol
energy, as in the AWGN channel. The impairment is not renormalized afterwards.
At $\phi=\pi$ the map preserves energy, so the comparison is fair; at other
phases the mean symbol energy can change, and that change is part of the
impairment rather than a separate gain control.

### 3.2 The three distinguished phases

| $\phi$ | $R(\phi)$ | Action | Invertible? |
|:---:|:---|:---|:---:|
| $0$ | $\mathrm{diag}(1,1)$ | identity | yes |
| $\pi/2$ | `[[1, -1], [0, 0]]` | Q collapses onto I | **no** |
| $\pi$ | $\mathrm{diag}(1,-1)$ | $y=\overline{x}$, Q reflection | yes |

At $\phi=\pi$,

$$
R(\pi)^2 = I,
$$

so the channel is an involution and preserves energy. At $\phi=\pi/2$, the
matrix is rank one: the Q component is lost and no inverse exists. The
singular case is $\pi/2$ (mod $\pi$), not $\pi$.

### 3.3 What a standard QPSK receiver sees

For Gray-mapped QPSK with bits $(b_I,b_Q)$, the reflection maps

$$
(b_I,b_Q) \longmapsto (b_I,1-b_Q).
$$

Every symbol is mapped to a different symbol, so the symbol error rate is
$1$. Exactly one of the two bits is wrong, so the bit error rate is
$1/2$, not $1$. For square QAM the reflection maps each Q-axis level to
its mirror; the Q-axis label is permuted and its sign bit flips, while the
I-axis bits survive. The exact BER depends on the Gray labeling, but it is
not zero.

The statement “the receiver decodes every QPSK or QAM symbol incorrectly” is
therefore true at the symbol level and false at the bit level. The experiment
reports BER, so it reports $1/2$.

### 3.4 Why $\pi$ is not simply “the worst case”

The determinant and the bit error rate measure different things:

- at $\pi/2$ the channel is rank deficient and loses the Q component, so
  the **information capacity** is reduced and the inverse does not exist;
- at $\pi$ the channel is invertible and energy preserving, but the
  reflection inverts the Q bit, so a slicer that assumes $\phi=0$ has
  BER $1/2$.

Which is “worse” depends on the metric. The claim that $\phi=\pi$ is the
worst case of I/Q phase imbalance is not supported without naming the metric.
The experiment prints both the determinant and the BER so the distinction is
visible.

---

## 4. Same matrix, different logical action

If the bosonic and fermionic sectors are identified with the I and Q axes,
then

$$
\tau = \mathrm{diag}(1,-1) = R(\pi).
$$

This is the precise algebraic correspondence: the same real $2\times2$
matrix. It is not the same complex-linear operator. On $\mathbb C$, the I/Q
map at $\pi$ is conjugation, $x\mapsto\overline{x}$, which is
anti-linear; DSIN's $\tau$ is complex-linear on
$\mathcal H_b\oplus\mathcal H_f$. The agreement is real-linear, after a
choice of real basis.

The logical action depends on the code:

- **DSIN encoding:** logical states are the $\sigma$ eigenstates
  $|+\rangle,|-\rangle$. The matrix $\mathrm{diag}(1,-1)$ swaps them,
  so the logical bit flips with probability one.
- **Standard QPSK:** logical states are the four I/Q sign combinations. The
  same matrix flips the Q bit only.
- **A rotated BPSK pair:** if the two logical states are
  $(|I\rangle\pm|Q\rangle)/\sqrt2$, the reflection swaps them, so the single
  logical bit flips with probability one.

Thus “the phase-shift involution is a bit flip” is true only relative to a
chosen encoding basis. The matrix is the same; the code is not.

---

## 5. Blindness is an identifiability theorem

### 5.1 The theorem

Let $\mathcal C\subset\mathbb C$ be a finite constellation, let $p(c)$ be
a prior on it, and let $R$ be the reflection $R(c)=\overline{c}$. Consider

$$
H_0:\ y = c+n,
\qquad
H_1:\ y = R(c)+n,
$$

with circular complex Gaussian noise of variance $\sigma^2$ per real
component. If

$$
R(\mathcal C)=\mathcal C
\qquad\text{and}\qquad
p(R(c))=p(c)
$$

for every $c$, then for every observation $y$,

$$
p(y\mid H_0) = p(y\mid H_1).
$$

*Proof.* Starting from $H_1$,

$$
p(y\mid H_1)
= \sum_{c\in\mathcal C} p(c)\,N(y-R(c);\sigma^2).
$$

Change variables to $c'=R(c)$. Since $R$ is an involution and
$R(\mathcal C)=\mathcal C$, the sum runs over the same set, and the prior is
invariant:

$$
p(y\mid H_1)
= \sum_{c'\in\mathcal C} p(c')\,N(y-c';\sigma^2)
= p(y\mid H_0).
$$

$\square$

The likelihood ratio is therefore identically zero. No threshold, no
“syndrome-like” function of the constellation, and no amount of data can
distinguish the hypotheses. For QPSK and square QAM with uniform data, the
conditions hold, so the blind ambiguity is exact.

### 5.2 What breaks the symmetry

The theorem fails as soon as the prior or the constellation is not invariant:

- **A pilot or known preamble.** A known symbol $c_p$ is not averaged over
  the constellation. The likelihoods under $H_0$ and $H_1$ differ by the
  distance between $c_p$ and $R(c_p)$.
- **An asymmetric constellation or prior.** If the data distribution is
  skewed, the reflection changes it.
- **A code constraint.** A CRC or FEC syndrome is an external constraint. The
  receiver can try both hypotheses and keep the one whose codeword is valid.
  This is legitimate side information, but it is not a property of Gray
  mapping or of the constellation's cyclic symmetry.
- **A second reference.** A pilot tone, a known phase reference, or a
  differential encoding across a known boundary can resolve the sign.

`tests/test_phase_imbalance.py::test_crc_constraint_breaks_the_blind_tie`
demonstrates the code-constraint case with CRC-8: at $\phi=\pi$, the
hypothesis $\phi=0$ leaves a nonzero CRC remainder and $\phi=\pi$ passes.
The code is doing the work; the constellation symmetry alone does not.

Gray mapping does not break the symmetry. It is a labeling of the
constellation; the reflection permutes the labels, and the likelihood is
unchanged. “Cyclic shift invariance” of the symbol set is likewise a symmetry,
not a syndrome.

### 5.3 The DSIN analogue

In DSIN, the parity meter observes $|\langle\sigma\rangle|$. The phase flip
maps $+$ to $-$, so the meter's statistic is invariant while the data are
inverted. In the I/Q problem, a blind detector observes only the symmetrized
likelihood, which is invariant under the reflection. The structural statement
is the same:

> An observable that is invariant under a transformation cannot detect that
> transformation. A reference state that is not invariant can.

This is why the engineering fix is a pilot or a code constraint, not a more
elaborate function of the invariant statistic. It is also why the DSIN
single-observable design cannot bound the eavesdropper's information from that
observable alone.

---

## 6. The experiment

### 6.1 What is implemented

- `modcode.channel.PhaseImbalanceChannel`: the forward map
  $y=I+jQe^{j\phi}$, optional AWGN, the real matrix $R(\phi)$, and the
  mean energy ratio.
- `demodulate_phase_imbalance(modulation, symbols, phi)` in
  `modcode.modulation.involution`: known-phase inverse.
- `demodulate_tau_aware(modulation, symbols, hypothesis=0|1)` in the same
  module: the two hypotheses; hypothesis $1$ conjugates the received
  symbols.
- `modcode.modulation.involution.tau_log_likelihood_ratio`: the per-symbol
  log-likelihood ratio for $\phi=0$ versus $\phi=\pi$.
- `block_tau_hypothesis`: sums the LLR over a block with one constant
  hypothesis, with an explicit numerical tie threshold.
- `pilot_tau_hypothesis`: selects the hypothesis from known pilot symbols.

### 6.2 Results at $\phi=\pi$

The noiseless run gives:

| Detector | BER | Why |
|:---|:---:|:---|
| standard coherent, assumes $\phi=0$ | $0.5000$ | Q bit inverted, I bit correct |
| oracle, knows $\phi=\pi$ | $0.0000$ | exact inverse |
| pilot-assisted | $0.0000$ | one known pilot resolves the sign |
| blind, tests $\phi=0$ vs $\pi$ | $0.5000$ | LLR identically zero; tie convention |

The finite-SNR run (40,000 bits per point) gives:

| Eb/N0 (dB) | coherent | oracle | pilot | blind | pilot hyp. error |
|---:|---:|---:|---:|---:|---:|
| 0 | 0.4975 | 0.0766 | 0.1201 | 0.4999 | 0.098 |
| 2 | 0.4989 | 0.0384 | 0.0557 | 0.5005 | 0.037 |
| 4 | 0.5003 | 0.0125 | 0.0184 | 0.4995 | 0.011 |
| 6 | 0.5000 | 0.0024 | 0.0030 | 0.4997 | 0.002 |
| 8 | 0.5001 | 0.0001 | 0.0010 | 0.5000 | 0.002 |

The pilot detector uses one pilot per 32 data symbols and decides the
hypothesis independently in each block. At low SNR a pilot decision can fail;
that block then contributes a full reflection error. The `pilot-hyp-err`
column in the run report is the fraction of blocks with the wrong hypothesis.
A global decision over all pilots would be more reliable for a phase that is
constant over the whole run, but it is less realistic for a block-varying
impairment; the per-block version is the conservative choice here.

### 6.3 What the figure shows

`examples/figures/phase_imbalance.png` has four panels:

1. the QPSK constellation at $\phi=\pi$, with arrows to the reflected points
   and the I/Q decision axes;
2. BER versus Eb/N0 for the four detectors;
3. $\det R(\phi)=\cos\phi$, marking the rank collapse at $\pi/2$ and the
   reflection at $\pi$;
4. the blind per-symbol LLR at 8 dB, which is zero for every symbol.

### 6.4 What the experiment does not show

- It does not show that a “τ-aware” detector can recover **blindly**. The
  blind detector is the control that fails.
- It does not show that $\phi=\pi$ is the worst impairment. It shows a
  reflection at $\pi$ and a rank collapse at $\pi/2$.
- It does not implement a blind I/Q imbalance estimator or a coded hypothesis
  test in the BER experiment; the CRC constraint is a unit test, and a
  production receiver would need a defined block/frame structure.
- It does not show any advantage over mature I/Q compensation. The point is
  the structural correspondence and the identifiability limit.

---

## 7. Literature positioning

The structural analogues are real, but they are analogues:

- **PT-symmetric Hamiltonians.** A PT-symmetric Hamiltonian can have a real
  spectrum when the symmetry is unbroken; broken PT symmetry gives
  complex-conjugate pairs. This is the closest spectral analogue of an
  involution whose unbroken sector gives a real spectrum. It is not the same
  as the DSIN phase flip, and it is not a proof of RH.
- **Kramers–Wannier duality.** The 2D Ising duality is an involution with a
  distinguished fixed point at the critical temperature. It shares the
  pattern “involution plus fixed locus,” not the phase-flip action.
- **The functional equation.** The completed zeta function satisfies
  $\xi(s)=\xi(1-s)$. A ratio of completed zeta functions can be arranged to
  satisfy $\varphi(s)\varphi(1-s)=1$. That is a functional-equation
  involution on a phase, not a relative phase between two quantum sectors.
  The two should not be identified without a construction.
- **I/Q imbalance.** I/Q gain and phase mismatch is a mature topic in
  communications. Standard references include Valkama, Renfors, and Koivunen
  (2001), Anttila, Valkama, and Renfors (2008), and Schenk (2008). Blind
  compensation algorithms exploit non-circularity, known constellation
  structure, or additional statistics; a circular QPSK/QAM signal with a
  uniform prior does not provide a second-order signature that distinguishes
  $\phi=0$ from $\phi=\pi$.

The contribution of this note is therefore modest:

1. a precise statement of what $\tau$ is, and what it is not;
2. the correct I/Q matrix and the $\phi=\pi$ identification;
3. the basis-dependence of the logical action;
4. the identifiability theorem that explains when a blind detector cannot
   work and what kind of side information can break the tie;
5. a runnable demonstration of all four cases.

It does not claim a new modulation scheme, a new blind estimator, a new
cryptographic protocol, or a connection to the zeros of zeta beyond the
structural analogy.

---

## 8. What this is not

- **Not a new modulation scheme.** I/Q imbalance is a hardware impairment;
  the experiment is a robustness study.
- **Not a security proof.** The DSIN discussion is about an invariant
  observable, not a composable security bound.
- **Not evidence for RH.** No statement here constrains the zeros of the zeta
  function.
- **Not an identity between the two models.** The real $2\times2$ matrices
  agree at $\phi=\pi$ after identifying sectors with I/Q axes. The complex
  structures and the encoding bases differ.
- **Not a blind recovery claim.** The blind detector fails by construction.
  The recovery uses an oracle or a pilot.

---

## 9. Reproduction

```bash
python -m pytest tests/test_phase_imbalance.py
python examples/phase_imbalance_lab.py
```

Files:

- `src/modcode/channel/phase_imbalance.py`
- `src/modcode/modulation/involution.py`
- `src/modcode/modulation/involution.py` (`demodulate_phase_imbalance`,
  `demodulate_tau_aware`)
- `src/modcode/viz/phase_imbalance.py`
- `examples/phase_imbalance_lab.py`
- `tests/test_phase_imbalance.py`

---

## References

1. C. M. Bender and S. Boettcher, “Real spectra in non-Hermitian Hamiltonians
   having PT symmetry,” *Physical Review Letters* **80** (1998), 5243–5246.
2. H. A. Kramers and G. H. Wannier, “Statistics of the two-dimensional
   ferromagnet. Part I,” *Physical Review* **60** (1941), 252–262.
3. M. Valkama, M. Renfors, and V. Koivunen, “Advanced methods for I/Q
   imbalance compensation in communication receivers,” *IEEE Transactions on
   Signal Processing* **49** (2001), no. 10, 2335–2344.
4. L. Anttila, M. Valkama, and M. Renfors, “Circularity-based I/Q imbalance
   compensation in wideband direct-conversion receivers,” *IEEE Transactions
   on Vehicular Technology* **57** (2008), no. 4, 2099–2113.
5. T. Schenk, *RF Imperfections in High-Rate Wireless Systems: Impact and
   Digital Compensation*, Springer, 2008.
6. The sibling `riemann-framework` checkout,
   `docs/phase_shift_involution.md`, for the DSIN-side statement and the
   cross-repository pointer.
