# Convolution versus involution, and RM(1,3) as a tau-invariant code

*Exploratory appendix. This adds a code whose defining structure is a
coordinate involution, not polynomial convolution.*

The words sound alike, but they are different operations on different objects.
This note separates them and then implements the smallest useful example of an
involution-invariant code: the Reed–Muller code RM(1,3).

---

## 1. Convolution

Convolution combines two functions or sequences:

$$
(f*g)(t)=\int f(\tau)g(t-\tau)\,d\tau,
\qquad
(f*g)[n]=\sum_k f[k]g[n-k].
$$

A convolutional encoder is a convolution of the input bit sequence with the
generator polynomials.  At each step,

$$
c[n]=\sum_k g[k]u[n-k].
$$

Convolution is a binary operation.  It is associative, has an identity, and is
not generally its own inverse.  Its Fourier transform turns convolution into
multiplication.

In `modcode-lab`, convolution appears in:

- `ConvolutionalCode` — generator-polynomial convolution;
- Hamming, BCH, and Reed–Solomon encoding — polynomial multiplication/division;
- the convolutional trellis in `TCM8PSK`;
- periodic autocorrelation in `modcode.waveforms` (correlation with a reversed
  sequence).

---

## 2. Involution

An involution is a unary map that is its own inverse:

$$
\tau(\tau(x))=x
\qquad\text{or}\qquad
\tau^2=\mathrm{id}.
$$

It acts on one element and returns one element.  It does not combine two
objects.  Its iterates form a group of order two, $\{\mathrm{id},\tau\}$, not
a monoid of sequences.

In `modcode-lab`, involutions appear in:

- `phase-shift-involution.md`: $\tau=\mathrm{diag}(I,-I)$ on a graded
  Hilbert space;
- the dimension-shift note: $\sigma(s)=1-\overline{s}$;
- coordinate flips of a code, as below.

A cyclic shift is a permutation, but it is an involution only when its order
is two.  For length $n>2$, the cyclic shift generates $C_n$, not $C_2$.

---

## 3. The distinction

| Property | Convolution $*$ | Involution $\tau$ |
|---|---|---|
| Arity | 2 | 1 |
| Structure | monoid (associative, identity) | group of order 2 |
| Input | two functions/sequences | one element |
| Output | a new function/sequence | another element of the same set |
| Inverse | only for special elements | always, $\tau^2=\mathrm{id}$ |
| Fourier | diagonalized to multiplication | no general analogue |
| Role in coding | encoding operation | symmetry constraint on the code |

Convolution is a product; involution is a symmetry.  A convolutional code is
built by multiplying the input with a generator polynomial.  An
involution-invariant code is constrained so that applying $\tau$ to a
codeword returns a codeword.

---

## 4. RM(1,3)

The Reed–Muller code RM(1,3) has length $n=8$, dimension $k=4$, minimum
distance $d_{\min}=4$, and rate $1/2$.  A generator is

$$
G=
\begin{pmatrix}
1&1&1&1&1&1&1&1\\
0&0&0&0&1&1&1&1\\
0&0&1&1&0&0&1&1\\
0&1&0&1&0&1&0&1
\end{pmatrix}.
$$

The code is the span of these four rows.  Its 16 codewords have weights
$0,4,8$.

### 4.1 The involution

Let

$$
\tau(i)=i\oplus1
$$

on the coordinate set $\{0,\dots,7\}$.  This swaps
$0\leftrightarrow1$, $2\leftrightarrow3$, $4\leftrightarrow5$, and
$6\leftrightarrow7$.  It is fixed-point-free and $\tau^2=\mathrm{id}$.

RM(1,3) is invariant under $\tau$: for every codeword $c$, the permuted
word $c\circ\tau$ is also a codeword.  The generator rows make this easy to
check:

- the constant row is fixed;
- the $x_2$ and $x_1$ rows are fixed by adjacent-pair swaps;
- the $x_0$ row maps to its complement, which is the sum of the constant row
  and the $x_0$ row.

The implementation verifies this over the full 16-word codebook, not just the
generators.

### 4.2 It is not cyclic

Take $c=11110000$.  Its cyclic shift is $01111000$, which is not in
RM(1,3).  The test `test_rm13_is_not_cyclic` pins this.  So RM(1,3) is
invariant under a non-cyclic involution.

### 4.3 A correction to the tempting slogan

It is tempting to say “cyclic codes are a proper subclass of
$\tau$-invariant codes.”  That is **not** correct for a fixed $\tau$.
Cyclic codes are invariant under cyclic permutations; RM(1,3) is invariant
under a particular coordinate flip.  These are different permutation-invariance
classes.  A code can be cyclic, $\tau$-invariant, both, or neither.  What is
true is the broader statement: both are families of codes defined by a
permutation symmetry, and cyclic symmetry is one choice among many.

---

## 5. Encoding and decoding

`TauInvariantCode` stores a binary generator and a coordinate involution.  It
enumerates the full codebook (small codes only), verifies that every codeword's
$\tau$-image is in the codebook, and raises if not.

- **Hard decoding:** nearest codeword in Hamming distance.
- **Soft decoding:** maximum correlation over the codebook; for RM(1,3) the
  default path uses the Fast Hadamard Transform.
- **Rate:** $k/n$.
- **Minimum distance:** exhaustive over the codebook.

### 5.1 Fast Hadamard decoding

For RM(1,3), the 16 codewords are the eight Walsh characters and their
complements.  Correlating the received LLR vector with all of them can be done
with the Walsh–Hadamard transform:

$$
W = H_8 L,
\qquad
H_8 = H_2^{\otimes 3},
$$

where $L$ is the length-8 LLR vector.  The largest $|W_k|$ identifies the
linear function, and the sign identifies the constant term.  The implementation
builds a table from the codebook so the Walsh-frequency-to-codeword mapping is
explicit, then verifies that FHT and exhaustive ML return the same decisions.

The transform is $O(n\log n)$, or 24 butterfly additions for $n=8$, versus
$O(2^k n)=128$ correlation operations for exhaustive ML.  At $n=8$ in Python
the vectorized codebook ML and the batched FHT are within a few milliseconds of
each other because overhead dominates; the asymptotic advantage is the point,
not a measured speedup.  A production decoder would use a hardware or compiled
FHT.

`decode_soft_ml` exposes the exhaustive path for comparison, and
`decode_soft_fht` is the explicit transform decoder.  `decode_soft` dispatches
to the FHT when the codebook has the required affine structure.

---

## 6. Measured comparison with cyclic Hamming(7,4)

`examples/tau_invariant_lab.py` compares BPSK/AWGN performance at 100,000 bits
per point:

| Eb/N0 | Hamming hard | Hamming soft | RM hard | RM soft |
|---:|---:|---:|---:|---:|
| 0 dB | 1.20e-1 | 8.33e-2 | 1.77e-1 | 9.10e-2 |
| 2 dB | 5.28e-2 | 2.90e-2 | 9.27e-2 | 3.02e-2 |
| 4 dB | 1.69e-2 | 5.62e-3 | 3.41e-2 | 4.94e-3 |
| 6 dB | 2.22e-3 | 2.20e-4 | 6.12e-3 | 1.20e-4 |
| 8 dB | 8.00e-5 | 0 | 3.70e-4 | 0 |

The soft-decision curves are close; RM is slightly better at 4–6 dB in this
run, while hard-decision RM is worse at this rate because its rate and distance
distribution differ from Hamming.  The exhaustive nearest-codeword decoder is
already maximum-likelihood for the hard metric; the FHT reduces the soft
decoder's complexity, not its BER.  This is a finite-length comparison, not a
claim that one symmetry class is better.

The structural distinction is the point: Hamming is cyclic; RM(1,3) is
invariant under a non-cyclic involution.  Their AWGN gains are similar because
both are short one-error-correcting codes; the symmetry does not by itself
create coding gain.

---

## 7. Limits

- RM(1,3) is tiny; exhaustive decoding is not a scalable method.
- The FHT path is implemented for the affine RM structure; no general fast
  decoder for arbitrary $\tau$-invariant codes.
- No code design or search over $\tau$-invariant codes.
- No claim that involution-invariance is better than cyclic-invariance.
- No claim that this closes a “gap” in the FEC core; it adds a structural
  example.
- The soft-decoding result is finite-length, at fixed seeds.

---

## 8. Reproduce

```bash
python examples/tau_invariant_lab.py
python -m pytest tests/test_tau_invariant.py -q
```

Figure: `examples/figures/tau_invariant.png`.

---

## 9. References

- F. J. MacWilliams and N. J. A. Sloane, *The Theory of Error-Correcting
  Codes*, North-Holland, 1977.
- R. E. Blahut, *Algebraic Codes for Data Transmission*, Cambridge University
  Press, 2003.
