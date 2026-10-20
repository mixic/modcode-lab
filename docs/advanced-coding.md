# Coding extensions: interleaving, fading, BCH/RS, and TCM

This note documents the coding and channel extensions that sit on top of the
core simulator in [`core-simulator.md`](core-simulator.md).  They are still
pedagogical implementations: the codes are short, the fading model is flat and
block-wise, and the TCM is a four-state example, not a standard.

---

## 1. Bit interleaving

`PermutationInterleaver(seed)` applies a deterministic permutation to a block
of encoded bits and inverts it before FEC decoding.  The permutation is
regenerated from the seed and block length, so it does not need to be
transmitted or stored.

```python
from modcode.pipeline import Link, PermutationInterleaver

link = Link(
    QPSK(),
    AWGNChannel(),
    Hamming74(),
    decoding="soft",
    interleaver=PermutationInterleaver(seed=9),
)
```

The interleaver is applied after FEC encoding and before modulation.  On the
receive side it is inverted before the hard or soft FEC decoder.  An
interleaver requires an FEC code; otherwise there is nothing to deinterleave.

This is a random permutation, not a standards-compliant block or convolutional
interleaver.  It is enough to spread one codeword's bits over different QAM
symbols.

### 1.1 Where the gain comes from

With QPSK, one symbol error flips two bits.  Without interleaving those two
bits can land in the same Hamming codeword, where a single-error-correcting
code cannot fix both.  A permutation separates them.

The AWGN panel of `examples/interleaver_lab.py` is a control: the gain is
small because symbol errors are independent.  The burst panel uses
`BurstChannel`, which erases contiguous runs of symbols.  The un-interleaved
Hamming curve floors around $1.3\times10^{-2}$, while the interleaved curve
continues down.  The burst channel is a modelling device, not a physical
fading model.

---

## 2. Block Rayleigh fading

`BlockRayleighChannel(block_length=10)` draws one circular complex Gaussian
coefficient per block, applies it to every symbol in the block, adds AWGN, and
zero-forcing equalizes using the known coefficient.

```python
from modcode import BlockRayleighChannel

channel = BlockRayleighChannel(block_length=10)
received = channel.transmit(symbols, noise_variance, rng)
coefficients = channel.last_coefficients
```

The receiver is assumed to know the channel state.  This is a **known-CSI**
model, not channel estimation.

For BPSK/QPSK the uncoded reference is

$$
P_b=\frac12\left(1-\sqrt{\frac{\bar\gamma}{1+\bar\gamma}}\right),
\qquad \bar\gamma=E_b/N_0,
$$

implemented as `theoretical_rayleigh_ber`.

After zero-forcing, the noise variance is $\sigma^2/|h|^2$ per symbol, which
is not a single scalar.  The current fading example therefore uses **hard**
decisions.  A soft fading path needs per-symbol variance supplied to the LLR
demodulator; that is a clearly marked extension, not a claim of the current
code.

`examples/fading_lab.py` compares the uncoded theory with BPSK/QPSK and then
shows Hamming (with and without interleaving), BCH, RS, and convolutional hard
decoding over block fading.

---

## 3. BCH(15,7)

`BCH15` is the primitive binary BCH code of length 15, designed distance 5,
correcting two bit errors.  The field is GF(16) with primitive polynomial
$x^4+x+1$; the generator is the product of the minimal polynomials of
$\alpha^1,\dots,\alpha^4$, of degree 8.

- **Encoding:** systematic polynomial division, message bits first.
- **Hard decoding:** syndromes at $\alpha^1,\dots,\alpha^4$,
  Berlekamp–Massey for the error-locator polynomial, Chien search for the
  positions, then flip the binary error positions.
- **Soft decoding:** maximum likelihood over the 128 codewords, using the same
  correlation metric as the Hamming decoder.

Tests cover noiseless round trips and **all** double-bit error patterns for a
codeword.  Soft decoding is practical only because the code is tiny; it is not
a scalable BCH soft decoder.

---

## 4. Reed–Solomon RS(15,11)

`ReedSolomon(m=4, t=2)` is a systematic RS(15,11) code over GF(16), correcting
two symbol errors.  Each symbol carries four bits.

- **Encoding:** generator with roots $\alpha^1,\dots,\alpha^4$, systematic
  polynomial division.
- **Hard decoding:** Berlekamp–Massey, Chien search, and a small GF linear
  solve for the error values.  The linear solve replaces a Forney formula and
  keeps the implementation readable.
- **Soft decoding:** not implemented.  `Link(decoding="auto")` falls back to
  hard decoding for this code.

Tests cover noiseless round trips and all double-symbol error patterns for a
block.  The code is a lab code, not a standards-compliant RS implementation.

---

## 5. TCM 8-PSK

`TCM8PSK` is a four-state, rate-2/3 trellis code mapped to Gray-labeled 8-PSK.
Two information bits enter each step; one parity bit is appended:

$$
p = b_1 \oplus b_2 \oplus s_1 \oplus s_2,
$$

where $(s_1,s_2)$ is the previous input pair.  The decoder runs soft Viterbi
on Euclidean distance.  Both TCM and uncoded QPSK carry two information bits
per symbol, so the comparison is at equal spectral efficiency.

```python
from modcode.modulation.tcm import TCM8PSK

tcm = TCM8PSK()
symbols = tcm.encode(bits)
bits_out = tcm.decode(received_symbols)
```

This is a pedagogical TCM, not a standardized Ungerboeck code.  It shows a
coding gain at moderate/high SNR but is worse than uncoded QPSK at very low
SNR; the crossing is visible in `examples/tcm_lab.py`.

---

## 6. Reproduce

```bash
python examples/interleaver_lab.py
python examples/fading_lab.py
python examples/bch_rs_lab.py
python examples/tcm_lab.py
```

Figures:

- `examples/figures/interleaver.png`
- `examples/figures/rayleigh.png`
- `examples/figures/bch_rs.png`
- `examples/figures/tcm.png`

Tests:

```bash
python -m pytest tests/test_interleaver.py tests/test_burst_channel.py \
  tests/test_rayleigh.py tests/test_bch_rs.py tests/test_tcm.py
```

---

## 7. Limits

- The interleaver is a random permutation, not a standard.
- The burst channel is synthetic.
- Rayleigh is block-flat and known-CSI; no estimation, Doppler, or multipath.
- BCH and RS are short pedagogical codes; the RS soft path is absent.
- TCM is four-state; no set-partition search or standard code.
- The fading example uses hard decisions because the equalized noise is not
  scalar.
