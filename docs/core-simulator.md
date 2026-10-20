# Core simulator: modulation, LLRs, FEC, and BER

This note is the technical reference for the **core** of `modcode-lab`:
the path from source bits through forward error correction, modulation, an
AWGN channel, soft or hard demodulation, FEC decoding, and bit-error
measurement.  The exploratory modules are deliberately not part of this
document; see [`exploratory-labs.md`](exploratory-labs.md).  The coding
extensions (interleaving, block Rayleigh, BCH/RS, TCM) are in
[`advanced-coding.md`](advanced-coding.md).

---

## 1. Signal and noise model

A complex baseband symbol is \(x=I+jQ\).  The AWGN channel adds circular
complex noise,

\[
y = x+n,
\qquad
n \sim \mathcal{CN}(0,2\sigma^2),
\]

so each real component has variance \(\sigma^2\).  The normalization used by
`AWGNChannel.noise_variance` is

\[
\sigma^2 = \frac{1}{2\,b\,R\,(E_b/N_0)},
\]

where \(b\) is the number of information bits per symbol and \(R\) is the FEC
code rate.  The factor \(2\) accounts for the two real dimensions.

---

## 2. Modulation

`modcode.modulation` provides a common interface:

```python
bits_per_symbol: int
constellation: dict[tuple[int, ...], complex]
modulate(bits) -> complex symbols
demodulate(symbols) -> hard bits
soft_demodulate(symbols, noise_variance) -> LLRs
```

The implemented constellations are BPSK, Gray-mapped QPSK, and square
Gray-mapped QAM (4, 16, 64, 256, ...).  Square QAM is sliced per axis: the
receiver scales each coordinate to the odd-integer level grid and rounds, so
hard slicing is O(N) and allocates no N-by-M distance matrix.  The QAM soft
LLR override uses the same per-axis decomposition.  APSK and the tail-biting
8-PSK example are exploratory; their constellation maps still work with the
generic LLR method, but they are not part of the core BER claims.

### 2.1 LLR convention

`soft_demodulate` returns

\[
L_k=\log\frac{P(b_k=0\mid y)}{P(b_k=1\mid y)}.
\]

Positive values favour bit 0.  This convention is used by every FEC
`decode_soft` method and by `run_link`.

The implementation uses the max-log approximation over the constellation:

\[
L_k \approx \frac{1}{2\sigma^2}
\left(
\min_{c:\,b_k=1}|y-c|^2
-
\min_{c:\,b_k=0}|y-c|^2
\right).
\]

For BPSK this reduces to \(L=2y/\sigma^2\).  For QPSK the two axes give
\(L_I=\sqrt2\,\Re(y)/\sigma^2\) and \(L_Q=\sqrt2\,\Im(y)/\sigma^2\).  For QAM
the max-log form is an approximation to the exact log-sum-exp LLR; it is the
standard practical choice and is adequate for the BER comparisons here.

---

## 3. Hamming(7,4)

`Hamming74` is the cyclic \([7,4,3]\) code with generator
\(g(x)=x^3+x+1\).  Data bits occupy the first four positions; the parity
checks are the three-bit remainder.

### 3.1 Hard decoding

Hard decoding computes the syndrome and flips the bit whose column matches
the syndrome.  It corrects any single bit error per codeword.  The 7 nonzero
syndromes map to distinct columns, so the implementation uses a 7-entry
syndrome-to-position lookup table and applies all corrections in one
vectorized step rather than looping over codewords.

### 3.2 Soft decoding

`Hamming74.decode_soft` performs maximum-likelihood decoding over the 16
codewords.  For codeword \(c\in\{0,1\}^7\), the correlation metric is

\[
M(c)=\sum_{i=0}^{6}(1-2c_i)L_i,
\]

and the decoder returns the data word of the codeword with the largest
\(M(c)\).  A positive LLR contributes positively to a 0 bit and negatively to
a 1 bit.  Soft decoding uses the reliability of every bit; hard decoding
discards it before the syndrome step.

The Hamming code is a teaching code, not a strong modern FEC.  It is useful
because the codebook is small enough to enumerate and because the hard/soft
difference is visible.

---

## 4. Convolutional code

`ConvolutionalCode` is a terminated rate-1/2 code with constraint length 3,
memory 2, and generators

\[
g_1=1+D+D^2 \quad(\text{octal }7),
\qquad
g_2=1+D^2 \quad(\text{octal }5).
\]

Each block of `block_bits` data bits is followed by two zero tail bits.  The
encoded block length is \(2(\text{block\_bits}+2)\), so the rate is slightly
below \(1/2\):

\[
R=\frac{\text{block\_bits}}{2(\text{block\_bits}+2)}.
\]

For the default `block_bits=64`, \(R=64/132\approx0.485\).

### 4.1 Trellis and Viterbi

The state is the last two input bits, so there are four states.  At each step
the encoder emits the two generator outputs.  The decoder runs Viterbi over
the trellis and forces the final state to zero, which is valid because the
tail bits are known.

For hard decisions the branch cost is the Hamming distance between the
received bit pair and the expected pair.  For soft decisions it is the
negative correlation between the expected bipolar pair and the two LLRs:

\[
\text{cost}(c)=-\sum_i (1-2c_i)L_i.
\]

Soft Viterbi is the main reason the core simulator is more than a Hamming
demo.  On the fixed test in `tests/test_convolutional.py`, at 2 dB the soft
decoder is roughly a factor of four lower in BER than the hard decoder.

---

## 5. Link pipeline

`Link(modulation, channel, fec=None, decoding="auto")` selects the decoding
path:

| `decoding` | Behaviour |
|---|---|
| `"auto"` | use soft decoding when the channel exposes a noise variance and the FEC implements `decode_soft`; otherwise fall back to hard |
| `"hard"` | force hard decisions followed by `fec.decode` |
| `"soft"` | require soft LLRs; raise if the channel or FEC cannot supply them |

`LinkResult.decoding` records which path actually ran.  This matters for
reproducibility: a BER curve is only comparable to another if the decoding
mode is stated.

`run_link` pads the source to whole FEC blocks, encodes, modulates, transmits,
decodes, and slices back to the requested number of bits.  The returned
`LinkResult` contains the source, received, transmitted, and symbol arrays.

The channel interface is intentionally narrow:

```python
noise_variance = channel.noise_variance(ebn0_db, bits_per_symbol, code_rate)
received = channel.transmit(symbols, noise_variance, rng)
```

The channel receives a scalar variance, not the modulation order or code rate.
`run_link` computes the variance once from the link parameters and passes the
same value to the channel and to the soft demodulator.  A custom channel that
does not expose `noise_variance` cannot be used by the pipeline; that is a
deliberate failure rather than a silent hard-decision fallback.

---

## 6. Theoretical BER references

`modcode.metrics.theory` supplies closed-form uncoded references:

- BPSK and Gray-coded QPSK:
  \[
  P_b=Q\!\left(\sqrt{2E_b/N_0}\right).
  \]
- Square Gray-coded \(M\)-QAM, nearest-neighbour approximation:
  \[
  P_b\approx
  \frac{4}{\log_2 M}
  \left(1-\frac{1}{\sqrt M}\right)
  Q\!\left(
  \sqrt{\frac{3\log_2 M}{M-1}\frac{E_b}{N_0}}
  \right).
  \]

The QAM expression is an approximation, not an exact bit-error probability.
`plot_ber(..., theory=...)` draws the references as dashed gray curves; the
measured curves use markers.  This makes a normalization error visible
immediately: if the simulation and theory do not overlap at high SNR, the
noise variance or the \(E_b/N_0\) convention is wrong.

### 6.1 Adaptive stopping and confidence intervals

`ber_sweep` keeps the historical fixed-`n_bits` behaviour by default.  For
serious curves use `ber_sweep_detailed`:

```python
from modcode.metrics import ber_sweep_detailed

curves = ber_sweep_detailed(
    {"BPSK": BPSK(), "QPSK": QPSK()},
    ebn0_db,
    min_errors=100,
    max_bits=1_000_000,
    batch_bits=50_000,
    seed=2026,
)
curve = curves["BPSK"]
print(curve.ber, curve.errors, curve.bits, curve.ci_low, curve.ci_high)
```

A point stops at `min_errors` errors or `max_bits` bits.  The interval is the
Wilson score interval at the requested confidence, which is well behaved at
zero errors: a zero-error point has `ber = 0` but `ci_high > 0`.  Plot the
upper bound, or report it explicitly, rather than treating a zero count as an
exact zero rate.

Each link gets a seed stream derived from its name, and each Eb/N0 point gets
an independent child stream.  Adding a modulation to a mapping therefore does
not change the results for the existing names, and reordering a sequence does
not silently share noise between links.

---

## 7. Reproduce the core results

```bash
python -m pytest tests/test_soft_decision.py tests/test_convolutional.py
python examples/run_ber_sweep.py
python examples/soft_decision_lab.py
```

The first script writes `examples/figures/ber_sweep.png`; the second writes
`examples/figures/soft_decision.png`.  Both use fixed seeds.  Monte Carlo
points are not pinned by tests; the codebooks, LLR signs, rates, and theory
values are.

---

## 8. What this core is not

- It is not a production modem.  The codes are small and the channel models
  are deliberately simple.
- It does not implement channel estimation, synchronization, equalization,
  fading, or RF nonidealities beyond the separate I/Q phase-imbalance
  experiment.
- It does not claim that Hamming or this convolutional code are competitive
  with modern LDPC/Polar/Turbo codes.
- The exploratory PQC/QFC/phased-array/CT-IM modules are not part of these
  BER claims; see [`exploratory-labs.md`](exploratory-labs.md).
