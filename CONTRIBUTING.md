# Contributing to modcode-lab

Thanks for considering a contribution.  This repository is a **pedagogical
reference implementation**, not a production modem or a research paper.  The
contribution contract is therefore stricter than “the code runs”: new material
must be tested, bounded in its claims, and placed in the core or the appendices
deliberately.

If you are unsure whether something belongs in the core, open an issue before
writing a large module.

---

## 1. Setup

Requires Python 3.10 or newer.

```bash
git clone https://github.com/mixic/modcode-lab.git
cd modcode-lab

python -m venv .venv
source .venv/bin/activate           # Windows: .venv\Scripts\activate
pip install -e ".[dev]"
```

Before sending a change:

```bash
ruff check .
pytest
```

The CI runs both on Python 3.10, 3.11, and 3.12.

---

## 2. The core/appendix contract

### Core

The core is the link-level signal chain:

- modulation and soft LLRs;
- FEC encoding and soft/hard decoding;
- AWGN, block fading, and burst channels;
- interleaving;
- BER and theoretical references;
- the `Link`/`run_link` composition.

A core addition must:

1. implement the relevant interface (`Modulation`, `FEC`, or `Channel`);
2. work with `Link` and `run_link`;
3. have focused tests, including failure paths;
4. have a document or a section in an existing document;
5. state what it does **not** claim;
6. avoid adding an appendix-only concept to a core interface.

The last point is important: `Modulation` does not carry phase-involution
methods, `Channel.transmit` receives a scalar noise variance rather than
coding parameters, and LDPC/RM are added through the FEC and appendix paths,
not by widening the core base classes.

### Appendix

Appendices are bounded experiments.  Examples: CT-IM, the I/Q phase
involution, lattice modulation, and RM(1,3).  An appendix must have:

- a scope box saying it is not core;
- a **tested claim**;
- a **bounded non-claim**;
- an **open experiment** or a clearly stated limitation;
- a runnable example;
- a figure indexed in `docs/figures.md`.

The template is the one used in
[`docs/exploratory-labs.md`](docs/exploratory-labs.md) and
[`docs/ct-im.md`](docs/ct-im.md):

> Here is what is claimed and tested.  Here is what is measured.  Here is what
> is not claimed.

Do not let an appendix drift into the core without the core-level tests,
interface fit, and documentation.  Do not let a core concept be described as a
result it has not earned.

---

## 3. Interfaces

### `Modulation`

Implement:

```python
bits_per_symbol: int
constellation: dict[tuple[int, ...], complex]
modulate(bits) -> complex symbols
demodulate(symbols) -> hard bits
```

`soft_demodulate(symbols, noise_variance)` has a generic max-log implementation
in the base class.  Override it only when the constellation has structure that
makes the generic N×M distance matrix unnecessary or numerically unsuitable.
The LLR convention is

```text
L = log P(bit=0 | y) / P(bit=1 | y)
```

positive for bit 0, and every FEC `decode_soft` must use that convention.

### `FEC`

Implement:

```python
data_bits_per_block: int
encoded_bits_per_block: int
rate: float
encode(bits) -> encoded bits
decode(bits) -> data bits
```

`decode_soft(llrs)` is optional.  If you implement it:

- positive LLR means bit 0;
- return exactly `data_bits_per_block` bits per block;
- raise `NotImplementedError` if the code cannot decode softly;
- test the sign convention, not only the noiseless round trip;
- state whether the decoder is ML, min-sum, Viterbi, or a heuristic.

`Link(..., decoding="auto")` uses soft decoding when both the modulation and
the FEC support it, and falls back to hard decoding otherwise.
`decoding="soft"` is strict and must fail loudly when soft decoding is
unavailable.

### `Channel`

Implement:

```python
noise_variance(ebn0_db, bits_per_symbol, code_rate) -> float
transmit(symbols, noise_variance, rng) -> received symbols
```

The channel receives a scalar variance and does not know the modulation order
or FEC rate.  The pipeline computes the variance once.  If your channel needs
per-symbol variance (for example after equalization in fading), do not silently
pass the scalar to a soft decoder; document the limitation or add a
channel-specific simulation path.

### Interleaver

`PermutationInterleaver` is optional in `Link`.  It requires an FEC code,
because there must be something to deinterleave.  The interleaver is applied
after encoding and before modulation, and inverted before decoding.

---

## 4. Tests

A focused test suite is the main review criterion.

### Documented parameters are pinned

If a note quotes a construction or a parameter — a block length, a degree, a
decoding factor, a pilot spacing — pin it in
`tests/test_documented_parameters.py` and add a row to the traceability table in
[`docs/scope-and-limitations.md`](docs/scope-and-limitations.md) §5.  A note is
allowed to describe a limited implementation; it is not allowed to describe code
that does not exist.

### One property per test

Name the property and assert it directly.  Prefer

```python
def test_hamming_corrects_every_single_bit_error():
    ...
```

over a single test that checks encoding, decoding, rates, and plotting at once.

### Pin numbers, not just trends

A test that only says “BER decreases with Eb/N0” cannot catch a wrong noise
variance or a wrong $E_b/N_0$ convention.  Pin at least one value against an
independent formula or a published reference.  Examples:

- `tests/test_ber_awgn.py::test_bpsk_ber_at_6db_matches_published_value`
  pins 6 dB BPSK to $2.388\times10^{-3}$;
- `tests/test_ber_awgn.py::test_theory_matches_simulation_at_high_snr`
  pins the 8 dB tail;
- `tests/test_tau_invariant.py::test_fht_matches_codebook_ml` pins exact
  decoder equivalence.

### Monte Carlo tests

Monte Carlo results are allowed only with:

- a fixed seed;
- an explicit sample size;
- a statistical tolerance justified by the sample size (for example a
  4-standard-error band), or a deterministic construction.

Do not write a test that passes only on a lucky seed.  Do not report a BER of
zero as an exact zero; use an upper confidence bound or a larger sample.

### Failure paths

Test invalid lengths, invalid parameters, non-binary inputs, wrong LLR signs,
and cases where an operation must raise.  A decoder that silently accepts a
malformed block is a bug.

### No hidden randomness

Use `np.random.default_rng(seed)` explicitly.  Do not rely on global state.
Do not reuse a generator across links in a way that makes one link's result
depend on another's presence; see the name-derived seed streams in
`metrics/ber.py`.

---

## 5. Style

- Python 3.10+; type hints on public functions.
- Ruff settings are in `pyproject.toml`: line length 88, rule families
  `E`, `F`, `I`, `UP`, `B`.
- NumPy for arrays; Matplotlib with the `Agg` backend in examples.
- No SciPy dependency; use `math.erfc` or a local implementation.
- Every Python file starts with the GPL-3.0-or-later SPDX header, after any
  shebang.
- Public modules should define `__all__`.
- Avoid private-attribute access across modules unless the module is an
  intentional internal companion (for example `metrics.lattice`).
- Documentation math uses GitHub-compatible delimiters: `$...$` inline and
  `$$...$$` for display.  Do not use `\(...\)` or `\[...\]`.

---

## 6. Documentation

A code change without a document is incomplete if it adds a claim.

- Core: update [`docs/core-simulator.md`](docs/core-simulator.md) or the
  relevant core note.
- Coding extensions: [`docs/advanced-coding.md`](docs/advanced-coding.md) or
  [`docs/ldpc.md`](docs/ldpc.md).
- Appendices: add the scope box, tested claim, non-claim, and open experiment
  to the appendix note; add a row to
  [`docs/exploratory-labs.md`](docs/exploratory-labs.md).
- Figures: generate them deterministically and add them to
  [`docs/figures.md`](docs/figures.md).  Do not add an unindexed figure to the
  README.
- Navigation: if the addition creates a new route, update
  [`docs/learning-paths.md`](docs/learning-paths.md).

Use precise language:

- “measured”, “tested”, “verified numerically” only when true;
- “should”, “may”, “expected” for unmeasured behaviour;
- “not claimed” for the boundary.

Do not describe a finite-length simulation as an asymptotic threshold.  Do not
describe a toy code as a standard.

### Math

Math is typeset on github.com, but the delimiters are parsed by CommonMark
first and the macros by KaTeX second.  Both stages can fail without an error in
the rendered page, so keep to these rules; `python scripts/check_md_math.py`
enforces them in CI.

- Inline math opens after a space, the start of the line, or `(`.  Write
  `base $b$`, not `base-$b$` or `"$a$`.
- Inline math closes before a space or punctuation, never before a letter,
  digit, or underscore.
- Display math puts `$$` alone at column 0, separated from the text by blank
  lines.  No body line may start with `-`, `+`, `*`, `1.`, `#`, or `>`, because
  CommonMark reads it as a new block and the math block collapses.  Keep a
  continuation operator at the end of the previous line.
- Prefer `\mathrm{ord}` to `\operatorname{ord}`: some renderers refuse
  `\operatorname` with “The following macros are not allowed”.

---

## 7. Adding a FEC code

1. Implement `FEC` in `src/modcode/fec/`.
2. Add `encode` and hard `decode`.
3. Add `decode_soft` only if it is genuinely soft; test the LLR sign.
4. Add a noiseless round trip.
5. Add a correction test (all single errors, or the code's guaranteed range).
6. Add a failure test for partial blocks.
7. Add a BER experiment if the code is meant to be compared.
8. Add a document with parameters, rate, distance, decoder, and limits.
9. Export it from `modcode.fec` and, if appropriate, `modcode`.
10. Run `ruff` and `pytest`.

If the code is an appendix (RM(1,3), for example), put the scope and non-claims
in its own document and keep the core FEC table honest about what is core.

---

## 8. Pull-request checklist

- [ ] The change is in the correct layer (core or appendix).
- [ ] Public interfaces are unchanged, or the breaking change is documented.
- [ ] Every Python file has the GPL header.
- [ ] New behavior has a focused test.
- [ ] At least one number is pinned independently, if the change makes a
      numerical claim.
- [ ] Monte Carlo tests have fixed seeds and a stated tolerance.
- [ ] `ruff check .` passes.
- [ ] `pytest` passes.
- [ ] Documentation states what is claimed and what is not.
- [ ] Figures are indexed.
- [ ] No appendix claim has leaked into the core README.

---

## 9. Scope boundaries

- Do not add a new module merely because it is interesting.  Add it because it
  has a tested claim, a bounded non-claim, and a user path.
- Do not claim a gain over a standard from a toy finite-length simulation.
- Do not claim security, cryptographic strength, or RH evidence from a toy.
- Do not use a Monte Carlo curve as proof of an asymptotic threshold.
- Do not silently change the noise normalization or the LLR sign convention.

The project is useful because it is honest about these boundaries.  Keep it
that way.

---

## 10. License

Contributions are accepted under GPL-3.0-or-later, matching the repository.
There is no separate CLA.
