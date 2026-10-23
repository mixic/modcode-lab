# Scope and Limitations: What This Lab Measures and What It Cannot

### A review of the project's scientific limits, and what the repository does about each

*Scope statement for [modcode-lab](../README.md)*

---

## Abstract

An external review of this repository raised five limits: the pedagogical code
lengths, the simplified channel models, the status of the exploratory
appendices, the risk that documentation drifts from code, and the positioning of
the project as a whole. This note records each one, puts a number on it where
the repository can, states what would have to change for the criticism to lapse,
and points at the artifact that answers it — a test, a table, or a scope
banner.

The short version: `modcode-lab` is a reproducible lab, not a
production-grade modem and not a source of peer-reviewed novelty. That
positioning is deliberate, and this note exists so that it is *checkable*
rather than merely self-asserted. Every claim of the form "the documentation says
$X$" is pinned by a test in `tests/test_documented_parameters.py`, so the
documentation and the code cannot silently diverge.

---

## 1. Positioning

| Reader | Is this repository for them? |
|---|---|
| Student who wants to see why soft decoding beats hard decoding, how Viterbi works, or what LDPC and polar BER curves look like | **Yes** — that is exactly what it does, with theory overlays |
| Engineer who needs numbers comparable to 5G NR or DVB-S2X | **No** — see §2; the block lengths are one to two orders of magnitude too small |
| Researcher looking for state-of-the-art implementations | **No** — see §3 and §4; the models are teaching models |
| Reviewer who wants claims that are tested and reproducible | **Yes** — see §5; every documented parameter is pinned |

The lab is honest about being a lab. The purpose of this note is to make the
honesty operational: a limit that exists only in a disclaimer is a limit that
readers will miss, and the review is right that a reader should not have to
discover these limits by inspection.

---

## 2. Code scale: the toy parameters, with numbers

### 2.1 Where the codes sit

| Family | This repository | $N$ | $k$ | $R$ | Standards reference |
|---|---|---:|---:|---:|---|
| LDPC | `LDPC36` | 1008 | 506 | 0.502 | 5G NR: up to 8448 information bits per block; DVB-S2X: $N=64800$ |
| Polar | `PolarAWGN` | 1024 | 512 | 0.5 | 5G NR: mother code $N \le 1024$, CRC-aided list decoding |

The two rows have very different relationships to their standards:

- The **polar** block length is at the scale of the 5G NR mother code, but
  everything around it is absent: no rate matching, no 3GPP sequence, no
  distributed CRC, no interleaving. The comparison here isolates the decoder,
  not the waveform.
- The **LDPC** block length is far below both standards. It is $8.4$ times
  shorter than the 5G NR maximum information block and $64$ times shorter than
  a DVB-S2X normal frame. Its parity-check matrix is a regular $(3,6)$
  quasi-cyclic construction with a deterministic shift, not a density-evolved
  degree distribution.

### 2.2 The crossing point does not extrapolate

The repository measures a crossing between the LDPC and polar curves, and the
review is right that this number cannot be carried to a real system. The reason
is worth stating precisely, because the repository already contains the
evidence for it.

**Measured.** In [`polar.md`](polar.md) §4, one polar decoder at $1.5$ dB gives
BER $3.67\times10^{-2}$ at $N=256$ and $5.68\times10^{-3}$ at $N=2048$ — a
factor of $6.5$ from block length alone, with rate, decoder, and design point
fixed. At the same rate, block length is not a nuisance parameter; it is one of
the largest terms in the comparison.

**Consequence.** A crossing point measured at $N \approx 1000$ is a property of
*this* LDPC code paired with *this* polar code at *this* length. It is not a
statement that one family overtakes the other, and it is not evidence about 5G
NR, where the LDPC code has a different degree distribution, a layered
schedule, and a block up to $8.4$ times longer.

**What would have to change.** Density evolution or EXIT-chart design of the
degree distribution; a lifting size and base graph matching a standard; layered
normalized min-sum with offset correction; rate matching and CRC interleaving
for polar; and block lengths in the thousands to tens of thousands. Until then
the honest reading of every LDPC/polar comparison in this repository is: *for
these two short codes, at this length, with these decoders.*

**Status: limitation, measured.**

---

## 3. Channel model simplifications

Every channel here is a teaching model. The table states what each one assumes
and what breaks when the assumption fails.

| Model | Implemented | Assumes | Breaks when |
|---|---|---|---|
| AWGN | real-valued LLRs, exact noise variance | perfect synchronisation, no phase noise, no I/Q imbalance, no interference | any of those exist; then the LLR scale is wrong and the decoder is mis-modelled |
| Block Rayleigh | one coefficient per block, zero-forcing with the true coefficient | **perfect CSI at the receiver**; block-constant fading; hard decisions only | estimation is required, or the fades change inside a codeword, or soft LLRs need a per-symbol variance |
| Multipath/Doppler OFDM | comb pilots, **least-squares** estimation, **zero-forcing** equalisation, per-subcarrier LLR variance | pilot grid resolves the channel; no residual ISI after equalisation; no ICI; perfect timing and carrier recovery | the channel is fast relative to the pilot spacing, or the equaliser should be MMSE, or detection should be iterative |

Two of the gaps have measurements in this repository:

- **CSI is assumed, not estimated**, in the Rayleigh path. The model is
  `known-CSI` by construction; there is no channel estimator to be wrong. A real
  link spends pilots and SNR on estimation, and the fading result here is
  therefore an upper bound on performance, not a predictor.
- **LS with zero forcing costs 1.5–2 dB** relative to genie CSI in the
  2-tap multipath case at $20$ dB (see [`multipath.md`](multipath.md) §3), and
  the 4-tap case shows an interpolation floor. MMSE equalisation would recover
  part of that at low SNR, and iterative detection/equalisation would recover
  more; neither is implemented, and neither is claimed.

**Status: limitation; measured where a comparison exists.**

---

## 4. The exploratory appendices

The appendices — optical IM/DD, phased arrays and QFC, toy post-quantum codes,
number theory, CT-IM, the I/Q phase involution, lattice modulation, and the
sphere-packing note — are **worked examples with unit tests**. They are not
peer-reviewed contributions, and the README's "no paper-level novelty claim"
applies to all of them.

What "tested" means in this repository, precisely:

- the implementation matches the model described in the note;
- the invariants that the note states (round trips, distances, energy
  normalisation, closure properties) hold under test;
- the numbers quoted were produced by the script in the note.

What it does **not** mean:

- that the model matches a standard, a field trial, or a published measurement;
- that the result is new, or that it would survive peer review;
- that the choice of parameters was optimised for anything but readability.

Each appendix therefore carries a one-line scope banner pointing here. A reader
who wants scientific weight should look at the references, not at the plots.

**Status: definition.**

---

## 5. Documentation and code: traceability

This is the review's strongest point, and the one that can be closed rather
than merely acknowledged. A reader without the repository cannot tell whether
the documented constructions and parameters are the implemented ones.

The answer is not more prose; it is a test. Every parameter that a note quotes
is pinned by `tests/test_documented_parameters.py`, which fails if the code and
the documentation disagree. The mapping:

| Documented claim | Code symbol | Pinned by |
|---|---|---|
| LDPC block size $z=84$, 6 block rows, 12 block columns | `LDPC36.block_size`, `.block_rows`, `.block_columns` | `test_ldpc_block_structure` |
| $N=1008$, $M=504$, $k=506$, $R=0.502$ | `LDPC36.n`, `.m`, `.data_bits_per_block`, `.rate` | `test_ldpc_block_dimensions` |
| Every variable degree 3, every check degree 6, 3024 edges | `LDPC36.H` row and column sums | `test_ldpc_is_regular_three_six` |
| Normalized min-sum factor $0.8$, 50 iterations | `LDPC36.alpha`, `.max_iterations` | `test_ldpc_decoder_defaults` |
| Polar $N=1024$, $k=512$, CRC-8 of $x^8+x^2+x+1$ | `PolarAWGN.block_length`, `.data_bits`, `.crc` | `test_polar_block_and_crc` |
| Default list size $4$, Gaussian-approximation design at $2$ dB | `PolarAWGN.list_size`, `.design_ebn0_db` | `test_polar_decoder_defaults` |
| OFDM $N_{\text{fft}}=64$, pilot spacing $4$, least squares | `OFDMLink.n_fft`, `.pilot_spacing`, `.estimation` | `test_ofdm_defaults` |
| Multipath taps $(0.7,0.3)$, delays $(0,1)$, 16 sinusoids | `MultipathChannel` constructor defaults | `test_multipath_defaults` |
| Rayleigh block fading with known CSI | `BlockRayleighChannel.block_length` | `test_rayleigh_block_length` |

Run it directly:

```bash
python -m pytest tests/test_documented_parameters.py -q
```

The three checks that CI performs on every push are therefore also the
verification procedure for this note: `ruff check .`,
`python scripts/check_md_math.py docs README.md CONTRIBUTING.md`, and `pytest`.
A documentation change that contradicts the code now fails the build.

**Residual risk.** A pinned parameter can still be *described* wrongly in prose,
and no test can read the argument of a note. Pinning catches the numbers, the
names, and the defaults; it does not catch a subtle mischaracterisation. That
residue is why the notes keep a claims-audit table.

**Status: mitigated by test, not eliminated.**

---

## 6. What would have to change for benchmark use

Six concrete items, in the order that would yield the most credibility:

1. **Standards-scale codes.** DVB-S2X-class LDPC (64800) with density-evolved
   degree distributions, and 5G NR rate matching for polar.
2. **Realistic receivers.** Estimation rather than genie CSI, MMSE rather than
   zero forcing, and optionally iterative detection and equalisation.
3. **Impairments.** Phase noise, I/Q imbalance, carrier frequency offset,
   timing offset, and nonlinearity, each with its own measurement.
4. **Decoder parity with hardware.** Layered schedules, fixed-point behaviour,
   early stopping, and a throughput measure on a defined platform.
5. **Statistical rigour.** Confidence intervals on every BER point, error floors
   measured with enough bits, and seeds recorded per point.
6. **External reproduction.** A container, a pinned dependency set, and a
   one-command run that regenerates every figure in the notes.

None of these is present, and the repository does not claim them. A reader who
needs them should use a standard-compliant simulator.

---

## 7. Review risk register

| Review point | Severity | What this repository now does | Residual risk |
|---|---|---|---|
| Toy code parameters | high | scale table in §2, crossing-point reasoning, `ldpc.md` and `polar.md` compare against standards explicitly | a reader may still quote a curve out of context |
| Simplified channels | high | assumption table in §3, measured LS-versus-genie loss, explicit "known-CSI" wording | MMSE and iterative detection remain absent |
| Exploratory self-consistency | medium | scope banner on every appendix, §4 defines what "tested" means | a tested model can still be a weak model |
| Documentation/code drift | medium | every documented parameter pinned by test, traceability table in §5 | prose mischaracterisation is not automatable |
| Positioning | low | §1 states who it is and is not for; README repeats it above the fold | readers skim |

---

## 8. References

1. 3GPP, *TS 38.212: Multiplexing and channel coding*, Release 17. LDPC base
   graphs and lifting sizes; polar coding with CRC-aided list decoding.
2. ETSI, *EN 302 307-2: Second generation framing structure, channel coding and
   modulation systems for Broadcasting, Interactive Services, News Gathering
   and other broadband satellite applications (DVB-S2X)*, 2021.
3. T. J. Richardson and R. L. Urbanke, *Modern Coding Theory*, Cambridge
   University Press, 2008. Density evolution and EXIT design, the tools this
   repository does not use.
4. E. Arıkan, “Channel polarization: a method for constructing capacity-achieving
   codes for symmetric binary-input memoryless channels,” *IEEE Transactions on
   Information Theory* **55** (2009), 3051–3073.
5. I. Tal and A. Vardy, “List decoding of polar codes,” *IEEE Transactions on
   Information Theory* **61** (2015), 2213–2226.
6. C. Douillard, M. Jézéquel, C. Berrou, A. Picart, P. Didier, and A. Glavieux,
   “Iterative correction of intersymbol interference: turbo-equalization,”
   *European Transactions on Telecommunications* **6** (1995), 507–511.
7. M. Tüchler and J. Hagenauer, “Iterative decoding of serial concatenated
   convolutional codes over frequency-selective channels,” *IEEE Transactions on
   Communications* **49** (2001), 1966–1975.
