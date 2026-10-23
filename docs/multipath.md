# Multipath and Doppler: OFDM with pilot-based channel estimation

*Core P2 extension. This is a frequency-domain OFDM model with a pilot-based
least-squares estimator, not a 3GPP waveform.*

The core channels in [`core-simulator.md`](core-simulator.md) are AWGN and flat
block Rayleigh.  This note adds frequency selectivity (multipath) and time
variation (Doppler), and the receiver work multipath forces: channel
estimation, equalisation, and per-subcarrier LLR scaling.

Two objects are added:

- `MultipathChannel` — a symbol-spaced tap-delay line with unit average power
  and AR(1) Doppler.  It produces inter-symbol interference, so it is **not** a
  drop-in for `Link`: the scalar-noise LLR path ignores the ISI.
- `OFDMLink` — the decodable path: one flat gain per subcarrier, comb pilots,
  least-squares estimation, zero-forcing equalisation, and LLRs with the
  per-subcarrier post-equalisation variance.  This is the "channel-specific
  simulation path" the contribution contract asks for when a single noise
  variance is not valid.

---

## 1. Channel model

The discrete-time channel is

$$
y[k] = \sum_{\ell} h_{\ell}[k]\,x[k-d_{\ell}] + n[k],
$$

with integer delays $d_{\ell}$ and unit average power
$\sum_{\ell}\mathbb{E}|h_{\ell}|^2 = 1$.  A new tap vector is drawn per
simulation block; consecutive blocks are correlated by

$$
h_{\ell}[k] = \rho\,h_{\ell}[k-1] + \sqrt{1-\rho^2}\,w_{\ell}[k],
\qquad w_{\ell}\sim\mathcal{CN}(0,p_{\ell}),
$$

where the one-step correlation $\rho$ is the Jakes coefficient evaluated on a
fixed sum-of-sinusoids grid,

$$
\rho = \frac{1}{N}\sum_{n=1}^{N}\cos\!\left(2\pi f_d T\cos\theta_n\right),
\qquad \theta_n = \frac{2\pi n}{N+1},
$$

with $T$ the update interval.  $\rho=1$ is a static channel and $\rho=0$ gives
independent blocks.  The AR(1) model matches the Jakes correlation at one step
and keeps the average tap power exact; it is not the full Jakes spectrum.

The OFDM link uses the frequency-domain equivalent.  With a cyclic prefix at
least as long as the delay spread, and taps constant within a block,

$$
y_m = H_m x_m + n_m,
\qquad
H_m = \sum_{\ell} h_{\ell} e^{-j2\pi m d_{\ell}/N},
$$

so each of the $N=64$ subcarriers sees a flat gain.  The CP itself is not
simulated; intra-block Doppler (ICI) is therefore not modelled, and Doppler
appears as block-to-block tap variation.

---

## 2. Receiver

1. **Pilots.** Comb pilots every $s$ subcarriers, transmitted as $+1$; the other
   $N-sN_p$ subcarriers carry data.  The default is $N=64$, $s=4$, i.e. 16
   pilots, 48 data subcarriers, and an effective rate of $0.75$ before FEC.
2. **Estimation.** Least squares per pilot, $\hat H_p = y_p / x_p$, then
   periodic linear interpolation over the band.  With `estimation="genie"` the
   true gains are used instead, which separates estimation loss from the
   fading penalty.
3. **Equalisation.** Zero forcing, $z_m = y_m/\hat H_m$, with post-equalisation
   variance $\sigma^2/|\hat H_m|^2$.  Zero forcing is chosen because its
   variance is exact and per-subcarrier; an MMSE equaliser would need a bias
   term and is left out.
4. **LLRs and FEC.** The modulation's max-log LLRs are computed at unit
   variance and divided by the per-subcarrier variance, then passed to
   `FEC.decode_soft`.  This is why the path cannot be `Link`.

$\mathrm{Eb}/N_0$ is per information bit and includes the pilot overhead: the
noise variance is
$1/(2\,b\,R_{\text{eff}}\,\mathrm{Eb}/N_0)$ with
$R_{\text{eff}} = R_{\text{FEC}}\cdot N_{\text{data}}/N$.

---

## 3. Measured results

`examples/multipath_lab.py` sweeps uncoded QPSK over AWGN and block fading,
300,000 bits per point with fixed seeds.  QPSK is used so that the curves show
the channel, not a code; the link accepts an FEC code as well.

| Eb/N0 | flat, genie | 2-tap, genie | 2-tap, LS | 4-tap, LS |
|---:|---:|---:|---:|---:|
| 0 dB | 1.72e-1 | 1.72e-1 | 2.37e-1 | 2.44e-1 |
| 5 dB | 8.22e-2 | 8.24e-2 | 1.21e-1 | 1.26e-1 |
| 10 dB | 3.00e-2 | 3.00e-2 | 4.76e-2 | 5.39e-2 |
| 15 dB | 8.86e-3 | 1.01e-2 | 1.61e-2 | 2.37e-2 |
| 20 dB | 4.17e-3 | 3.15e-3 | 5.16e-3 | 1.20e-2 |
| 25 dB | 9.67e-4 | 1.03e-3 | 1.61e-3 | 8.77e-3 |
| 30 dB | 2.97e-4 | 3.50e-4 | 6.43e-4 | 7.64e-3 |

The flat and 2-tap **genie** curves agree, as they should: summing Gaussian taps
with power one gives the same Rayleigh marginal per subcarrier, and uncoded QPSK
does not exploit frequency diversity.  The flat genie points sit on the pilot-
adjusted flat-Rayleigh closed form (dashed curve in the figure), which validates
the noise and Eb/N0 accounting.

**Estimation loss.**  At 20 dB the 2-tap LS curve is at $5.16\times10^{-3}$
against $3.34\times10^{-3}$ with genie CSI — about 1.5–2 dB of estimation loss
from 16 pilots and linear interpolation.  Denser pilots lower the estimation
error but raise the noise per data bit at fixed Eb/N0:

| pilot spacing | effective rate | NMSE | BER at 20 dB |
|---:|---:|---:|---:|
| 2 | 0.500 | 7.37e-3 | 7.30e-3 |
| 4 | 0.750 | 4.53e-3 | 5.16e-3 |
| 8 | 0.875 | 4.60e-3 | 4.73e-3 |
| genie | 0.750 | 0 | 3.34e-3 |

The interpolation error is not monotone in spacing because the channel is
smooth at these delays; that is a property of this two-tap profile, not a
general rule.

**Interpolation floor.**  The 4-tap profile has delays up to five symbols.  With
pilot spacing four, linear interpolation leaves an error that dominates at high
SNR: the 4-tap LS curve flattens near $8\times10^{-3}$ instead of following the
genie slope.  This is an interpolation-limited floor, not a code floor.  A
denser grid or a DFT-based estimator would move it; neither is implemented.

**Doppler.**  The AR(1) one-step correlation across one 64-symbol block:

| $f_d$ (cycles/symbol) | 0 | 0.002 | 0.005 | 0.020 | 1.000 |
|---|---:|---:|---:|---:|---:|
| $\rho$ | 1.000 | 0.854 | 0.258 | 0.183 | 0.000 |

The OFDM link draws one tap vector per block, so $\rho$ sets how fast the
realisation changes; the BER curves use $\rho\approx0$ (independent blocks), so
the points average over many realisations.  Static channels
($\rho=1$) are a single realisation and must be averaged over drops if used for
BER; the figure therefore shows the fast-fading case.

---

## 4. Complexity and limits

- The estimator is $O(N)$ per block (one division per pilot, one interpolation
  pass); the equaliser is $O(N)$.
- Equalisation is zero forcing, not MMSE, and detection is not iterative: no
  turbo equalisation, no soft interference cancellation.  At $20$ dB the
  LS-plus-ZF receiver is 1.5–2 dB behind genie CSI, and the gap widens at
  low SNR, where zero forcing amplifies noise.  See
  [`scope-and-limitations.md`](scope-and-limitations.md) §3.
- The model is frequency-domain: no CP, no time-domain equaliser, no ICI, no
  timing or carrier synchronisation.
- Taps are drawn per block; there is no continuous-time Jakes process, only its
  one-step correlation.
- The 4-tap floor is an interpolation artefact, not a statement about 4-tap
  channels in general.
- Measured BER used fixed seeds; the tables are reproducible, and the points
  are not pinned by tests.

---

## 5. What is claimed, measured, and not claimed

**Tested**

- tap power is one on average and the static process has constant taps;
- the delay line places an impulse at the expected taps;
- the noise variance matches the requested moment;
- the subcarrier gains equal the DFT of the tap-delay line;
- interpolation is exact for a flat channel;
- a noiseless flat link decodes with genie and with LS estimation;
- an FEC-coded multipath link decodes below $10^{-2}$ at 20 dB;
- LS estimation costs SNR relative to genie CSI on a fixed sample.

**Measured**

- the BER table, the estimation loss, the pilot-density trade-off, the
  interpolation floor, and the Doppler correlation table.

**Not claimed**

- a 3GPP or Wi-Fi waveform, or any standard's numerology;
- optimal pilots, MMSE equalisation, iterative estimation, or ICI mitigation;
- a diversity or coding gain from multipath: uncoded flat and 2-tap curves
  coincide, and coded gains are a separate experiment;
- a hardware or latency claim.

---

## 6. Reproduce

```bash
python examples/multipath_lab.py
python -m pytest tests/test_multipath.py -q
```

Figure: `examples/figures/multipath.png` (Figure 7 in [`figures.md`](figures.md)).

---

## 7. References

- T. S. Rappaport, *Wireless Communications: Principles and Practice*,
  2nd edition.
- J. G. Proakis and M. Salehi, *Digital Communications*, 5th edition,
  chapters on fading and equalisation.
- W. C. Jakes, *Microwave Mobile Communications*, 1974.
- S. Coleri, M. Ergen, A. Puri, and A. Bahai, “Channel estimation techniques
  based on pilot arrangement in OFDM systems,” *IEEE Transactions on
  Broadcasting* **48** (2002), 223–229.
