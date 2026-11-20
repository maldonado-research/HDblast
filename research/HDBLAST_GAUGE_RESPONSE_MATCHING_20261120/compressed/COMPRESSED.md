# Can a few oscillators stand in for 144 photon shells?

Two rolling scalar backgrounds provide an independent test of an approximation suggested by the static probe. Keep every photon shell from $L=2$ to 48. Replace shells 49 through 192 by a smaller collection of oscillators whose frequencies and weights are determined entirely on a static sphere. Then let the scalar, radiation and geometry evolve together. This tests whether the static calibration predicts the energy transferred during a finite change in the electromagnetic coupling.

## Which quantities the approximation preserves

Write $w_L=L(L^2-1)$ for the dimensionless energy weight and $t_L=1/L^2$. In the static weak-field probe, the response of the upper band below its first pole is

$$
\chi(\omega)=\sum_{L=49}^{192}\frac{\mu_L}{1-(\omega^2/4)t_L},\qquad \mu_L=w_L/L^2=L-1/L.
$$

Gaussian quadrature with $m$ nodes for this discrete positive measure gives abscissae $t_j$ in the physical band and positive weights $b_j$. Its exactness for powers through degree $2m-1$ matches the first $2m$ coefficients of the Taylor series in $\omega^2/4$. The implementation uses a reorthogonalized Lanczos recurrence on $t_L/t_{49}$ and diagonalizes its tridiagonal matrix. Each effective oscillator has conformal frequency $L_j=1/\sqrt{t_j}$ and energy weight $W_j=b_j L_j^2$. For example, a single oscillator at $L_1=112.19745$ exactly reproduces the static inertia and the slope of the response in $\omega^2$.

For the dynamical test, the same weight $W_j$ multiplies the occupation in the energy equation and the anomalous covariance in the scalar force. With conformal time $\eta$, $s=d\log I/d\eta$ and complex covariance $M_j$, the solver evolves $n_j'=2s\,\mathrm{Re}\,M_j$ and $M_j'=-2iL_j M_j+s(2n_j+1)$. Vacuum initial data and $|M_j|^2=n_j(1+n_j)$ give nonnegative excitation energy. Photon work and scalar work enter with opposite signs. The effective oscillators approximate the combined upper-band energy. Their individual occupations have no identification with particular integer photon shells. The scalar kinetic coefficient is fixed to its bare value $Z=1$ throughout this test. [Static susceptibility derivation](../probe/PROBE.md) · [Evolution equations](../matched/matched.py)

## Full-model targets

The self-consistent calculations use the archived scalar potentials, vacuum photon data and smooth coupling change of log-amplitude $\delta$. The comparison event is the first passage through the potential minimum after the coupling has stopped changing. The physical low-band fraction and the total photon fraction use the same denominator, the combined scalar and photon energy density. The last column uses the total photon energy alone as its denominator.

| Scalar branch $\phi_m$ | $\delta$ | Photons as share of matter energy | Physical $L\le48$ as share of matter energy | Upper band as share of photons |
|---|---:|---:|---:|---:|
| $+0.5$ | 0.1 | 0.68495% | 0.53188% | 22.35% |
| $-0.5$ | 0.1 | 1.02039% | 0.62314% | 38.93% |
| $+0.5$ | 3 | 21.56277% | 8.28639% | 61.57% |
| $-0.5$ | 3 | 20.39901% | 5.17047% | 74.65% |

## What static matching predicts

The table gives the signed relative error in the total photon fraction, $100(F_m/F_{\mathrm{full}}-1)$, in percent. The frequencies and weights were selected before these four comparisons.

| Effective upper shells | Weak $+0.5$ | Weak $-0.5$ | Strong $+0.5$ | Strong $-0.5$ |
|---:|---:|---:|---:|---:|
| 1 | −13.4 | −14.6 | −29.7 | −26.7 |
| 2 | −5.92 | −2.99 | −28.4 | −19.6 |
| 4 | −0.974 | +1.84 | −11.1 | −19.5 |
| 8 | −0.331 | −0.0318 | +1.70 | −0.183 |
| 16 | +0.00449 | −0.00274 | −1.02 | +4.56 |
| 32 | −7.04e−9 | +2.96e−9 | +3.60 | −1.73 |
| 48 | −1.05e−10 | +5.11e−10 | +0.989 | +1.76 |
| 64 | −1.32e−10 | +1.88e−10 | +0.000243 | +0.00546 |
| 96 | −8.46e−11 | +2.43e−10 | +1.59e−8 | +2.94e−8 |

**Thirty-two effective shells reproduce both weak total photon fractions within $7.1\times10^{-11}$ relative error. The strong fractions still differ by several percent.** At this same order the 64 matched static moments agree within $4.7\times10^{-14}$ relative error and the static response at $\omega=96$, just below the first pole at 98, agrees to machine precision. Matching the low-frequency response alone gives no comparable accuracy guarantee for energy that remains in the field after a strong roll. The physical low-band fraction itself is more stable in the two strong cases, with relative errors of +0.116% and −0.113% at 32 nodes. The error in the added band's outgoing comoving energy is +5.79% and −2.28%.

Among the orders tested, 64 is the first to place all four total fractions within 1% relative error simultaneously. The largest error at that order is 0.00547% and at 96 nodes it is below $3.0\times10^{-8}$ percent. The error can change sign or grow before shrinking. On the negative branch the 8-node strong calculation happens to agree to 0.183%, while 16 nodes give a 4.56% error. Because a finite pulse samples a range of frequencies and the scalar responds to the emitted radiation, a test against the resolved trajectory matters even when many static Taylor coefficients match.

## Numerical checks and scope

The fully resolved $K=192$ integration reproduces the archived total photon fractions to $9.8\times10^{-11}$ absolute. At eight nodes, integrating the canonical field quadratures with tighter solver tolerances reproduces the covariance calculation's total photon fractions within $8.0\times10^{-10}$ absolute. Across the reduced runs, the maximal sampled work residual is below $3.7\times10^{-15}$ and the largest sampled constraint, scalar balance or covariance-purity residual is below $2.8\times10^{-10}$ in the solver's dimensionless units. Outgoing comoving photon energy is constant after the coupling stops. Entries close to the integration error floor should be interpreted as numerical agreement.

All weights and nodes come from an explicitly specified finite spectrum. The exercise quantifies the accuracy of a compressed representation for these two initial geometries and coupling histories. It neither determines the physical spectrum beyond $L=192$ nor extends the separate kinetic-counterterm model past its positive-inertia limit. Observational consequences of a higher-dimensional origin require additional physics linking such a source to the initial data used here.

From the study directory, run `OPENBLAS_NUM_THREADS=1 python -IB compressed/compressed.py --output compressed/COMPRESSED.json` with the archived parent inputs at `prior/`. The JSON includes the quadrature rules, the full-target and compressed endpoint values, the band energies at switch-off, the held-out static responses and the independent-coordinate checks. [Calculation](compressed.py) · [Numerical results](COMPRESSED.json)
