# Finite-sample pressure noise in one thermal snapshot

Suppose that a single charge stream is truly isotropic in its rest frame and that its representative momenta are independent classical Maxwell–Jüttner draws. How often would a snapshot slip outside the five-percentage-point pressure band by sampling noise alone? We use $T/m=2500$, close to the fitted temperature at radius doubling in a 512-particle run. We fit the same axial Landau rest frame to each artificial snapshot that is used in the particle diagnostic.

| Representatives present | Fitted-frame probability of missing $1/3\pm0.05$ | Known-rest probability in the ultrarelativistic limit |
|---:|---:|---:|
| 256 | 2.054% ± 0.037% | 1.992% |
| 450 | 0.1964% ± 0.0041% | 0.20405% |
| 512 | 0.0999% ± 0.0022% | 0.10054% |
| 1,024 | 0.0003440% ± 0.0000095% | 0.0003480% |

The uncertainties are one sampling standard error. The 450-particle example reflects a 512-representative run before all births have occurred. The live simulation continuously reuses particles, creates new ones and scatters pairs. Successive snapshots have time correlations. These probabilities calibrate an isolated thermal snapshot and cannot assign a probability to the measured duration of a live interval. [Calibration results](ANGULAR_NOISE.json) · [Coupled evolution](RESULTS.md)

## Expected variance

For a known rest frame, let $u_i$ be the axial direction cosine, $w_i=p_i^2/E_i$ and $R=\sum_iw_i u_i^2/\sum_i w_i$. Isotropy makes $u_i$ uniform on $[-1,1]$, with $\langle u_i^2\rangle=1/3$ and $\langle u_i^4\rangle=1/5$. Conditional on the weights,

$$
\operatorname{Var}(R\mid w)=\frac{4}{45N_{\rm eff}},\qquad
N_{\rm eff}=\frac{(\sum_iw_i)^2}{\sum_iw_i^2}.
$$

At high temperature $y_i=p_i/T$ follows a gamma distribution of shape three. The normalized weights then have a Dirichlet distribution with all parameters three. This gives the unconditional standard deviation $4/\sqrt{45(3N+1)}$, equal to 0.01521 at 512 particles and 0.01076 at 1,024. The reported tail probabilities use the full distribution. To evaluate the known-frame tail at threshold $q$, invert the characteristic function of $\sum_i y_i(u_i^2-q)$. A single draw has characteristic function $\int_0^1[1-it(u^2-q)]^{-3}du$.

For the finite-mass fitted-frame calculation, the sample axial stress transforms to $P_{zz,L}=[\sqrt{(e+P_{zz})^2-4S_z^2}-(e-P_{zz})]/2$. Transverse stress is invariant under the axial boost. To resolve the rare tails at 1,024 particles, the sampler tilts the massless distribution by $\exp[\lambda y(u^2-q)]$ and restores the exact massive thermal distribution with likelihood weights. Each tail uses 4,096 independent sample sets with a recorded seed. The finite-mass known-frame calculation independently agrees with the characteristic-function inversion within the reported sampling uncertainty.

## The bounded lab-energy diagnostic

An independent ultrarelativistic calibration covers the live bounded dispersion measure. In a known rest frame, for a fixed gas-energy fraction $f_g$ of the added wall energy, its massless limit is $f_g[1-|\sum_i y_i u_i|/\sum_i y_i]$. The same characteristic-function calculation with the first angular moment gives the probability of dropping below 90%.

| Representatives present | Gas supplies all added energy | Gas supplies 97.5% |
|---:|---:|---:|
| 256 | 1.6204% | 6.4650% |
| 512 | 0.06711% | 0.89556% |
| 1,024 | 0.0001505% | 0.02182% |

These examples fix the gas fraction and the physical frame across draws. Both can vary during actual evolution. The charge-reflected signs and the symmetry-related axes share the same representatives in the simulation and provide no additional independent samples.

Reproduce with `python -I -B live/angular_noise.py --output /tmp/ANGULAR_NOISE.json`. The calculation uses NumPy and SciPy.
