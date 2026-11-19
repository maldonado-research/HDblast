# Matching the static restoring force across cutoffs

Two finite models can have the same frozen static inverse susceptibility at every sampled parent momentum and still predict opposite changes in the shell's expansion. We imposed this equality and evolved the full five-dimensional shell separately with each cutoff. At proper time $s=2.10$, the fine quiet shell responds more slowly at pair cutoff 60 and more rapidly at 120. Both unperturbed and perturbed shells are expanding. The imposed static response is an exploratory benchmark chosen from earlier calculations. No measurement of a physical shell fixes its value. [Registered prediction](https://github.com/maldonado-research/HDblast/issues/157#issuecomment-7001899995) · [Sixteen new native evolutions](GRAVITY.json)

## What was matched

Set $a=1$ for the frozen response measurement, keep the free mass fixed, and write $\lambda=y^2$. The positive spectral bath shifts the inverse susceptibility at zero frequency by $-\lambda F^2 C_\Lambda(q)$, where the finite oscillator sum evaluates $C_\Lambda(q)=\sum_j g_{qj}^2/\omega_{qj}^2$ with the original 64-node rule. We denote the residual coefficient after adding a local quadratic interaction by $K(q)$. The matched inverse response is $q^2+m^2+\lambda F^2 K(q)$. This defines the static susceptibility through the order in $\lambda$ tested here, for the same parent normalization. On the evolving shell the additional Hamiltonian relative to the static-subtraction model is $\lambda a F^2 K(q)\chi_q^2/2$, with the same scalar clock on both cutoffs. [Finite Hamiltonian and spectral weight](https://doi.org/10.5281/zenodo.24358059)

The earlier choice $\kappa=1$ imposes $K(q)=0$. It cancels the frozen bath-induced change in the restoring force and constitutes one possible normalization. More generally, a measured inverse response could supply $K(q)$ after the free mass, gate and field normalization had also been specified. For the held-out test we set $K(q)=1.76 C_{120}(q)$. Thus $\kappa_{120}=2.76$ and $\kappa_{60}(q)=1+1.76 C_{120}(q)/C_{60}(q)$. The latter runs from 8.0400 to 8.5089 across the 24 momenta. Fitting a single constant at the upper central node $q=4.2562$ would miss the residual by between −4.11% and +2.27% elsewhere. The varying coefficient matches every native mode. It defines a momentum-dependent finite interaction in the specified comoving regulator. [Momentum-by-momentum values](GRAVITY.json)

## A prediction from the earlier expansion data

Let $S_\Lambda(K)=\partial_\lambda(H/H_0)$ at fixed $s=2.10$ on one grid. At zero coupling, the geometry, free parent, clock and initial data coincide across the two radiation cutoffs. The variation of their coupled equations is linear. An additional local $K(q)$ therefore induces the same parent, stress and geometry variation on both cutoffs. In particular,

$$
S_\Lambda(K)=S_\Lambda(0)+\mathcal L[K],\qquad
S_{60}(K)-S_{120}(K)=S_{60}(0)-S_{120}(0).
$$

Here $S_\Lambda(0)$ **uses the earlier static subtraction**. The bare model has $K=-C_\Lambda$. The linear functional $\mathcal L$ depends on grid, initial state, observable and gate, all fixed in this comparison. Earlier runs fix its value on two momentum profiles, $G_d=\mathcal L[C_d]=S_d^{\rm static}-S_d^{\rm bare}$ for $d=60,120$. Thus for $K=\alpha C_{120}+\beta C_{60}$ both cutoff outputs receive the same $\alpha G_{120}+\beta G_{60}$. Those two scalar measurements determine the functional on the span of the two profiles. The cancellation of a common response in the cutoff difference holds for any common $K$ with the same time-dependent Hamiltonian completion and no additional operators.

On the fine grid at step 0.005 the static slopes are approximately $-2.90\times10^{-5}$ and $-2.85\times10^{-5}$. The measured gains are $G_{60}=3.92598\times10^{-6}$ and $G_{120}=1.63821\times10^{-5}$, so either response can be computed directly for any combination of these two profiles. On the line $K=(\kappa_{120}-1)C_{120}$ the fine-grid zero crossings inferred from the old data are 2.77070 at cutoff 60 and 2.74173 at 120. Their coarse-grid counterparts are 2.78300 and 2.75321. The overlap of the opposite-sign intervals across both grids and both tested steps contains 2.76. We selected that value after examining the archived expansion slopes, then registered it before running the new native trajectories. The two baseline static models themselves both have negative slopes. [Earlier gravitational response](https://doi.org/10.5281/zenodo.24360096) · [Prediction protocol](https://github.com/maldonado-research/HDblast/issues/157#issuecomment-7001899995)

## Direct test and resolving power

The entries are absolute slopes $\partial_\lambda(H/H_0)$ in units of $10^{-7}$ per unit $\lambda$. Predictions use only the archived static and bare readouts. The observations use the metric velocity in new native runs with positive and negative source strengths. Each trajectory is interpolated to its own $s=2.10$.

| Spatial grid | Strength step | Predicted at 60 | Measured at 60 | Predicted at 120 | Measured at 120 |
|---|---:|---:|---:|---:|---:|
| Coarse | 0.010 | −3.7453 | −3.7460 | +1.1025 | +1.1018 |
| Coarse | 0.005 | −3.7434 | −3.7465 | +1.1049 | +1.1019 |
| Fine | 0.010 | −1.7277 | −1.7365 | +3.0118 | +3.0028 |
| Fine | 0.005 | −1.7534 | −1.7397 | +2.9927 | +3.0002 |

For the fine $\epsilon=0.005$ pair, the cutoff-60 minus cutoff-120 slope is $-4.740\times10^{-7}$, compared with $-4.746\times10^{-7}$ predicted from the earlier static data. At illustrative $y=0.1$ the individual fine effects correspond to about −0.00985 and +0.01698 parts per million of the baseline expansion. The energy in the outgoing pairs stays at its own cutoff's leading value, 21.57789 and 22.51057, in units of $y^2H_0^4$ per normalized comoving volume. These yields differ between cutoffs. Changing $K$ at a fixed cutoff affects the outgoing yield only at a higher order.

The independent source calculation propagates the 24 free covariances and their local-contact variations with DOP853. It derives their density, pressure, scalar current and all three proper-time derivatives directly from those covariances. The force depends on $F_\phi$ and its derivative includes $F_{\phi\phi}$. Against the independently propagated earlier RK4 difference at cutoff 120, all six quantities agree to $4.31\times10^{-8}$ or less relative to their respective peak magnitudes. The absolute source Ward residual is below $6.9\times10^{-12}$ over both grid histories. The native solver adds the new source at each stage and propagates the bulk and free parent from the original fork. Every forecast coefficient came from the preceding study.

Switching from cubic to linear matched-time interpolation changes the measured slopes by at most $1.24\times10^{-9}$. Halving the signed step changes them by at most $3.18\times10^{-10}$. The largest absolute disagreement with the predeclared forecasts is $1.38\times10^{-9}$. Moving from coarse to fine shifts the individual slopes by about $2.0\times10^{-7}$, comparable to their values near the tuned crossings. Both grids give the predicted pair of signs at both steps. A continuum claim requires further spatial refinement. A physical normalization, dynamical susceptibility and a covariantly matched total stress would be needed for a cutoff-independent cosmological inference. [All signed comparisons and source diagnostics](GRAVITY.json)

Run the complete test with the pinned upstream files present under the repository root and Python, NumPy and SciPy installed. The script verifies the prior inputs, reconstructs the four quiet fixed-wall sectors, independently integrates the matched contact and reruns all sixteen native trajectories.

~~~sh
OPENBLAS_NUM_THREADS=1 python -IB gravity.py --repo /path/to/repo --jobs 2 --output /tmp/GRAVITY.json
~~~

The optional `--work` directory retains the actual source tables and native trajectories. [Pinned solver and input manifest](gravity.py)
