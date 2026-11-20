# What do two measured response coefficients determine?

Two positive finite models can give the same response to a slowly changing field and emit different amounts of light when that field changes quickly. We construct such a pair on the two archived closed-universe first descents. The models have the same radius history, field history, incoming vacuum, and low harmonics. Only the coupling to higher harmonics changes. This inverse problem complements the calculation that varies the cutoff for a single prescribed local interaction.

## A finite positive oscillator model

Let $\ell=\log I$ be the archived smooth switch with endpoint change $\delta=0.1$, measured in conformal time $\eta$. For every transverse oscillator in shell $L$, choose the Hamiltonian

$$ H_L=\frac12e^{-2h_L\ell}\Pi_L^2+\frac12L^2e^{2h_L\ell}A_L^2 . $$

Every coefficient is positive throughout the switch. Both polarizations have the same real $h_L$, and all states within a shell have the same coupling. We take $0\le h_L^2\le1$ nonincreasing with $L$, set $h_L=1$ through $L=48$, and set $h_L=0$ above $L=192$. A coefficient that depends on harmonic number is a spatial form factor. These are different finite microscopic Hamiltonians. Their rotational symmetry and positive oscillator energies leave open whether a specified local five-dimensional action could realize their form factors. In particular the freely specified $h_L$ are additional model assumptions.

To first order about the incoming vacuum the real part of the anomalous covariance solves $R_L''+4L^2 R_L=h_L\ell''$. Each shell contributes to the linear force with one further factor $h_L$. Removing the common normalizations, the high-band susceptibility below the first possible pair frequency satisfies

$$ T(z)=\sum_{L=49}^{192}\frac{W_L h_L^2/L^2}{1-z/L^2},\qquad z=\omega^2/4,\qquad W_L=L(L^2-1). $$

The two matched quantities are $U_0=T(0)$ and $U_1=T'(0)$. These are the first two coefficients of the frequency-squared expansion. Higher coefficients determine further frequency dependence. The retarded function has the continuation $\omega\to\omega+i0$. The comparison on the real axis below threshold avoids the poles. The common low band can be added to either system.

Any monotone squared form factor in this class can be written $h_L^2=\sum_{k=48}^{192}p_k\mathbf{1}_{L\le k}$ with $p_k\ge0$ and $\sum p_k=1$. Here $p_k$ parameterizes the deterministic Hamiltonian. No random cutoff is sampled. For four illustrative support points, take the signed cofactors of a three-by-four rational matrix. Each cutoff contributes one column with entries $1$, $U_0$, and $U_1$. Normalize the positive and negative parts of that null vector separately to construct models A and B.

| Cutoff step, $k$ | 48 | 64 | 96 | 192 |
| --- | ---: | ---: | ---: | ---: |
| Model A, $p_k$ | 0 | 0.95071768 | 0 | 0.04928232 |
| Model B, $p_k$ | 0.50730013 | 0 | 0.49269987 | 0 |

Both have $U_0=1714.25659536$ and $U_1=0.338882482218$. The equality is exact before rounding. The next coefficient $\sum W_L h_L^2/L^6$ is $9.75882\times10^{-5}$ in A and $7.82438\times10^{-5}$ in B. We selected these support points during the investigation to illustrate nonuniqueness.

## What positivity guarantees below the band

The normalized susceptibility is an average of $f_z(t)=1/(1-zt)$ with positive weights $W_L h_L^2/(L^2 U_0)$, support $t=1/L^2\in[l,r]=[1/192^2,1/49^2]$, and known mean $m=U_1/U_0$. Convexity below the first pole gives

$$ f_z(m)\ \le\ T(z)/U_0\ \le\ \frac{r-m}{r-l}f_z(l)+\frac{m-l}{r-l}f_z(r),\qquad 0\le z<49^2. $$

The lower bound uses Jensen's inequality. The upper bound follows from the secant line through the two endpoints. At $\omega=20$ the interval is $[1.020167,1.020566]$. At $\omega=90$ it expands to $[1.66753,3.39179]$, and the two models give 2.18693 and 1.90182. The analytic enclosure uses the two positive moments and the support interval alone. Discrete support and monotonicity can tighten it.

## Energies after the coupling has settled

We evolved the full Gaussian covariance separately for each oscillator through each *prescribed* rolling history. The initial state is the exact conformal vacuum before the compact switch. The final constant coefficient supplies the energy measurement. With the unit-sphere degeneracies and natural units, $D=a_{\rm physical}^4\rho_\gamma=\pi^{-2}\sum W_L n_L$. Feedback on the trajectories is excluded from this comparison. The archived unit form factor is independently reproduced at a maximum absolute occupation difference below $6\times10^{-15}$.

| Field minimum | $D$, model A | $D$, model B | Ratio A/B for all shells | Ratio A/B above 48 |
| --- | ---: | ---: | ---: | ---: |
| $+0.5$ | 351.078 | 331.078 | 1.06041 | 1.48389 |
| $-0.5$ | 549.431 | 510.671 | 1.07590 | 1.32509 |

Thus the same two measured low-frequency coefficients permit at least these two distinct finite-amplitude outcomes. The discrepancy in energy above shell 48 is 48.4% and 32.5% relative to model B on the two branches. Normalized by the original parent density at first minimum, the total photon shares are 0.7110% versus 0.6705%, and 0.9967% versus 0.9264%. These examples leave the parent fixed during the comparison.

## Sharp Born bounds for the declared monotone class

At leading order in $\delta$, the energy after a prescribed pulse is linear in $h_L^2$. Its value at each sharp cutoff is the cumulative sum of the archived Born spectrum. All attainable energies at specified response moments therefore lie in the convex hull of the 145 sharp-cutoff columns. The minimum and maximum follow from a finite linear program in $p_k$. The program checks normalization, moment agreement, positivity and primal/dual agreement at the precision of the numerical spectrum. Each endpoint requires at most two or three sharp-cutoff steps for one or two fixed moments, respectively.

| Field minimum | Born $D$ with $U_0$ fixed | Born $D$ with both $U_0,U_1$ fixed |
| --- | ---: | ---: |
| $+0.5$ | 297.598 to 359.186 | 329.812 to 349.562 |
| $-0.5$ | 417.836 to 587.038 | 508.908 to 547.188 |

Here the constraints use the moments of the explicit pair. The solver considers **all** monotone profiles between the stated cutoffs. Both branches have A and B among the numerical Born maximizers and minimizers, respectively. The finite-amplitude energies in the preceding table include higher powers of the coupling and can fall outside these Born intervals. The construction establishes nonuniqueness under two calibrations for this finite model class. It supplies no bound on arbitrary local or nonlocal theories and no observational identification of the photon interaction. [Source](completion.py) · [Numerical results and checks](COMPLETION.json)
