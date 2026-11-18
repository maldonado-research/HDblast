# Which part of the scalar’s work becomes radiation?

Consider one material volume expanding with the shell. The modeled scalar can do work in two ways. It creates particles at the energies they have at production, and it changes the rest masses of particles already present. Decay then moves energy from the particles to radiation. While the volume grows, both components do pressure work. Consequently, the radiation gained by the end of an interval can be substantially smaller than the energy that passed through the decay channel.

We checked the two archived $b=6.49\times10^{-6}$ runs at $\phi_*=0.5$, $G=100$, $y=1$, $\mu_0=0$, with 24 momentum bins, model detuning $\delta=0.1$ and seed offset 0.01. Use the forward clock $s=H_0\tau$, $A=\exp(A_b-A_b(0))$, $h=H/H_0$, $V=d\phi_b/ds$ and $r_b=H_0^{-1}$. The selected 1–4 interval starts with zero particle and radiation density on both grids. It lies inside the archived reliability interval. These are calculations on the model's stored trajectories. [Archived B1 model](https://github.com/maldonado-research/HDblast/blob/a011dee2b490faa56ae81d3165b4b5fa5bbfe636/research/HDBLAST_CHECKPOINT_20260930/B1_derived_chi_source/README.md)

## How the work is distributed

Integrate in $s$ with exactly one factor of $A^3$ throughout. The particle rest-mass work and filling work rates are $A^3b\widehat jV$ and $A^3b\widehat P_{\rm fill}$. The gross energy transfer to radiation is $D=\int A^3 b\widehat Q ds$. Compare this with the retained radiation gain $\Delta E_\gamma=\Delta(A^3r_bR)$ over the same interval. The combined physical-volume energy budget is
$$
W_{\rm mass}+W_{\rm fill}
=\Delta E_\gamma+\Delta(A^3b\widehat\rho_\chi)
+\int h A^3r_b R\,ds+3\int h A^3 b\widehat p_\chi\,ds .
$$
The integrals use all 153 original interior outputs on each grid and linear interpolation of derived columns at the two boundaries. Each recorded force was used in its original normalization. Integrating $A^3r_bJ_0V$ gives the same incoming work as the sum of its recorded components to floating-point precision. Integrating $A^3 b\widehat j\,d\phi$ and $A^3 b\widehat j_{\rm prod}\,d\phi$ instead uses actual field increments between saved outputs. That alternative quadrature shifts the total work down by about 0.030% on both grids. [Archive column definitions](https://github.com/maldonado-research/HDblast/blob/a011dee2b490faa56ae81d3165b4b5fa5bbfe636/research/HDBLAST_CHECKPOINT_20260930/B1_derived_chi_source/evolve_b1.py#L365-L395)

| On $1\le s\le4$ | Coarse | Fine |
| --- | ---: | ---: |
| Particle mass work | 0.0695106400 | 0.0695328637 |
| Particle filling work | 0.0459553423 | 0.0459700268 |
| Signed mass work / filling work | 1.51256930 | 1.51256957 |
| Gross decay energy / incoming work | 0.98567360 | 0.98567276 |
| Retained radiation gain / incoming work | 0.79089958 | 0.79090730 |
| Retained particle gain / incoming work | 0.00999303 | 0.00999520 |
| Radiation pressure work / incoming work | 0.19477418 | 0.19476567 |
| Particle pressure work / incoming work | 0.00430906 | 0.00430907 |
| Total gains and pressure work minus incoming work, ppm of incoming | −24.15 | −22.76 |

The two signed mass-work integrals include tiny negative contributions, $-2.89\times10^{-11}$ on the coarse grid and $-7.07\times10^{-13}$ on the fine grid. We integrated the positive and negative pieces of the linearly interpolated rates separately, splitting at any zero between nodes. The sign-reversed pieces vary strongly with resolution and carry no resolved physical return of energy. The fine result says that 60.2001% of the incoming work increased existing particle energies through the changing mass. This contribution is easy to miss if one charges the scalar only for creating quanta. The transfer fractions apply to this chosen interval and model. Determining a thermalization efficiency additionally requires particle interactions and distributions that predict a temperature.

## What would an instantaneous insertion cost?

The paid branch fills a cohort with the ramp $f(u)=u^3(10-15u+6u^2)$ up to $u=1$. Production consumes work from the moving field through $j_{\rm prod}V$. The implementation opens its cohort at the *end* of the solver step that crosses $\phi_*$. It gives the cohort an initial occupation fraction $f(u_0)$ with $u_0=|\phi_{\rm after}-\phi_*|/\Delta\phi_{\rm ramp}$. This discrete initialization has a small recorded positive energy jump. The integrated normal source then pays for ongoing production. [Cohort initialization and force](https://github.com/maldonado-research/HDblast/blob/a011dee2b490faa56ae81d3165b4b5fa5bbfe636/research/HDBLAST_CHECKPOINT_20260930/B1_derived_chi_source/shell_matter.py) · [Recorded crossing step](https://github.com/maldonado-research/HDblast/blob/a011dee2b490faa56ae81d3165b4b5fa5bbfe636/research/HDBLAST_CHECKPOINT_20260930/B1_derived_chi_source/evolve_b1.py#L302-L346)

For a matched counterfactual, freeze every field, momentum and crossing parameter at that same step and fill all the modes at once. Density depends linearly on occupation, so the cost is exactly the recorded jump divided by $f(u_0)$ for these initial zero-matter states, up to floating-point rounding. These are **local density jumps** in units $\kappa_5^2\rho/H_0$. Subsequent evolution of that alternative state has not been computed.

| Crossing-step quantity | Coarse | Fine |
| --- | ---: | ---: |
| $s_*$ | 1.48601507 | 1.48612162 |
| Initial $u_0$ | 0.00102358 | 0.000595981 |
| Initial occupation fraction $f(u_0)$ | 1.07077 × 10⁻⁸ | 2.11499 × 10⁻⁹ |
| Actual uncompensated step jump | 2.87626 × 10⁻¹² | 5.68282 × 10⁻¹³ |
| Hypothetical full occupation at the same step | 2.68617 × 10⁻⁴ | 2.68693 × 10⁻⁴ |
| Full occupation at the exact massless crossing, 24 nodes | 2.6880071 × 10⁻⁴ | 2.6880077 × 10⁻⁴ |
| Final stored $u$ | 1.00055576 | 1.00005496 |

Exactly at the crossing, integrating the massless continuum spectrum analytically gives $b(1+D/q)q^2/(4\pi^4)$, about $2.68754\times10^{-4}$ on both grids. This independently checks the scale of the full crossing cost. The continuum and 24-node numbers differ slightly through finite momentum quadrature. Cohort filling took a finite interval during which the particles were gaining mass. The paid creation work exceeds the full-occupation comoving jump evaluated exactly at the crossing by a factor of about 2.589 on each grid. The actual archived instant-control run used the much larger coupling $b=0.01$. Its recorded local energy jump, 0.41375, belongs to that different model cell. [Instant-control run summary](https://github.com/maldonado-research/HDblast/blob/a011dee2b490faa56ae81d3165b4b5fa5bbfe636/research/HDBLAST_CHECKPOINT_20260930/B1_derived_chi_source/runs/ctl/ctl_ps0.5_G100_y1_bneed_instant_dzf1e-3_summary.json)

The final numerical progress variable overshoots one on both grids. The code sets the filling derivative to zero once $u\ge1$. An independent test counts particles produced using the **number** continuity equation. For each occupied momentum, $m^2/\omega=\omega-p^2/\omega$. The modeled time-dilated decay rate gives
$$
N_{\rm births}=\Delta(A^3\widehat n)+\frac{y^2}{8\pi}
\int A^3(\widehat\rho_\chi-3\widehat p_\chi)\,ds .
$$
This checks the full cohort count without using $u$, the reported production force, or the energy transmitted to radiation. Comparing it with $A_*^3\widehat n_{\rm fin}$, inferred births are lower by 2.44 ppm on the coarse grid and 2.55 ppm on the fine grid. The endpoint overshoot in $u$ therefore gives a misleading estimate of the effect on occupancy. The number difference includes saved-cadence quadrature and evolution error.

The historical run classifier found no qualifying A1 radiation plateau, and the roughly 79% retained fraction measures the destination of the *selected scalar work*. It does not measure the radiation fraction of the total cosmic energy. The scalar work on the boundary alone does not close a five-dimensional global energy budget. [Original B1 data interpretation](https://github.com/maldonado-research/HDblast/blob/a011dee2b490faa56ae81d3165b4b5fa5bbfe636/research/HDBLAST_CHECKPOINT_20260930/B1_derived_chi_source/AUDIT.md)
