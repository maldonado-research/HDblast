# Which term sets the radiation contrast on the rolling B1 wall?

At the same wall proper time $s=H_0\tau=2.15$, the positive and negative rolling pulses change the fine-grid radiation density by **+3.67% and −3.35%** relative to the quiet run. The quiet wall already produces a massive-scalar cohort and its decay radiation. At this late cut the remote compensating disturbance has also had time to arrive. The comparison concerns the entire prepared perturbation on an already rolling wall. [Independent results](ENERGY.json)

## Crossing yield and dependence on the production rule

The B1 event prescription targets a dimensionless cohort density $N_*=q^{3/2}F/(8\pi^3)$, where $q=100|d\phi_b/ds|_*$ and $F=1+D/q$. The original event routine fits the correction $D$ to the preceding scalar and expansion histories. The quiet, positive and negative fine-grid events occur at $s=1.486122$, $1.485484$ and $1.486754$. The modeled radiation and massive gas both vanish for all three cases at the earlier common cut $s=1.45$. Percentages near the opening of the first cohort are ill-conditioned.

| Fine-grid change relative to quiet | Positive pulse | Negative pulse |
|---|---:|---:|
| Leading target yield from $q^{3/2}$ | −0.5190% | +0.7810% |
| Fitted multiplier $F$ | +4.3216% | −4.2503% |
| Target yield including $F$ | +3.7802% | −3.5025% |
| Recorded radiation at $s=2.15$ | +3.6677% | −3.3489% |
| Remaining massive gas energy at $s=2.15$ | +3.6326% | −3.4010% |
| Radiation with $F=1$ on the same wall trajectories | −0.6268% | +0.9415% |

The conditional $F=1$ comparison follows from linearity. For a prescribed single crossing and prescribed wall fields, the initial cohort is empty, the injection is proportional to $F$, the decay rate is independent of the occupation, and the radiation source is linear in that occupation. Consequently division of the two stored densities by each trajectory's $F$ removes this multiplier exactly in the fixed-background transport problem. Both late radiation contrasts reverse sign on both grids. Performing the same ablation directly in the independent transport solver gives −0.6270% and +0.9390% on the fine grid. A fully coupled change of prescription would also allow the metric and wall scalar to respond to the modified matter source.

Of the total positive-pulse change $\Delta(D/q)=+0.04235$, the term proportional to $\phi'''/\phi'$ contributes +0.03290, and the mixed $H\phi''/\phi'$ term contributes +0.00856. On the negative branch the corresponding contributions are −0.03285 and −0.00869 out of −0.04165. The positive event requires $G\ge216.1$ under the historical local-approximation screen. The model uses $G=100$. Quiet and negative events pass that particular screen. Matching the energy budget tests the numerical transport given the source prescription. It cannot validate the source's local expansion at the positive crossing.

## Mass-dependent decay and scalar work

Use densities in units of $H_0^4$, with $\widehat R=R\rho_b/b$, $X=\widehat\rho_\chi$, $h=H/H_0$, $v=d\phi_b/ds$, and $a=\exp(A_b-A_{\rm restart})$. The radiation balance is

$$
\frac{d(a^4\widehat R)}{ds}=a^4Q,\qquad Q=\frac{m_\chi^2}{8\pi H_0^2}\,\widehat n_\chi .
$$

The total matter budget in a unit comoving volume at restart is

$$
a^3(X+\widehat R) = \int a^3(j_\chi+j_{\rm prod})v\,ds
  -\int a^3 h(3p_\chi+\widehat R)\,ds ,
$$

starting before production. The first integral is work done by the rolling scalar. Expansion does pressure work on both components. Decay transfers energy between components and cancels in the total. Reconstructing $Q$ from the recorded number and the scalar-derived mass, then integrating it independently, matches all six stored late radiation densities to better than $10^{-8}$ relative. The largest total-energy residual is $4.2\times10^{-7}$ of accumulated scalar work. The quiet, positive and negative fine-grid scalar work integrals are 398.2863, 412.7977 and 384.7923. At this cut the respective physical radiation densities are 66.0643, 68.4873 and 63.8519. Remaining massive gas densities are 67.3641, 69.8112 and 65.0731. Radiation is present in the model, and a substantial massive component persists. The calculation neither measures thermalization nor establishes an origin for the quiet crossing.

## Independent reconstruction

The transport calculation takes only the sampled wall fields and their proper-time derivatives as trajectory input and uses the stored fitted $F$ for the corrected branch. It recovers $s_*$ from the zero of a cubic Hermite interpolant of $\phi_b-0.5$ and recovers $q$ from its derivative. It then evolves the cohort and radiation from empty initial data. The new implementation uses Gauss–Legendre quadrature in $z=\sqrt{\pi/q}\,k_*$ on $[0,8]$, with 48 independent momentum nodes, and a DOP853 integrator for the ramp, time-dilated decays and radiation. It uses no saved matter history in that transport solution. Increasing the momentum order from 32 to 48 changes the late solutions by less than $3.4\times10^{-8}$ relative. The reconstructed late component densities differ from the recorded ones by at most $5.9\times10^{-5}$ relative over both grids and three branches. These checks share the modeled masses, injection rule and decay law with the evolution. The audit uses the full raw matter histories separately to test energy conservation. [Reconstruction code](energy_budget.py)

With NumPy and SciPy installed, run from this directory using the included rolling histories.

~~~sh
python -I -B energy_budget.py ../classical/runs > regenerated-energy.json
~~~
