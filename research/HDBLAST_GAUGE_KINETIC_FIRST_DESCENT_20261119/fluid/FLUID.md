# How much conversion would make radiation dominant?

A source that draws energy from the rolling parent field can produce a radiation-dominated interval in either reference universe. Both very weak and very strong damping miss that target in this scan. For the assumed smooth conversion profile, the largest radiation shares found were 74.82% and 58.24%. These are conditional results for a relativistic fluid with pressure equal to one third of its energy density. The calculation specifies how an assumed energy-transfer rate would need to behave. It supplies no microscopic derivation of that rate, particle spectrum, thermalization time, or measurement of primordial radiation.

## The coupled model

The two launches start at rest on the same closed round geometries used in the preceding photon study, with zero radiation. Write $x=H_0t$, $h=H/H_0$, $v=dq/dx$, $a=H_0a_{\rm phys}$, and measure energies in units of $3M_{\rm Pl}^2H_0^2$. The parent density and pressure are $e_\phi=u+v^2/6$, $p_\phi=-u+v^2/6$. Let $r$ be the relativistic component's density. In this benchmark the dimensionless rate is $G=\Gamma/H_0=\gamma b(z)$, with

$$
s=\frac{q-q_i}{q_{\min}-q_i},\qquad z=\frac{s-0.2}{0.6},\qquad
b(z)=e^{4-1/[z(1-z)]}\quad(0<z<1).
$$

The window vanishes elsewhere and reaches one at $s=0.5$. It acts on every traverse of the same field interval, including a return before the first maximum radius. Positive $\gamma$ is an input. The model couples the dilution and conversion of energy to the evolving geometry,

$$
q'=v,\quad v'=-3u_{,q}-(3h+G)v,\quad Q=Gv^2/3,\quad
r'=-4hr+Q,\quad a'=ha,\quad h'=a^{-2}-v^2/2-2r.
$$

The constraint is $h^2+a^{-2}=e_\phi+r$. The sign and factor of three in $Q$ follow from $e_\phi'=-h v^2-Q$. We set $h_i=v_i=r_i=0$ and $a_i=1/\sqrt{u_i}$. Source-free radiation with these initial conditions stays exactly zero. The nonzero cases therefore test energy actually transferred from the parent in this model.

## Where the fraction crosses one half

Each trajectory is followed through the first maximum radius, $h=0$ after expansion. The “first descent” column restricts the same trajectory to its first arrival at the potential minimum. The table gives the two threshold rates detected by a 26-rate scan over $0\leq\gamma\leq352$ and continuous root-finding in each adjacent pair that straddled one half.

| Fitted branch | Half reached before first turnaround, lower and upper $\gamma$ | Half reached by first descent, lower and upper $\gamma$ | Peak efficiency near the best sampled rate |
|---|---:|---:|---:|
| $\phi_m=+0.5$ | 11.4521, 292.3114 | 11.5730, 285.7461 | 0.748188 at $\gamma\simeq46.73$ |
| $\phi_m=-0.5$ | 22.1759, 199.7142 | 23.2241, 194.3301 | 0.582433 at $\gamma\simeq62.28$ |

Every sampled rate between its detected full-interval thresholds reaches the target and every sampled rate outside misses it. The sweep does not establish global uniqueness or rule out additional behavior between tested rates or beyond 352. All 52 coarse-grid trajectories reach the first minimum and then the first turnaround within the calibrated field domain. The latest turnaround is at $x=11.203$. At $\gamma=352$, the first descents take $x=3.341$ and $2.111$, compared with $0.517$ and $0.541$ at zero coupling. The peak radiation shares fall to 0.445 and 0.399. Increasing the conversion coefficient also slows the field that supplies the energy. The source depends on its squared velocity and previously created radiation dilutes during the longer expansion.

On the positive branch at its local optimum the first half crossing occurs at $x=0.4423$, first arrival at the minimum at $0.8483$, and largest share at $0.9317$. On the negative branch the corresponding times are $0.6304$, $0.7755$, and $0.8096$. At first half, the integrated comoving injection divided by $a_i^4$ is 0.6642 and 0.8373, respectively. In both cases it equals the amount required by the constraint, $\tfrac12(a/a_i)^4(h^2+1/a^2)$, to numerical precision.

## Why the peak can follow the end of injection

The radiation fraction $F=r/(r+e_\phi)$ obeys

$$
F'=\frac{Q}{r+e_\phi}+h F(1-F)(3w_\phi-1),\qquad
w_\phi=\frac{p_\phi}{e_\phi}.
$$

Once the field leaves the conversion interval, radiation can keep gaining *share* while its physical density falls. The parent redshifts faster than radiation when $w_\phi>1/3$. At each of the four full-interval threshold peaks the source has already switched off. At those interior maxima the numerics give $w_\phi=1/3$, equivalently $v^2/6=2u$. This criterion distinguishes an increase in radiation's fraction from fresh energy injection.

We locate fraction extrema by bracketing zeros of its continuous derivative on subdivided adaptive integrator steps and including the endpoints. The primary integration evolves $\log(a/a_i),h,r$ and tracks $\int Q\,dx$, the expansion loss of each component, and $\int (a/a_i)^4Q\,dx$. Across the scan, the maximum absolute constraint or energy-budget error is below $1.6\times10^{-9}$. An independent tighter RK45 integration evolves $a,a',C=(a/a_i)^4r$ and the acceleration equation directly at both full-interval thresholds, both optima and the scan endpoints. The largest absolute change in the reported times, radii, fractions, densities or integrated inputs across these checks is below $1.3\times10^{-8}$ in model units. The zero-coupling first minima and first turnarounds agree with the archived scalar-only reference trajectories within $5.8\times10^{-10}$.

[Numerical data](FLUID.json) · [Reproducible calculation](fluid.py)
