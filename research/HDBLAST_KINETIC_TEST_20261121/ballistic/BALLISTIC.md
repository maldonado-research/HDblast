# First collisions on the accelerating wall

At the time the fluid calculation has received half of its startup absolute electrical work, can newly born particles still be moving freely in the electric field? We replayed the field, source and geometry from the fixed-pole fluid solutions. For each charge sign, a particle starts at its actual accepted birth time with a momentum selected from the directional decay law. The particle accelerates and its momentum redshifts until its first collision. We estimate that first-collision probability against a target population whose momenta also follow free trajectories. A separate upper bound on the pair collision rate covers arbitrary target momenta at the same density.

| Frozen field history | Unscattered at 50% of startup absolute work | At 90% of startup absolute work | At first hot crossing |
|---|---:|---:|---:|
| Dedicated vector pole | 34.5% | 8.1% | 7.0% |
| Shared vector, inclusive channel | 61.3% | 21.5% | 3.9% |

Each percentage is a number-weighted average over **all particles born by the indicated time**. Startup ends at the first hot crossing of the original fluid history. We mark time with the cumulative integral of $|D Q X v_{\rm fluid}|$, so its successive fractions reflect the fluid's electric work in both directions. The event times and the birth populations differ between the two histories. A comparison at their respective hot crossings therefore combines different field, density and birth histories. [Raw results and convergence](BALLISTIC.json)

## Ballistic momenta and collision probability

Energies and momenta below are in units of the daughter mass $m$. For a positive charge born at $t_b$, write $R_b=M_{\rm parent}(t_b)/(2m)$, $p_b=\sqrt{R_b^2-1}$, and let $\xi$ be its direction cosine along its axis. The accepted comoving rest increment $dD(t_b)$ sets its birth probability. Given $t_b$, the density in $\xi$ is proportional to $\min(1,|c-d\xi|/f)$, with the angular crossing speeds $c,d$ and requested flux $f$ taken from the pinned source. A constant acceptance gives isotropic directions. Let $I(t)=Q\int^t a X\,dt'$; then

$$
p_z(t)=\frac{a_b p_b\xi+I(t)-I(t_b)}{a(t)},\qquad
p_\perp(t)=\frac{a_b p_b\sqrt{1-\xi^2}}{a(t)}.
$$

A local magnetic field along the same axis rotates the transverse gyrophase. The independent uniform azimuth of the opposite charge removes the common reference direction from the pair average. The negative population is the momentum reflection of an independently sampled positive one. This approximation uses the locally homogeneous axisymmetric field and neglects connected two-particle correlations.

The invariant two-body flux supplies $v_{\rm M}=s\beta/(2E_1E_2)$ with $s=2(1+E_1 E_2-\mathbf p_1\cdot\mathbf p_2)$, $\beta=\sqrt{1-4/s}$. Set $M=8300$, $g_s=0.30(s+2)$. Multiplying the elastic pole by the flux gives the event kernel

$$
K=m^2v_{\rm M}\sigma=\frac{6\pi\beta g_s^2}
 {E_1E_2[(s-M^2)^2+\beta^2 g_s^2]}.
$$

The number density in the wall frame of the allowed opposite-sign target flavor satisfies $n_{\rm lab}L^3=S_0D/(6mL\,a^3)$. The hazard per unit dimensionless wall time is $\lambda_i(t)=b\,S_0D\langle K_i\rangle/[6(mL)^3a^3]$, with $b=1$ or $1/3$ for the dedicated or inclusive pole. Using lab density here accounts for particles occupying a common wall-frame volume. A tagged particle survives with probability $\exp[-\int_{t_b}^t\lambda_i(t')dt']$ on the prescribed target bath. The survival probabilities must be averaged individually. At the dedicated hot crossing the median integrated depth is about 59, but late births still contribute to the roughly 7% unscattered fraction.

At the half-work marker in the dedicated history, the unscattered tags carry approximately 48% of the kinetic energy those same cohorts would have if all particles continued freely. By the first hot crossing that fraction drops below 0.01%. These kinetic-energy fractions use a hypothetical ballistic denominator. Real scattered particles would acquire different momenta and contribute to both the field current and the target bath.

## How a collision deflects a particle

For the stated pure vector pole, the spin trace in the pair center-of-momentum frame has spatial factors proportional to $\delta_{ij}-\beta^2 n_i n_j$ for the incoming and outgoing directions. Contracting them gives $3-2\beta^2+\beta^4 z^2$, where $z=\cos\vartheta_*$ is the outgoing direction relative to the incoming axis. The azimuth is uniform. With $A=3-2\beta^2$ and $B=\beta^4$,

$$
p(z)=\frac{A+Bz^2}{2(A+B/3)},\qquad
\langle z\rangle=0,\qquad
\langle P_2(z)\rangle=\frac{2B}{15(A+B/3)}.
$$

This even law gives $\sigma_T=\sigma_{\rm total}$. Its second Legendre moment rises from zero at threshold to $1/10$ in the ultra-relativistic limit. The distribution of the outgoing charge in an individual event retains an angular quadrupole even when the ensemble-average transferred directed momentum equals the incoming momentum in the pair frame. Neither moment specifies the distribution of particle energies after many collisions.

## An upper bound valid for any target distribution

At fixed $s$, $E_1E_2\ge s/4$. A boost from the pair rest frame gives $E_1E_2=s/4+\gamma_b^2 v_b^2(E_*^2-p_*^2\cos^2\chi)$ with $E_*^2-p_*^2=1$. Thus $K\le2\beta m^2\sigma$. For $s\le M^2/2$, bound the squared denominator from below by $M^4/4$ and maximize $(s+2)^2/s$ at $s=M^2/2$. For $s\ge M^2/2$, put $x=s/M^2$, $\epsilon=.30$ and $d=\epsilon^2(1-8/M^2)$. The kernel is at most

$$
\frac{24\pi}{M^2}\epsilon^2(1+4/M^2)^2
\frac{x}{(x-1)^2+dx^2}.
$$

The rational factor peaks at $x=1/\sqrt{1+d}$ with value $(1+\sqrt{1+d})/(2d)$. Taking the larger bound gives $K<1.118569\times10^{-6}$. Integrating this maximum over the same density and birth histories gives an unscattered fraction of at least **0.206%** (dedicated) and **0.431%** (shared) at half work with cell-wise upper exposure and earliest births. A simple endpoint quadrature gives 0.210% and 0.434%, respectively. This lower bound still applies if targets redistribute among momenta and the same available target density, birth source, and specified collision channel are retained. Other interactions or a change in the field and density history require a new calculation. The transverse support of those unscattered particles also provides a direct test of the thermal distribution. [Transverse momentum test](TUBE.md)

## Resolution and scope

The dense run uses 512 tagged trajectories, 1,024 target draws at each time and 2,048 midpoint intervals. A pair of independent scrambled Sobol runs with 256 tags and 256 targets gave dedicated half-work survival of 33.4% and 33.9%, versus 34.5% in the dense run. Corresponding results for the shared case were 60.6%, 61.2% and 61.3%. At first hot crossing all three estimates were in 6.94–7.04% for the dedicated pole and 3.60–3.96% for the shared pole. Repeating the same intermediate sampling on the half-time-step field histories changed the six reported survival estimates by less than 0.02 percentage point. The archive includes the coarser 128-particle tests and individual energy and depth statistics. These convergence comparisons supply numerical sensitivity checks rather than formal confidence intervals.

The calculation tests free-flight memory and first encounters on an imposed environment. Once many particles have scattered, the freely moving target bath increasingly differs from a gas with collisions. It does not evolve the scattered target distribution or update the electric current and gravity from it. A first encounter may redirect a particle while leaving the full distribution far from thermal equilibrium.

Run `python -I -B ballistic/ballistic.py --traces /path/to/native --output /tmp/BALLISTIC.json --jobs 2` after generating the binary histories with `history/export.py`. The independent event-kernel controls include stationary charges, co-propagating beams, a head-on resonant pair, and partial-time exposure with electric translation and redshift.
