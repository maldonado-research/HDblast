# Can a vanishing current establish a thermal spectrum?

A charged gas can absorb energy from the HDblast wave and return it all before the wall finishes expanding. It can also lose 99% of its current while its distribution over momentum magnitudes retains the entire initial mismatch with thermal equilibrium. The two tests below follow particle distributions on the dense background of the published bounce. They measure different physical processes.

## Reversible acceleration on the actual field history

Take two equally populated massless fermion species of opposite charge and an initially isotropic Fermi–Dirac distribution with zero chemical potential. This is a **test population** in a classical kinetic continuum. Its initial temperature is an assumed parameter, and its stress and current do not feed back on the stored history. Only a small contribution to the background gravity and wave current would justify treating an actual species this way. A microscopic application on a compact sphere would also require thermally populated wavelengths short compared with its radius and a conversion between model units, ℏ and measured energies.

A single charge vector couples to one linear combination of the three Abelian fields. We call its coupling to the archived amplitude $q_{\rm eff}$. This includes the projection and common field normalization. In the homogeneous orthonormal frame of the round sphere, the electric and magnetic fields of this combination point along one fixed frame axis. The antisymmetric spatial connection gives no geodesic force proportional to the square of momentum. The magnetic force rotates the transverse momentum, which leaves an axisymmetric distribution invariant. In this frame define $T_c=a_iT_i$, comoving momentum magnitude $y=a|\mathbf p|/T_c$, and $\mu$ as the cosine to the chosen field axis. With the charge sign convention in the program,

$$
s(\tau)=\frac{q_{\rm eff}}{T_c}\int_0^\tau aX\,d\tau'
=\frac{q_{\rm eff}}{T_c}[A(\tau)-A(0)],\qquad A=a^2B_{\rm mag}/2.
$$

The collisionless Lorentz equation gives $f_\pm(y,\mu)=[1+\exp(\sqrt{y^2+s^2\mp2ys\mu})]^{-1}$. Comoving particle number and fine-grained Fermi entropy are invariant under the translation of momentum. For an initial dimensionless momentum $u$ and $v=|s|>0$, the mean final dimensionless energy over initial directions is

$$
M(u,v)=\begin{cases}u+v^2/(3u),&u\ge v,\\v+u^2/(3v),&u<v.\end{cases}
$$

Write $F_j(\eta)=\int_0^\infty u^j/[1+\exp(u-\eta)]\,du$. The normalized comoving energy $L=(a^4\rho)/(a^4\rho)_i$ is the integral of $u^2M(u,v)/(1+e^u)$ divided by $F_3(0)$. The charge current relative to $q_{\rm eff}n_{\rm total}$ follows from the same integral with $M$ replaced by its signed shift derivative and $F_3$ by $F_2$. Differentiating $L$ gives back the field work including the factor $F_2/F_3$.

For the illustrative $q_{\rm eff}/T_c=100$, the particles carry **1.68580 times their original comoving energy** at wall time 0.2 and **1.17964 times** at 8.1. The corresponding normalized currents are 0.80241 and 0.51985. The pulse-induced potential changes sign between 10.4805 and 10.481, with a linearly interpolated zero at 10.48052. At the zero, before the first expansion maximum, the original distribution in comoving momentum is exactly recovered in this collisionless model. The local kinetic energy has also redshifted with expansion. The same kinematic identity holds for the other charge/temperature ratios 1, 10 and 300 in the data. The exact return applies to this test on the published field and would have to be recomputed with appreciable backreaction.

Adding pairs introduces a further unknown distribution. In a dilute collisionless model an injection at time $u$ enters the solution as its birth distribution translated by the subsequent impulse $q_{\rm eff}[A(t)-A(u)]$. Summing these contributions needs the momentum-resolved source and its energy budget. A count of created pairs alone cannot determine the final spectrum. The kinetic runs here set the birth source to zero.

## A controlled relaxation experiment

Start at $\tau_*=8.1$ with the shifted distribution just computed. In this separate experiment turn off the force on this test gas after preparation and retain the expanding scale factor until the first turnaround at 12.42163. The original electromagnetic field continues to act in the collisionless experiment above. The following analytic experiment isolates the effects of collisions. It does not model those collisions on the continuously driven trajectory.

For each species let $\bar f$ be the angular average and $f_*$ the isotropic Fermi–Dirac shape with the same comoving number and energy. The species have opposite drifts, hence their total momentum vanishes. Elastic redistribution between the two species may relax that relative drift. In comoving momentum the model equation is

$$
\partial_\tau f=\nu_{\rm ang}(\bar f-f)+\nu_{\rm kin}(f_*-f).
$$

Both rates are nonnegative inputs. The first term randomizes direction at each momentum magnitude. The second also redistributes its magnitude. Both conserve comoving energy and the number of each species for these symmetric distributions. The target is time independent in comoving momentum. With $\chi_i=\int\nu_i d\tau$, the exact solution is a convex combination of the initial distribution, its angular average and the target, with weights $e^{-\chi_{\rm kin}}e^{-\chi_{\rm ang}}, e^{-\chi_{\rm kin}}(1-e^{-\chi_{\rm ang}}), 1-e^{-\chi_{\rm kin}}$. This preserves the Fermi occupation bounds. Concavity of the binary entropy makes angular averaging entropy increasing. At fixed number and energy, $f_*$ maximizes that entropy. These facts give nonnegative entropy production for the two operators.

We measure the distance of the radial **particle number spectrum** from the target by total variation, $D=(2F_2(0))^{-1}\int y^2|\bar f-f_*|\,dy$. The current decays by $e^{-\chi_{\rm ang}-\chi_{\rm kin}}$, whereas $D$ decays by $e^{-\chi_{\rm kin}}$. A common remaining current of 1% gives the following distinct outcomes.

| Integrated angular rate | Integrated kinetic rate | Final radial distance | Fermi entropy per particle in units of $k_B$ |
|---:|---:|---:|---:|
| $\ln100$ | 0 | 0.08406 | 4.69080 |
| $\ln10$ | $\ln10$ | 0.00841 | 4.71672 |
| 0 | $\ln100$ | 0.00084 | 4.71696 |

The preparation entropy is 4.20183 per particle and the matched equilibrium entropy is 4.71702. The total entropy can increase substantially through angular scattering alone. Its increase and a small current are individually insufficient tests of a thermal radial shape. Over the entire period of electromagnetic dominance the integrated Hubble rate is $\int H\,d\tau=0.12140056$. In the undriven relaxation law, setting $\nu_{\rm kin}=H$ throughout that duration leaves **88.57%** of the initial radial distortion. Removing 99% requires a logarithmic-expansion-weighted mean $\nu_{\rm kin}/H=37.93$. This duration comparison depends on when the distortion is prepared and on the relaxation law.

## A separate check for chemical equilibrium and temperature

Matching $n\propto T^3F_2(\mu/T)$ and $\rho\propto T^4F_3(\mu/T)$ at $\tau_*$ gives $aT/T_c=1.20276$ and $\mu/T=-0.59597$. Both charge signs here have the same abundance chemical potential because their separate numbers are conserved. For a massless initially zero-chemical-potential gas, conservation of pair number during energy gain gives $n/n_{\mu=0}(\rho)=L^{-3/4}$, independent of the subsequent elastic redistribution rates. At 8.1 the ratio is 0.88346. Attaining zero chemical potential at that energy would require **13.19% more pairs** relative to the available count. The fitted scale is a temperature parameter of the maximum-entropy target. A physical gas temperature requires demonstrating that the actual distribution approaches that target. Predicting blackbody photon radiation would additionally require a photon distribution, the exchange of energy with it, and suitable photon emission or absorption.

Run `python -I -B theory/spectrum.py --output /tmp/SPECTRUM.json` from the study root. The only dynamical input is `background/trace.csv.gz`. The program checks the work identity, number and energy using quadrature in initial and final momenta, angular resolution, the entropy interval and the force reversal. [Program](spectrum.py) · [Numerical results](SPECTRUM.json) · [Parent electromagnetic model](https://doi.org/10.5281/zenodo.24387565)
