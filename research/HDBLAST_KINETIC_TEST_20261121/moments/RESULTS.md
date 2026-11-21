# Can the wall's fluid measurements determine the collision force?

Two streams can have identical measured particle density, current, momentum, energy and pressure tensor and still lose directed momentum at different rates. At an early point on the preceding HDblast trajectory, two angular distributions with the same fluid measurements give rates differing by **19.9%**. With rest-frame rotational symmetry imposed and the exact equilibrium fluid moments used as targets, two different radial spectra give rates differing by a factor of **3.96**. These are local force calculations for the same elastic interaction. Following either spectrum as the wall expands would require its own evolution. [Numerical results](MOMENTS.json) · [Preceding HDblast calculation](https://doi.org/10.5281/zenodo.24408924)

## What the fluid measures

Measure individual energies $E$ and momenta in units of the common rest mass $m$. In a stream's rest frame let $n$ be its number density, $e$ its energy per particle, $\Theta=P/(nm)$ its pressure per particle, and $h=e+\Theta$. For any rotationally isotropic distribution, or any distribution with isotropic second angular moments on each energy shell, the relations are

$$
e=\langle E\rangle,\qquad 3\Theta=\langle p^2/E\rangle
=\langle E\rangle-\langle 1/E\rangle,\qquad
(h-\Theta)(h-4\Theta)=\langle E\rangle\langle 1/E\rangle\ge1 .
$$

Cauchy's inequality gives equality precisely when all particles have the same energy. The Taub–Mathews equation used in the earlier wall example saturates this inequality, with $E_*=3\Theta/2+\sqrt{1+9\Theta^2/4}$. An **exact** kinetic realization of these particular rest moments must lie on that one energy shell. A fluid approximation may of course represent a broader gas to finite accuracy. The distribution over directions on a shell can still vary. If two angular laws have $\langle n_i\rangle=0$ and $\langle n_i n_j\rangle=\delta_{ij}/3$, they share the full number current and stress-energy tensor. Applying identical opposite Lorentz boosts to the two streams preserves the equality between candidate models in the wall frame.

Three sets of angular laws illustrate the remaining freedom.

- One law puts equal weight on the cones $\cos\vartheta=\pm1/\sqrt3$, uniformly in azimuth. Another puts one sixth of the weight at each pole and two thirds uniformly around the equator. Their mean and covariance agree exactly.
- Smooth polar densities $(1\pm P_4(\cos\vartheta))/(4\pi)$ also share the same moments, because $P_4$ is orthogonal to all angular harmonics through degree two. Each law uses independent draws for the incoming positive and negative particles.
- Six equally weighted directions on the Cartesian axes, or eight at the vertices of a cube, give the same moment identities. Their collision forces can be found by finite sums over 36 or 64 pairs. This gives a separate example without angular quadrature.

## What the interaction measures

Give each sign the same rest-frame distribution and boost the signs by opposite rapidities $\pm r$. Let $p_i$ denote the boosted four-momenta and $E'_i$ their energies **before** boosting. For the parity-even unpolarized vector pole in the previous study, the expected outgoing four-momentum of either particle in a collision is $(p_1+p_2)/2$. Thus the normalized collision rate inferred from the force on one stream is

$$
C=\frac{m^2\nu}{n}=
\frac{1}{h\cosh r\sinh r}
\left\langle \frac{G(s)(p_{1z}-p_{2z})}{2E'_1E'_2}\right\rangle,
\qquad G(s)=\frac{6\pi\beta[0.30(s+2)]^2}
{(s-8300^2)^2+[0.30(s+2)\beta]^2},
\qquad \beta^2=1-4/s .
$$

Here $s=(p_1+p_2)^2$ in units of $m^2$, and $G$ is the invariant incoming flux times the dimensionless transport cross section. The rest-frame energy denominators arise from the invariant momentum measure. Exchanging both charge labels lets us use $r=|\operatorname{asinh}(u)|$, where $u$ is the signed proper velocity in the wall solver. Different energies and directions sample different parts of the resonance, so equality of the stress tensor leaves the average of this kernel open. [Force and angular law](../SOURCES.md)

## Two spectra with completely isotropic directions

The exact Maxwell–Jüttner equation of state gives $B=K_1(1/\Theta)/K_2(1/\Theta)$ and $A=3\Theta+B$. Any isotropic population with $\langle E\rangle=A$ and $\langle 1/E\rangle=B$ reproduces its energy and pressure. For any $1\le l<1/B$, choose two energy shells

$$
L=l,\qquad U=\frac{A-l}{1-Bl},\qquad
w_L=\frac{U-A}{U-l},\qquad w_U=1-w_L .
$$

These positive probabilities give both desired moments exactly. We compare $l=1$ to $l=1+0.9(B^{-1}-1)$. Both prescriptions use the same $n$, $\Theta$ and $u$ at the stated trajectory point. The radial comparison matches the **exact Maxwell–Jüttner** stress within the pair of examples. At startup its specific energy exceeds the earlier Taub–Mathews value by 0.0192%. Thus the angular and radial comparisons refer to two precisely specified equations of state.

| Sample on the preceding flavor-specific trajectory | $\Theta$ | $u$ | Polar/equatorial rate divided by conical rate | $(1+P_4)/(1-P_4)$ rate ratio | Second radial rate divided by first |
|---|---:|---:|---:|---:|---:|
| $\tau=6.255$ | 16.901 | −30.112 | 1.199 | 1.042 | 3.963 |
| $\tau=6.30$ | 1526.97 | −7.012 | 0.815 | 1.089 | 1.705 |
| First hot crossing, $\tau=8.43896$ | 4765.47 | 0.00009397 | 0.903 | 1.000 | 2.109 |

The ordering of the angular examples reverses between the first and second sample. The ratios in this table give specific counterexamples, with no extremization over allowed distributions. Neither witness is claimed to be the population produced by the pulse. At exactly zero drift the mean axial force vanishes by symmetry; its linear response can still depend on the energy distribution.

### Finite-width dilute spectra

An ideal energy shell concentrates probability on a surface of zero thickness. We therefore repeated the startup test with continuous isotropic spectra. Give each band center $E_i$ a half-width $d_i=0.04(E_i-1)$ and let $x=(E-E_i)/d_i$ have density $3(1-x^2)/4$ on $[-1,1]$. The two lower centers use $l=1+0.02(B^{-1}-1)$ and $l=1+0.9(B^{-1}-1)$. In each case the upper center is adjusted so that the mixture retains the same $A$ and $B$ as above. The resulting dimensionless rates are $1.94348\times10^{-6}$ and $7.50134\times10^{-6}$, a **factor of 3.860**. Refining the radial, polar and azimuthal integrals changed the rates by less than $2.7\times10^{-7}$ relative.

At the same trajectory point the model has $n/m^3=3.718\times10^{-8}$ per charge and flavor. If $w(E)$ is the unit-normalized energy density in probability, the per-spin phase-space filling of a spin-degenerate Dirac species is $f(E)=\pi^2(n/m^3)w(E)/(E\sqrt{E^2-1})$. Bounding the numerator and denominator independently on each band gives $f<8.53\times10^{-8}$ for both distributions. These finite-width examples realize the ambiguity with dilute, nonnegative occupations, fully isotropic momenta in each rest frame, and the same complete stress and number current. They also quantify how close two dynamical predictions can appear in fluid measurements while differing in their instantaneous exchange of momentum. [Energy bands and occupancy bounds](MOMENTS.json)

## Numerical checks

The code interpolates the three sample states from the replayed flavor-specific native trace and integrates over both rest-frame polar directions and the relative azimuth. It pairs directions related by inversion to condition the small-drift limit. The three-dimensional integral uses 192 Gaussian nodes in each polar cosine and 384 midpoint azimuths, compared with 128 and 256. The largest fractional change across all three states is 0.022%. At the first hot crossing an independently differentiated zero-drift integral for the two isotropic radial spectra agrees with their finite-drift results to better than $1.2\times10^{-7}$ relative. An independent reduction for the isotropic single-energy shell agrees with the three-dimensional calculation within $2.7\times10^{-6}$ relative at these three states. [Moment calculation](counterexample.py) · [Single-shell reduction](../closure/shell.hpp)

From the capsule root use the native trace generated by the common history exporter.

```sh
python -I -B history/export.py --prior prior --outdir /tmp/hdblast-native
python -I -B moments/counterexample.py --history /tmp/hdblast-native --output /tmp/MOMENTS.json
```
