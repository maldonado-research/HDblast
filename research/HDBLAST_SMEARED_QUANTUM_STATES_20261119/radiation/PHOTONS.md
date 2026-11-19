# Could light already present at the bounce take over?

Consider a round, closed universe in four-dimensional Einstein gravity containing a canonical scalar field and a uniform gas of free photons. Use $c=\hbar=k_B=1$, $x=H_0t$, $a=H_0a_{\rm physical}$, $h=H/H_0$, $u=V/(3M_{\rm Pl}^2H_0^2)$, and write $r$ for photon density in the same units as $u$. Primes below denote derivatives with respect to $x$. The photon excess is prepared explicitly. Its own gravitational effect is included. We take time-symmetric initial data with the field momentarily at rest.

## An upper bound that holds at any later radius

Photons redshift as $r=d/a^4$ in the absence of sources. The Friedmann equation for positive curvature fixes the total density to $h^2+1/a^2$. Their ratio is therefore exactly

$$
F_\gamma(x)=\frac{d}{a^2(1+a^2h^2)}\le F_{\gamma,i}\left(\frac{a_i}{a}\right)^2 .
$$

At a regular initial bounce, $h_i=0$ and $a_i''/a_i=u_i-v_i^2/3-r_i>0$, where $v=\Theta'$. This requires $F_{\gamma,i}<1/2$. At every point with $a\ge a_i$, the photons remain a minority. The proof permits a rolling field after the bounce and uses no special shape for its potential. The fraction need not decrease monotonically. At the next point with $h=0$ the bound is saturated for that radius.

If a mechanism feeds light at a rate $Q=r'+4hr$, let $C=a^4r$, so $C'=a^4Q$. To reach half the total energy at a proposed time, the integrated source must satisfy

$$
\frac{1}{a_i^2}\int_{0}^{x}a^4 Q\,dx'\ \ge\ \frac12\left(\frac{a}{a_i}\right)^2(1+a^2h^2)-F_{\gamma,i}.
$$

The right-hand side uses the radius and expansion rate of the trajectory **with** the source. The tabulated values on the unsourced runs give reference targets. Introducing a source generally changes the trajectory and requires evolving the parent and gravity together. For example, at $a/a_i=1.10$ and $F_{\gamma,i}=0.25$, at least $0.355a_i^2$ of additional comoving radiation energy is required, with a larger requirement if $h\ne0$. This bound gives a target for a conversion law. It does not supply one.

## How much light fits at launch?

For a free Maxwell field on the round unit three-sphere the physical frequencies are $L/a_{\rm physical}$, with $2(L^2-1)$ states for every integer $L\ge2$. The proper volume is $2\pi^2a_{\rm physical}^3$. At temperature $T$ the thermal excess over the conformal vacuum is determined by $y=a_{\rm physical}T$,

$$
D(y)=a_{\rm physical}^4\rho_{\gamma,\rm excess}
=\frac{1}{\pi^2}\sum_{L=2}^\infty\frac{L(L^2-1)}{e^{L/y}-1},\qquad
d=\frac{(H_0/M_{\rm Pl})^2}{3}D(y).
$$

This is an isotropic thermal occupation assigned by hand. It assumes no equilibration process. We exclude vacuum Casimir energy, the trace anomaly and any separate gravitational counterterms from both the thermal density and its gravitational source. The numerical sum checks its exponential tail, its low-temperature leading term $6e^{-2/y}/\pi^2$ and the high-temperature expansion $\pi^2y^4/15-y^2/6+y/(2\pi^2)-11/(240\pi^2)$. The linear and constant terms follow by expanding the cubic and linear sums separately.

At rest the constraint reads $u_i a_i^4-a_i^2+d=0$. For $u_i>0$ the larger positive root gives the regular bounce when $0\le d<d_{\rm crit}=1/(4u_i)$. The smaller positive root for $d>0$ gives a maximum. At equality the acceleration vanishes. See the expansion at the degenerate root in the [critical-point calculation](CRITICAL.md). The regular launch fraction is $(1-\sqrt{1-d/d_{\rm crit}})/2$, so $d/d_{\rm crit}=0.75$ corresponds to a **25%** photon fraction.

Using the archived field displacement and normalization $H_0/M_{\rm Pl}=0.00881344$, the limiting $y$ is 11.01206 for the positive-minimum branch and 11.01200 for the negative-minimum branch. At the limiting radius, $T/H_0$ is about 15.57 in either case. A state with $y=1$ starts with a photon fraction about $1.39\times10^{-5}$. These are dimensionless model scales with an assumed starting temperature.

## Evolution to the first maximum of expansion

We followed the field, photons and radius until the first maximum of the radius. In the archived model the location of that maximum changes with the photon loading. Bisection around the sampled reversal gives the photon fraction for which the maximum coincides with the first field crossing of its potential minimum.

| Location of field minimum | Initial photon fraction at coincidence | First reaches the minimum at 25%? | First reaches the minimum at 38.82%? |
| --- | ---: | --- | --- |
| $+0.5$ | 27.8906% | Yes | Radius turns first |
| $-0.5$ | 38.1433% | Yes | Radius turns first |

These thresholds apply to the tested family of time-symmetric launches at the fixed archived displacements. A different potential, release velocity or energy source calls for a new calculation. The method tests the change of event order in a sampled bracket and the field residual at its boundary. No global uniqueness proof is used. The time integration stops before leaving the domain of the fitted potential. Separate integrations evolve the logarithmic scale factor and conserved radiation coefficient, or directly evolve the radius and radiation density.

[Executable calculation](photons.py) · [Initial states, trajectories and solver checks](PHOTONS.json)
