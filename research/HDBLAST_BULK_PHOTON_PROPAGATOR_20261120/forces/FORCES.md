# How far must two wall charges be apart to follow Coulomb’s law?

Two small conserved test charges probe the bulk photon without shaking the wall. On each flat static leaf, normalize the four-dimensional Coulomb force to the massless photon pole measured at large separation. Massive modes in the fifth dimension add a positive correction to its magnitude. With the specified Maxwell weight, the correction to the force falls below 1% at separations between 2.69 and 4.15 times the wall length scale, depending on where the wall sits. A physical length in metres would require an independently calibrated wall scale.

| Wall field | Potential excess at $R=2L_0$ | Force excess at $R=2L_0$ | Potential reaches 1% at $R/L_0$ | Force reaches 1% at $R/L_0$ |
|---:|---:|---:|---:|---:|
| 0 | 1.24% | 5.45% | 2.131 | 3.271 |
| +0.5 | 3.57% | 15.76% | 2.804 | 4.151 |
| −0.5 | 0.60% | 2.65% | 1.711 | 2.687 |

## From the brane spectrum to a force

These ratios assume conserved external charges coupled minimally to the boundary gauge field, a flat wall with warp factor one, a fixed geometry, and no independent brane kinetic term. Write $R$ in units of $L_0$, $m$ in inverse $L_0$, $N=\int G\,dy$ for the one-copy mode norm, and $\rho(m)$ for the positive continuum density per unit mass. For static conserved sources the longitudinal polarization decouples. Three-dimensional Fourier transformation of $1/(Q^2+m^2)$ gives $e^{-mR}/(4\pi R)$, while the massless pole has residue $1/N$. Thus the relative excesses in the potential and force magnitudes are

$$
\eta_V=\frac{|V|}{|V_0|}-1=N\int_0^\infty \rho(m)e^{-mR}\,dm,\qquad
\eta_F=\frac{|F|}{|F_0|}-1=N\int_0^\infty \rho(m)(1+mR)e^{-mR}\,dm .
$$

Here $|V_0|$ and $|F_0|$ follow $1/R$ and $1/R^2$. Both excesses decrease with $R$. The force correction is larger because differentiation also samples the exponential attenuation. The overall coupling and the common factor from a second identical copy cancel in each ratio. Within this model an independently normalized 1%-precision Coulomb test at separation $R_{\rm exp}$ would bound $L_0$ by $R_{\rm exp}$ divided by the appropriate force crossing in the table. No such experiment or unit conversion is used here. [Spectral density and its normalization](../retarded/RETARDED.md) · [Euclidean response](../euclidean/EUCLIDEAN.md)

## Distance limits

The low-mass density obeys $\rho(m)\sim C m^{6.2}$ for the specified weight. With $a=7.2$, the Laplace integrals give

$$
\eta_V \sim N C\,\Gamma(a) R^{-a},\qquad \eta_F \sim (1+a)N C\,\Gamma(a)R^{-a}.
$$

The absolute long-distance correction to the Coulomb force consequently decays as $R^{-9.2}$. The numerical force integral at $R=256L_0$ agrees with its leading asymptote within 0.053% on all three leaves. This power refers to the assumed squared-wall-gradient Maxwell weight. Changing its tail changes the power.

At the short-distance end, the smooth wall looks locally five dimensional. The half-line propagator approaches $1/(G_b Q)$ and its corresponding continuum density tends to $2/(\pi G_b)$. The resulting leading potential magnitude is proportional to $1/R^2$, with $|V|/|V_0|\sim 2N/(\pi G_b R)$. Any microscopic brane kinetic term, finite source size or ultraviolet limitation on the five-dimensional field theory can change this ideal short-distance prediction.

## Numerical checks

The calculation integrates the continuum on panels from $m=0$ to $M=16$, with Gauss–Legendre orders 8 and 16. The threshold is integrable and a second radial evolution checks the spectrum in both integration directions, including a differently placed outer boundary. For the high-mass tail we use the finite Euclidean spectral moment $N\int \rho(m)/m^2\,dm=J/N$, derived independently from the zero-frequency radial equation. When $MR\ge 1+\sqrt{3}$, both $m^2e^{-mR}$ and $m^2(1+mR)e^{-mR}$ decrease throughout the omitted mass range. The omitted force excess is therefore bounded by $(J/N)M^2(1+MR)e^{-MR}$, and omitting the $(1+MR)$ factor gives the potential bound. The JSON gives the tail bounds, integration refinement and the separate one-percent radii from the truncated integrals and their maximal omitted tails. The analytic tail bound does not certify floating-point integration of the retained band. [Numerical results](FORCES.json) · [Replay](forces.py)
