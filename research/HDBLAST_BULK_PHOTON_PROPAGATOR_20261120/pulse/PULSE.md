# What happens to light when the wall moves?

A massless light wave traveling along the tuned wall can excite waves that travel in the extra dimension. We compute this effect for a wall that makes a small, prescribed excursion and returns to its original position. The incident light, spacetime and wall are specified in advance. The results below are classical and leading order in the wall displacement. Photon production from the vacuum, the energy budget of a blast and a self-gravitating wall trajectory require different calculations.

## The moving boundary is the source

Work in the fixed tuned metric $ds^2=dy^2+e^{2A(y)}(-dt^2+d{\bf x}^2)$, with $A(b)=0$ at the unperturbed wall. All lengths, including time with light speed one, are in wall-width units $L_0$. A transverse polarization with spatial wave number $k>0$ has radial weight $w=G$ and stiffness $p=G e^{2A}$. The normalized massless mode is the constant $u_0=1/\sqrt N$, with $N=\int_b^\infty w\,dy$. The other modes have continuum mass $m>0$ and frequency $\Omega=\sqrt{k^2+m^2}$. Their squared boundary amplitudes per unit mass on one retained half-line define $\rho(m)$.

Write the perturbed wall position $y=b+h(t)$ with $h=d\exp[-t^2/(2T^2)]$. Assume both the displacement and the wall speed are small. The normal Maxwell condition for this polarization is $\partial_y A_\perp+e^{-2A}\dot h\,\partial_t A_\perp=0$ at the wall. Since the incident wave is independent of $y$, the first-order boundary source is

$$p(b)\partial_y A_\perp^{(1)}(b,t)=-w_b\dot h(t)\,\partial_t A_\perp^{(0)}(b,t).$$

Projection on the fixed Neumann modes gives a forced oscillator with drive $w_b u_0 u_m(b)\dot h\,\dot a_0$, where $a_0$ is the incident massless-mode amplitude. A real oscillator initially at rest ends the drive with energy $\tfrac12|\int dt\,f(t)e^{i\Omega t}|^2$. Averaging over the phase of the incoming wave $a_0(t)\propto\cos(kt+\theta)$ gives the leading continuum energy coefficient

$$
\mathcal E(k,T)=\lim_{d\to0}\frac{E_{\rm continuum}}{d^2 E_{\rm in}}
=\frac{w_b^2}{4N}\int_0^\infty \rho(m)
\big[X_T(\Omega-k)+X_T(\Omega+k)\big]\,dm,\qquad
X_T(\nu)=2\pi(T\nu)^2e^{-(T\nu)^2}.
$$

The numerator is the energy in the excited massive sector. The moving wall can do work, so the loss from the original light wave need not equal that energy. A two-copy even sector gives the same normalized coefficient when both copies and both boundary terms are included. An overall rescaling of the gauge kinetic weight cancels out.

## Two interactions on the same geometry

The table uses $k=1$ and the central cut $\phi_b=0$. Each entry multiplies $d^2$ to give the leading fraction of incident energy. For example, at $T=1$, $d=0.01$ the leading fractions are approximately 44 and 48 parts per million. Higher orders in displacement have not been calculated.

| Pulse width $T$ | Earlier weight $G=(1-\phi^2)^2$ | Weight selected by the fixed-profile null wave, $G\propto e^{2A}$ |
| ---: | ---: | ---: |
| 0.5 | 1.521923 | 1.428530 |
| 1 | 0.440110 | 0.475633 |
| 2 | 0.100404 | 0.237554 |

For a very slow excursion at fixed $k>0$, write $\rho(m)\sim C m^{p_*}$ at threshold. The sum-frequency term is exponentially suppressed. In the difference-frequency term, $\Omega-k\sim m^2/(2k)$. Setting $x=Tm^2/(2k)$ and integrating the resulting Gaussian gives

$$ \mathcal E(k,T)\sim\frac{w_b^2}{4N}\frac{\pi C}{2}
   \Gamma\!\left(\frac{p_*+5}{4}\right)
   \left(\frac{2k}{T}\right)^{(p_*+1)/2} . $$

Both interactions admit arbitrarily low-mass waves. The earlier weight has $p_*=6.2$ and the excited continuum energy decreases as $T^{-3.6}$. The weight selected by the null wave has $p_*=1$ and this energy decreases as $T^{-1}$. At $k=1$, $T=4096$, the integrated coefficients for all three wall positions differ from their respective leading formulas by less than 0.12% and 0.62%. This limit holds at fixed nonzero incoming wave number.

## Connecting a static charge test and a moving-wall test

The fractional increase in the magnitude of the static charge potential is $D(r)=N\int_0^\infty \rho(m)e^{-mr}\,dm$. At large $r$ it approaches $NC\Gamma(p_*+1)r^{-(p_*+1)}$ to leading order. Comparing the two probes at $r=\sqrt{T/(2k)}$ cancels the threshold coefficient $C$,

$$ \lim_{T\to\infty}\frac{\mathcal E(k,T)}{D(\sqrt{T/(2k)})}
=\left(\frac{w_b}{N}\right)^2
\frac{\pi\,\Gamma((p_*+5)/4)}{8\,\Gamma(p_*+1)} .$$

Here $(w_b/N)^2$ is the square of the logarithmic derivative of the zero-mode norm with respect to the wall's proper position for a fixed global weight. At the central wall the predicted ratios are 0.00141034 and 0.324397 in the stated units. Evaluating both integrals at $T=4096$ gives ratios within 1.72% of their limits across the six cases. This relation connects two idealized measurements on the same model. Measurements of these quantities for a physical HDblast wall remain open.

## Numerical controls

The primary integration interpolates the outgoing boundary density on a logarithmic mass grid and checks grid and Gaussian quadrature refinement. Below the grid it uses the analytic threshold law. The upper mass limit $M$ satisfies $T(\sqrt{k^2+M^2}-k)=8$. For the first weight $\int\rho/m^2$ is finite, bounding the omitted contribution by this moment times the maximum of $m^2$ times the pulse kernel on $m\ge M$. For the selected weight we bound $\int\rho/(1+m^2)$ by the full Euclidean propagator at $Q=1$, itself bounded by the complementary trial flux $j(y)=\max(1-y+b,0)$. The static integral ends at $mr=50$ with a corresponding tail bound. The archive reports each high-mass tail bound separately from interpolation error.

A second code discretizes the radial energy on finite conformal intervals with reflecting boundaries, identifies and removes the massless mode by its eigenvector, and varies both length and spacing. The selected weight needs larger boxes and a correction for the finite-volume massless-mode residue. The data include the uncorrected box results. Direct time integration of the forced oscillators checks their final energy against the Gaussian frequency integral. The extrapolated box and outgoing half-line values agree within one part in a million at all 36 shared parameter points. [Pulse calculation](pulse.py) · [Finite-box calculation](box.py) · [Data](PULSE.json) · [Spectrum](../retarded/RETARDED.md) · [Static charge test](../forces/FORCES.md)
