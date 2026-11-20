# Can incoming radiation populate the massless mode?

A compact pulse launched with its nearest edge at conformal distance six reaches a static reflecting wall and exerts a force in the model. The portion already in the massless mode keeps exactly the same energy throughout. A pulse prepared to have zero overlap reaches the wall too, then reflects without creating a massless excitation. This result applies to either gauge coupling on the fixed central slice. It gives a control for the [responding wall calculation](../response/).

## A bound on what the source prepares

Use units with $c=L_0=1$, the central cut $b=0$, $A(0)=0$ and one copy of the exterior. The proper distance is $y$ and the conformal distance obeys $dz=e^{-A}dy$. For the earlier gauge law $G=(1-\phi^2)^2=\mathrm{sech}^4 y$ or the selected law $G=e^{2A}$, set $w=e^A G$. A transverse photon of tangential wave number $k$ obeys

$$w(\ddot q+k^2q)=(wq_z)_z,\qquad q_z(t,0)=0.$$

Its energy per unit three-volume along the wall, using a spatial Fourier basis of unit mean square, is $E=\tfrac12\int_0^\infty w(\dot q^2+q_z^2+k^2q^2)dz$. The field and its velocity have compact radial support at each finite time. Integrating the equation therefore gives an isolated oscillator for the weighted mean $M=\int wq\,dz$,

$$\ddot M+k^2M=0,\qquad E_0=\frac{\dot M^2+k^2M^2}{2N},\qquad N=\int_0^\infty w\,dz=\int_0^\infty G\,dy.$$

Here $E_0$ is the energy in the normalized constant radial eigenfunction $1/\sqrt N$. If the source prepares initial field and velocity only on $S$, Cauchy–Schwarz bounds both moments and gives

$$\frac{E_0}{E}\le\frac{\int_S w\,dz}{N}.$$

The bound is sharp as a supremum for general compact data of a given support interval. Set the initial field to zero and approximate constant velocity on that interval by smooth functions. Directional launch and smoothness can reduce the fraction attained by any specified example. In a compact pulse with nonzero $M$, the continuum contribution initially cancels the massless mode outside the source. A nonzero mode projection is compatible with a strictly retarded wall signal.

| Gauge coupling | $N$ | Maximum possible initial pole fraction for $6<z<10$ | Upper bound for any launch beyond $z=32$ |
| --- | ---: | ---: | ---: |
| Earlier $\mathrm{sech}^4y$ | $2/3$ | 0.008796% | $2.0863\times10^{-9}$ |
| Selected $e^{2A}$ | 1.03577126 | 2.790413% | $2.4951\times10^{-3}$ |

The last column uses fractions of one. At large launch radius $R$, $A\sim-5y/9$ and $z\propto e^{5y/9}$. The exterior fraction decays as $R^{-36/5}$ for the earlier coupling and $R^{-2}$ for the selected coupling. Both limits exclude finite static capture from a pure continuum input.

## A compact incidence experiment

Set $k=1$, $x=(z-8)/2$ and $B(x)=\exp[1-1/(1-x^2)]$ for $|x|<1$, extended by zero. The raw field is $q=aB$. The balanced field is $q=aB(x^2+c_1x+c_0)$, with $c_0,c_1$ determined by $\int wq\,dz=\int wq_z\,dz=0$. In both cases set $\dot q=q_z$ initially and set $a$ by $E=1$. The instantaneous radial energy current is $-w\dot q q_z\le0$ everywhere at launch. This prescribes inward flux for a broadband disturbance. The pulse can disperse because the background is inhomogeneous and $k$ is nonzero. The coefficients, initial moments and complete sampled wall traces are included in the [numerical record](SCATTERING.json).

At the fixed boundary the conformal normal pressure per unit tangential volume is $P_b=w(0)(\dot q_b^2-k^2q_b^2)/2$. Positive $P_b$ pushes the wall towards decreasing $z$. Define the signed pressure impulse by $I=\int_0^{26} P_b\,dt$. A static reflecting boundary has zero energy flux across it at all times.

| Gauge coupling | Launch | Initial massless energy fraction | Peak $P_b$ | Time of sampled peak | $I$ |
| --- | --- | ---: | ---: | ---: | ---: |
| Earlier | Raw | $2.80152\times10^{-5}$ | 0.61371 | 6.85 | 0.65334 |
| Earlier | Balanced | 0 | 1.98990 | 9.85 | 1.51678 |
| Selected | Raw | 0.00820510 | 0.65408 | 10.70 | 1.06937 |
| Selected | Balanced | 0 | 1.48471 | 9.90 | 1.76405 |

The continuum front can first reach the wall at $t=6$. For the balanced packets, the maximum absolute wall field reaches 0.44843 and 0.69329 respectively despite zero massless energy. The semi-infinite continuum also has unit flux reflectance and a coupling-dependent return phase. A separate calculation compares identical asymptotic incoming spectra on both branches, which removes differences from the branch-dependent preparation of the balanced pulses. [Reflection phases and delays](PHASE.md)

The inhomogeneous background also exchanges radial momentum with the pulse. Write $p=-\int w\dot q q_z\,dz$. On a domain ending at $L$ the integrated balance is

$$\Delta p=\int dt\,[P_b-P_L]+\frac12\int dt\,dz\ w'(\dot q^2-q_z^2-k^2q^2).$$

This explains why twice the initial radial momentum magnitude alone does not give the force impulse. For the selected balanced packet $p$ changes from $-0.92555$ to $+0.92866$ and the distributed background contributes $+0.09017$.

## Numerical controls and scope

The calculation uses positive quadrature weights on degree-four Legendre–Gauss–Lobatto elements, the variational Neumann stiffness and adaptive eighth-order time integration. It compares three nominal node spacings, $0.04,0.02,0.01$, and samples the wall every 0.025 up to $t=26$. The far surface at $z=42$ stays outside the light cone of the entire initial support. Diagonalizing the full symmetric mass-scaled wave operator independently reproduces the coarsest semidiscrete wall pressure within $7\times10^{-9}$. On the two finest meshes the peak-normalized pressure curves differ by at most 0.312% and the wall impulses differ by less than $2.32\times10^{-6}$. On the finest mesh the integrated momentum-balance residual is at most $2.49\times10^{-6}$, total energy drift is below $5\times10^{-13}$ at stored precision, and the oscillator moment differs from its predicted trace by at most $1.5\times10^{-11}$. Residuals printed as zero are rounded to 12 decimal places. The stored near-wall energy is a quadrature estimate over $0\le z\le1$. [Replay code](scattering.py)

The experiment tracks one tangential Fourier component. Finite *total* probe energy can also be achieved. For example, choose a spatial vector profile $\mathbf e(\mathbf x)=\nabla\exp[-|\mathbf x|^2/(2R^2)]\times\hat{\mathbf x}_3$. With temporal and radial components of the gauge potential zero, choose $\mathbf A=q(z)\mathbf e$ and $\dot{\mathbf A}=g(z)\mathbf e$ at launch. The vector is divergence-free and satisfies the transverse Maxwell constraint. Its squared norm on three-space is $\pi^{3/2}R$ and $\|\nabla\mathbf e\|^2/\|\mathbf e\|^2=5/(2R^2)$. The total energy is therefore $\pi^{3/2}R\int w[g^2+q_z^2+5q^2/(2R^2)]dz/2$. Radial balancing applies to every Fourier component. This establishes finite-energy initial data directly. The numerical impulses above are for the single $k=1$ component.

The selected gauge weight gives the same radial wave operator as the transverse-traceless tensor sector on this fixed background. The pressure used here comes from the Maxwell action. Gravitational backreaction, wall motion, emission by wall matter and an expanding universe each require their own dynamical equations.
