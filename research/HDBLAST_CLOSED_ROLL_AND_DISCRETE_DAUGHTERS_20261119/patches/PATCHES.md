# How much of the slice passes the minimum?

Four radial solutions cross the minimum of the negative-minimum potential. Reversing the initial quadrupole changes the dimension of the set that first touches the minimum. Two poles touch for positive amplitude. An entire equatorial two-sphere touches for negative amplitude. The distinction produces different powers of elapsed proper time in the volume that has rolled beyond the minimum.

| Initial amplitude | First site | First contact $x_*$ | Volume in the well window at contact | Volume beyond the minimum at $x_*+0.01$ | Leading coefficient $C$ |
|---:|---|---:|---:|---:|---:|
| +0.001 | Both poles | 0.48839541 | 3.4072% | 0.90153% | 9.57209 |
| +0.0005 | Both poles | 0.55185054 | 3.4226% | 0.89498% | 9.49830 |
| −0.001 | Equator | 0.58913424 | 31.2245% | 20.27595% | 2.08655 |
| −0.0005 | Equator | 0.65307973 | 31.2757% | 20.24488% | 2.08296 |

Each volume fraction uses the proper measure of its own synchronous slice. The four times start at the respective time-symmetric initial slices. The well window is the interval between the two archived field markers $\phi_m\pm0.1$. Going beyond the minimum means $q>q_m=0.2969393579$ in the canonically normalized offset field. [Numerical results](PATCHES.json)

## First-contact volume law

Let $\xi$ be angular distance from a first-contact site. In these solutions $v=\partial_x q>0$ and $c=-\partial_\chi^2q>0$ there. Reflection symmetry gives

$$
q-q_m=v\,\delta x-\tfrac12c\,\xi^2+O(\delta x^2+\delta x\,\xi^2+\xi^4).
$$

The level surface therefore has angular distance $k\sqrt{\delta x}[1+O(\delta x)]$ with $k=\sqrt{2v/c}$. Write $J=\exp(3m+\sin^2\!\chi\,d)$ for the determinant factor after removing the round spatial volume density. Let $I=\int_0^{\pi/2}J\sin^2\chi\,d\chi$ at contact. The factors from the azimuthal area and the two hemispheres cancel in the volume fraction. For a polar cap $\sin^2\chi\sim\xi^2$. For an equatorial belt $\sin^2\chi\sim1$. This gives

$$
f_{q>q_m}=C(\delta x)^{p/2}[1+O(\delta x)],\qquad
C=\frac{J_* k^p}{p I},\qquad
p=\begin{cases}3 &\text{polar contact},\\1 &\text{equatorial contact}.\end{cases}
$$

The code evaluates $v,c,J_*,I$ from the contact slice. Separately it finds the field-level root on each subsequent slice, integrates the proper volume on either side and compares the measured fractions to this formula. From $\delta x=10^{-6}$ to $10^{-5}$ the measured log slopes are 1.49998 and 0.49999, respectively. At $10^{-6}$ the coefficients inferred from the volumes differ from the local prediction by at most 6.1 parts per million. The relation applies to the specified geodesic slicing and a smooth first contact with nonzero velocity and transverse curvature. Simultaneous contact over the equator uses the imposed symmetry. The calculation concerns the geometry of a scalar level set.

The same local calculation has a coordinate-independent form on a specified foliation. Suppose the set of simultaneous first contacts is a smooth submanifold $S$ of codimension $p$ on the contact slice. Assume the initial-slice normal observers see $v>0$ there and the negative spatial Hessian $B=-D^2q$ is positive definite on the directions normal to $S$ within the slice. In proper orthonormal coordinates $z$ along those directions, the crossed region has leading section $z^TBz<2v\,\delta x$. Integrating the volume of these ellipsoids yields

$$
f_{q>q_m}\sim
\frac{\omega_p(2\delta x)^{p/2}}{V_*}
\int_S \frac{v^{p/2}}{\sqrt{\det B_\perp}}\,dA_S ,
$$

where $\omega_p$ is the Euclidean unit $p$-ball volume and $V_*$ the full proper volume at contact. At isolated points the integral is a sum. This is the local volume law for a nondegenerate transverse maximum. It explains why a symmetry-breaking perturbation that reduces the first-contact set from a whole equatorial sphere to isolated points can change the earliest exponent from $1/2$ to $3/2$. Testing such nonspherical data requires a separate evolution.

## Can the field already be confined to the well window?

With $s=\partial_{\hat r}q$ in the orthonormal frame of the initial-slice normal observers, the scalar $X=-\nabla_\mu q\nabla^\mu q=\Pi^2-s^2$ determines whether a scalar rest frame exists. At first contact, $X$ ranges from 4.005 to 4.778 across the four occupied well regions. In the scalar rest frame the energy above the potential minimum is $(u-u_m)+X/6$ in the archive units. The lower of the two edge energies above the minimum is 0.11349372. Rest energy throughout the occupied part of the window exceeds that level by factors ranging from 6.667 to 7.333. After another 0.01 proper time the corresponding range is 6.458–7.259. The old energetic criterion for remaining within the window is unsatisfied on both slices. Further motion and energy transfer determine whether confinement occurs later. The window edges are finite field markers in the previous study.

Spatial gradients vanish at the symmetry sites, but the scalar Laplacian there is appreciable. It ranges from −1.309 to −1.105 at the first crossings, where the local speeds are 2.131–2.183. Even at those sites the scalar acceleration receives a spatial-derivative contribution. A timelike gradient means the passing minimum level surface is spacelike wherever $X>0$ on it. The positions of successive crossings therefore describe the phase of the field.

## Regional energy balance

Take the fixed comoving subdomains consisting of two polar caps with $\chi\le\pi/3$ or $\chi\ge2\pi/3$ and their complementary belt. The initial quadrupole changes sign on their shared interfaces. Let $\alpha,\beta$ be the local radial and tangential expansion rates. With $K=\Pi^2/6$ and $G=s^2/6$, the normal-frame energy, directional pressures and outward-oriented radial flux are

$$
\rho=K+G+u,\quad p_r=K+G-u,\quad p_t=K-G-u,\quad j_{\hat r}=-\Pi s/3.
$$

Here the flux component is defined along increasing $\chi$ and the domain boundary supplies the outward orientation. Stress-energy conservation on a fixed comoving domain in unit-lapse, zero-shift coordinates gives

$$
E_D(x)-E_D(0)
+\int_0^x dx'\int_D(\alpha p_r+2\beta p_t)\,dV
+\int_0^x dx'\int_{\partial D} j_{\rm out}\,dA=0.
$$

Up to first contact, integrated cap flux is inward for positive seeds, about $2.00$–$2.02\times10^{-5}$, and outward for negative seeds, about $5.13$–$5.18\times10^{-4}$. Units are the initial stationary round-sphere volume times the hilltop energy density. The largest absolute residual across the eight regional identities is $2.18\times10^{-12}$ in those units. Time integration uses 64-point Gauss–Legendre quadrature separate from the evolution. The stress-energy identity checks regional energy transport independently of the Einstein constraint residuals.

The spatial and temporal comparison uses maximum even degrees 64 and 96, projection node counts 256 and 384, maximum time steps 0.004 and 0.003, and 48 and 64 time quadrature points. The largest relative shifts in the predicted coefficient and measured crossed fractions are $4.45\times10^{-10}$ and $4.36\times10^{-10}$. For the degree-96 evolutions through $x_*+0.01$, the largest sampled off-grid Hamiltonian, physical momentum and Raychaudhuri errors in the archived normalizations are $1.96\times10^{-10}$, $5.35\times10^{-10}$ and $4.04\times10^{-10}$. [Calculation](patches.py) · [Full records](PATCHES.json)
