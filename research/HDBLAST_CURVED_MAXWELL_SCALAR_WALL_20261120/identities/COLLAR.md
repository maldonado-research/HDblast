# Reconstructing the source just inside the wall

The spherical junction calculation specifies a compact Maxwell characteristic at the wall. This calculation asks which fields a distinct emitting surface must supply to produce the same wall history. It reconstructs a narrow self-gravitating annulus from the timelike wall trace. For the dose $0.123$ example, **the required incident field at the finite-distance source persists well after the wall pulse ends**. Its integrated energy after the direct pulse is small. The numerical reconstruction covers the bounce, the radiation-rich phase and the beginning of recontraction. [Timelike junction and input](../loading/loading.py) · [Spherical field equations](IDENTITIES.md)

## Coordinates, wave direction and data

Use $\ell=\lambda=1$, $\kappa_5^2=6$ and $ds^2=e^{2L}(-dt^2+dx^2)+r^2d\Omega_3^2$. The wall is at $x=0$, $t$ is wall proper time, and $x$ increases inward. On each side the outward normal toward the wall points toward **decreasing** $x$. The boundary values and first normal derivatives are

$$
L=0,\quad L_x=-(1-2e-3p),\quad r=a,\quad r_x=-a(1+e),\quad
\mathcal A=a c,\quad \mathcal A_x=-aY,\quad \mathcal A_t=aX .
$$

Here $c$ is the tangential potential amplitude divided by the area radius. The two sides have identical data. Before the pulse the vacuum history near $t=0$ has $\mu=4$, $a=\sqrt{4-t^2}$ and $\mathcal A=0$. At $x=d=0.004$, the proper elapsed time is $d\tau_d=e^{L(t,d)}dt$ and the areal radius changes. Coordinate depth $d$ specifies a timelike tube in this chart.

For $\Box_h f=e^{-2L}(f_{xx}-f_{tt})$ the radius and conformal geometry obey

$$
\Box_h r=4r+\frac{2[1-(Dr)^2]}r-\frac{\kappa_5^2r B^2}{3},\qquad
\Box_h L=4-\frac{3\Box_h r}{2r}-\frac{\kappa_5^2(D\mathcal A)^2}{6r^2}
+\frac{\kappa_5^2 B^2}{3},\qquad B=2\mathcal A/r^2.
$$

The Maxwell equation is $\mathcal A_{xx}-\mathcal A_{tt}=(r_t\mathcal A_t-r_x\mathcal A_x)/r+4e^{2L}\mathcal A/r^2$. The second-order system advances in the inward coordinate, shrinking its available interval in $t$ by $d$ at each end. In two dimensions the same null lines carry data in this inverse construction. The emitting history at $(t,d)$ depends on boundary data up through $t+d$. This dependence makes the calculation an inverse boundary reconstruction for the symmetric mode. It supplies no prediction for an independently chosen source.

For any field $f$, define $J_\pm=f_x\pm f_t$. The code transports each $J_\pm$ from the appropriate null footpoint with trapezoidal source integration. At the emitting tube,
$I_d=e^{-L}(\mathcal A_x+\mathcal A_t)/(2r)$ travels toward the wall and
$W_d=e^{-L}(\mathcal A_t-\mathcal A_x)/(2r)$ returns. Per side and unit sphere volume, integrate their positive energies with $r^3e^L I_d^2\,dt$ and $r^3e^L W_d^2\,dt$. These quantities use the local tube clock and area. Comparing an energy budget at a moving surface at another radius also requires gravitational redshift and pressure work.

## Emission, return and late drive

The wall pulse is supported on $0<t<0.0039783904$. The direct outward characteristic at the emitting tube uses $t_{\rm wall}=t+d$, so its direct support ends at $t=-0.0000216096$. Starting a post-pulse budget at $t=0$ omits a short part of the tail and gives a conservative estimate of its energy.

| Event or time interval | Reconstructed result per side |
|---|---:|
| Peak absolute incident amplitude at tube | 8.6259, at $t=-0.001175$ |
| Peak absolute return during the launch calculation | 0.55625 |
| Integrated emission, $-0.008\le t\le0.024$ | 0.493916 |
| Integrated return, same interval | 0.0196330 |
| Further return, $0.024\le t\le8.738$ | 0.00491052 |
| Integrated emission from $t=0$ through $8.738$, divided by total emission over the run | $3.74\times10^{-7}$ |
| Minimum sampled mass anywhere on the annulus | 0.06712 |
| Minimum sampled $(Dr)^2$ anywhere on the annulus | 3.4545 |

The source mass starts at 4, reaches its minimum near $t=-0.000040$, and recovers to about 0.311 near the wall bounce as returning radiation and pressure work contribute. At the tube the incoming amplitude is $5.36\times10^{-5}$ near the bounce, $1.18\times10^{-5}$ at the wall’s 90% material-radiation crossing and $6.51\times10^{-6}$ at maximum radius. It rises again during recontraction to $3.87\times10^{-5}$ by $t=8.524$. The wall's specified incident amplitude is exactly zero at all of these late times. The largest change on the common late time grid under the final twofold refinement is $1.32\times10^{-8}$. Source and wall have different ingoing traces even when the correction carries little energy.

## Independent constraints and scope

The integrator advances the three mixed field equations. For $D_\pm=\partial_t\pm\partial_x$ it independently checks both unused null constraints $D_\pm^2r-2D_\pm L D_\pm r+\kappa_5^2(D_\pm\mathcal A)^2/(3r)=0$ using fourth-order time differences. It also differentiates the reconstructed mass and compares it to the local flux plus pressure-work law. At the finest launch spacing $5\times10^{-6}$, the worst normalized null residual is below $4\times10^{-4}$, the worst absolute null residual below $0.004$ and the worst absolute mass-derivative residual below $0.017$. The normalized residual divides by one plus the sum of absolute terms of its numerical identity. At late times with spacing $10^{-4}$ the normalized null residual is below $10^{-6}$. In an independent prelaunch vacuum interval, both wave traces vanish and the mass differs from four by at most $1.1\times10^{-9}$.

This inverse calculation determines the two characteristics on one particular nearby timelike tube for the selected wall solution. It does not establish a regular centre or an emitter mechanism at smaller radius. It supports a separate forward calculation with a specified compact tube input as the direct test of whether removing the reconstructed late drive changes the history. [Reconstruction](collar.py) · [Numerical controls](IDENTITIES.json)
