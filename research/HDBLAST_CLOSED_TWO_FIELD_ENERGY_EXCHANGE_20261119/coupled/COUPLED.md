# A rolling field can give energy to a classical wave or take it back

At the same energy scale for the weak second field, reversing which part of its oscillation is excited at the start reverses the direction of direct energy transfer by first passage. Reversing the ripple in the rolling field also changes the sign. These calculations follow both fields and their gravity, with no prescribed decay rate. They start with a classical wave that already carries energy.

| First contact of parent | Daughter launch | Net direct transfer to daughter by first contact, % of daughter’s launch energy | Net change of daughter energy including pressure work, % of launch energy | Direct transfer over the 0.01 time window centered on first contact, % of launch energy |
|---|---|---:|---:|---:|
| Poles | Velocity | +0.8511 | −4.3809 | +0.1529 |
| Poles | Displacement | −1.4915 | −13.5603 | −0.0987 |
| Equator | Velocity | −7.1698 | −9.4587 | −0.8472 |
| Equator | Displacement | +19.4677 | +5.3892 | +2.5860 |

The table uses a daughter seed $s=0.01$ and signed parent amplitude $A=\pm0.001$. All four launch energies are approximately $1.6667\times10^{-5}$ in units of the initial stationary round volume times the hilltop energy density. Thus the largest positive entry corresponds to $3.245\times10^{-6}$ in those reference energy units. At all four first contacts the daughter contributes only $1.04$–$1.21\times10^{-5}$ of the instantaneous total energy of both scalar fields. Positive exchange can accompany decreasing daughter energy because the wave also does pressure work. All energies and rates below are whole-sphere quantities in the proper-time slicing of the initial normals. The unseeded classical daughter stays exactly zero.

## Action and local evolution

In the archive’s dimensionless canonical variables, the field Lagrangian is proportional to $-[(\nabla q)^2+Z(q)(\nabla X)^2]/2-3u(q)$. The parent is $q=\Theta-\Theta_{\rm top}$. The daughter $X$ is massless and has shift symmetry. The physical calculation uses $Z=f(\Theta_{\rm top})/f(\Theta)$, and the control uses $Z=1$. The constant in the numerator fixes the normalization of $X$ at the hilltop. The calibrated potential and matter kinetic factor are copied from the preceding study. We use the potential whose minimum is at $\phi_m=-0.5$ in the original field coordinate. Its canonical offset is $q_m=0.29693936$.

Write $\pi=\dot q$, $v=\dot X$, and let $g=D_{\hat r}q$, $c=D_{\hat r}X$ be physical radial gradients. If $\theta$ is the volume expansion and $\Delta$ the spatial Laplacian, the field equations are

$$
\dot\pi=\Delta q-\theta\pi-3u_{,q}+\tfrac12Z_{,q}(v^2-c^2),\qquad
\dot v=\Delta X-\theta v+(\ln Z)_{,q}(gc-\pi v).
$$

The Einstein evolution has the archived synchronous regular radial form. Its gradient anisotropy contains $q_y^2+ZX_y^2$ in $y=\cos\chi$. Its angular evolution contains the parent potential. The energy and momentum constraints contain the additional positive energy $Z(v^2+c^2)$ and the momentum proportional to $ZvX_y$. The independent Raychaudhuri check uses $\pi^2+Zv^2$. Thus the daughter changes the parent acceleration, anisotropic stress and constraint data from the outset.

## Matched initial stress

The initial slice has $K_{ij}=0$, $q=A U_2(\cos\chi)$, and $\pi=0$, where $U_2=4\cos^2\chi-1$. We tested two daughter phases.

* For a velocity seed, $X=0$ and $v=s h U_2/\sqrt Z$, with $h$ the hilltop Hubble scale.
* For a displacement seed, $v=0$ and $dX/dU_2=s/(\sqrt{8Z})$. The number 8 is the round-sphere Laplace eigenvalue of this harmonic, so the two phase choices have the same leading round-slice launch energy. An irrelevant integration constant is chosen by $X(U_2=0)=0$.

Each phase pair at constant and physical $Z$ has identical pointwise stress on its initial slice. Between the two phases the spatial distribution and the kinetic and gradient shares differ. No classical daughter momentum flux is present initially. The velocity seed is compatible with vanishing geometric extrinsic curvature because its initial spatial gradient vanishes. It is not a stationary field configuration. We solve the Hamiltonian constraint afresh in each case. In the conformal equation from the parent study replace $\Delta u/u_{\rm top}$ by $\Delta u/u_{\rm top}+s^2U_2^2/6$ for the velocity seed. For the displacement seed add $(s^2/8)(\partial_\chi U_2)^2$ to the dimensionless round-gradient source. The initial metrics of each constant- and physical-coupling pair coincide.

## Energy, work and shift charge

For either scalar let $K$ and $G$ be its local kinetic and gradient contributions. For the parent they are $(\pi^2,g^2)/6$, and for the daughter they are $Z(v^2,c^2)/6$. For potential $V_i$ equal to $u$ or $0$, the density and directional pressures are $\rho_i=K+G+V_i$, $p_{r,i}=K+G-V_i$, and $p_{t,i}=K-G-V_i$. Write $E_i=\int\rho_i\,dV$, and $W_i=\int dx\int(\alpha p_{r,i}+2\beta p_{t,i})\,dV$, where $\alpha,\beta$ are local radial and tangential expansion rates. The whole three-sphere has no boundary flux. The field equations and integration by parts then give

$$
T=\int dx\int dV\,\frac{Z_{,q}\pi}{6}(c^2-v^2),\qquad
\Delta E_X+W_X=T,\qquad
\Delta E_q+W_q=-T.
$$

The sign of the integrand depends both on the motion of the parent and on the difference between gradient and kinetic energy in the daughter. A single nonnegative friction coefficient cannot reproduce all four integrated signs. The integrated daughter shift current $Q_X=\int Zv\,dV$ is separately conserved. The pressure-to-energy ratio of the daughter, $(K-G/3)/(K+G)$ using volume-integrated energies, is $-0.303$ or $+0.844$ at polar first contact and $-0.320$ or $+0.992$ at equatorial first contact for velocity or displacement launch respectively. These are instantaneous classical coherent-field diagnostics. None of these runs measures a temperature or includes scattering, a quantum vacuum, or spontaneous emission.

## Backreaction, convergence and reproduction

To compare the two couplings at a common proper time we use $x=0.505$ for polar and $0.605$ for equatorial data. In the equatorial velocity-seed run the matter coupling shifts the volume expansion $H_V$ by $+4.16\times10^{-5}$ relative to the constant-coupling solution. In the equatorial displacement-seed run the difference is $-3.83\times10^{-5}$. With the half-amplitude seed the changes are $+1.04\times10^{-5}$ and $-9.58\times10^{-6}$. At both amplitudes the relative transfer as a fraction of seed energy stays close to its tabulated value. The largest change in that fraction is 0.0025 percentage points. These effects pertain to this classical seed and normalization.

The calculation uses 18 primary configurations at maximum even degree 64, with 192 angular integration nodes and maximum time step 0.004. For each of the eight nonzero full-seed choices, a separate degree-96 run with 256 nodes and a separate time-step-0.002 run check angular and temporal convergence. Every trial right-hand-side call encloses the angular parent polynomial within the fitted field range using an upper bound on its second angular derivative. A pole-inclusive grid brackets the first minimum crossing between accepted time steps. The $x_*\pm0.005$ records lie on the same uninterrupted solution. A finer spatial quadrature evaluates the energies. Independent Gauss–Legendre time quadratures compare integrated pressure work and transfer with endpoint energy changes.

Across the spatial checks the maximum first-crossing time change is below $5\times10^{-14}$ and the largest fixed-time transfer change below $1.4\times10^{-18}$ reference energy units. The degree-96 checks have normalized Einstein constraint and Raychaudhuri residuals below $9\times10^{-10}$ on 31 time slices. Across primary cases the parent and daughter energy balances close to better than $2.6\times10^{-12}$ and $4\times10^{-18}$ respectively in reference energy units. The shift-charge drift is below $8\times10^{-16}$ in the reported units. The two zero-seed first crossings reproduce the archived single-field event times to better than $10^{-12}$.

From the study root run `OPENBLAS_NUM_THREADS=1 python -IB coupled/coupled.py --prior prior --output coupled/COUPLED.json`. [Results and numerical controls](COUPLED.json) · [Source](coupled.py)
