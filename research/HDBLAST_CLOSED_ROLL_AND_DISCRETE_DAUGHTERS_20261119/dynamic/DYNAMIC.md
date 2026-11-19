# Finite-band daughter propagation on the evolving slices

The evolving inhomogeneous geometry produces off-diagonal mixing in a finite daughter band. The matter kinetic weight changes that mixing substantially. These results use a specified finite-time vacuum convention: the instantaneous Hamiltonian of the original daughter field in the archived Einstein-frame geodesic proper time. They do not define a unique particle number in curved spacetime, a thermal radiation energy, a UV-complete yield, or a backreaction calculation.

The source and full numerical controls are [spectator.py](spectator.py) and [DYNAMIC.json](DYNAMIC.json). Reproduce from the study root with:

```sh
OPENBLAS_NUM_THREADS=1 python -IB dynamic/spectator.py --gravity gravity --output /tmp/DYNAMIC.json
```

## Field, time and spatial projection

The supplied background has unit lapse and zero shift,

$$
ds_E^2=-dx^2+a(x)^2[e^{2L_d}d\chi^2+e^{2M_d}\sin^2\chi\,d\Omega_2^2],
\quad L_d=m+\sin^2\chi\,d,\quad M_d=m.
$$

Here $x$ is the dimensionless proper-time coordinate of the gravity capsule. The fields $m,d,q$ come from its full nonlinear Einstein--scalar evolution, with $\Theta=\Theta_{\rm top}+q$. For a real, massless scalar minimally coupled to the matter metric, $g_J=f(\Theta)^{-1}g_E$, the Einstein-frame action is

$$
S_\chi=-\frac12\int d^4x\sqrt{-g_E}\,B(\Theta)g_E^{\mu\nu}
\partial_\mu\chi\partial_\nu\chi,
\qquad B(\Theta)=\frac1{f(\Theta)}.
$$

The implementation uses this corrected $B$, without replacing it by $f$, and checks its positivity. A second propagation on exactly the same evolved geometry sets $Z=1$; the physical matter calculation sets $Z=B$. Each calculation starts in its own initial Hamiltonian vacuum under the same stated prescription.

We retain one $\ell=2$, fixed-$m$ block of the real scalar, using reflection-even round radial harmonics with $n=2,4,\ldots,n_{\max}$. All five angular $m$ blocks evolve identically. For a fixed radial basis $R_i$, normalized on the physical round $S^3$ together with unit-normalized angular harmonics, its projected action is

$$
S=\frac12\int dx\,[\dot X^T M\dot X-X^T KX],
$$
$$
M_{ij}=a^3\int Z e^{L_d+2M_d}R_iR_j\,d\mu,\qquad
K_{ij}=a\int Z[e^{2M_d-L_d}R'_iR'_j+
6e^{L_d}R_iR_j/\sin^2\chi]\,d\mu.
$$

Here $d\mu=\sin^2\chi\,d\chi$; the harmless constant normalization of the basis is accounted for by the normalized Gauss--Chebyshev quadrature. The momentum is $P=M\dot X$. We solve the Hamiltonian equations directly,

$$
\dot X=M^{-1}P,\qquad \dot P=-KX.
$$

This formulation retains spatial gradients of the kinetic weight through the weak spatial operator. No homogeneous oscillator potential is substituted for the actual inhomogeneous field equation. The metric and parent field are prescribed by the gravity solution and receive no daughter feedback.

## Explicit vacuum and Bogoliubov convention

At each endpoint, real generalized eigenvectors obey
$Kv_j=\omega_j^2Mv_j$, $v_i^TMv_j=\delta_{ij}$.
The incoming positive-frequency mode matrices are

$$
U(0)=V(2\Omega)^{-1/2},\qquad
P_U(0)=-iMV(\Omega/2)^{1/2}.
$$

At the final time, let $Q=V_f^TM_fU$ and $\Pi=V_f^TP_U$. The reported coefficients are

$$
\alpha=\sqrt{\Omega_f/2}\,Q+i(2\Omega_f)^{-1/2}\Pi,
\qquad
\beta=\sqrt{\Omega_f/2}\,Q-i(2\Omega_f)^{-1/2}\Pi.
$$

Thus $\beta_{ji}$ is the negative-frequency component of initial positive mode $i$, expressed in final mode $j$. The occupation diagnostic is $N_j=\sum_i|\beta_{ji}|^2$ over the stated initial band. Ordered radial eigenmodes inherit the labels $n=2,4,\ldots$. Changing time, the canonical field, or the finite-time positive-frequency prescription can change these occupations. Off-diagonal $\alpha$ measures mixing within positive frequencies; it is recorded separately from $\beta$.

## Results

The primary band is $n_{\max}=24$. At $x=0.12$, the physical kinetic-weight calculation gives the following dimensionless occupations per angular block. The last column is precisely $\sum_{i,j\in\{2,4,6,\ldots,24\},\ i\ne j}|\beta_{ji}|^2$: the 132 off-diagonal elements of a 12-by-12 radial matrix, for one fixed $\ell=2,m$ block. It is a sum, with no division by the number of modes and no factor of five for the angular multiplicity.

| Potential minimum | Initial quadrupole | $N_2$ | $N_4$ | $N_8$ | $\sum_{i\ne j}|\beta_{ji}|^2$ |
|---|---:|---:|---:|---:|---:|
| $\phi_m=+0.5$ | $+0.001$ | 0.0000482339 | 0.0000504885 | 0.0000440057 | 0.00000158548 |
| $\phi_m=+0.5$ | $-0.001$ | 0.0000525284 | 0.0000454935 | 0.0000364349 | 0.00000153766 |
| $\phi_m=-0.5$ | $+0.001$ | 0.0000489840 | 0.0000493737 | 0.0000422848 | 0.000000578219 |
| $\phi_m=-0.5$ | $-0.001$ | 0.0000516382 | 0.0000462910 | 0.0000376228 | 0.000000595494 |

For example, the $+0.5$ minimum and $+0.001$ seed give $(N_2,N_4,N_8)=(0.0000502873,0.0000477720,0.0000398470)$ with the same geometry and $Z=1$, and an off-diagonal sum $2.7250\times10^{-11}$. The JSON also records the intermediate $x=0.20$ results for every sign, branch and kinetic choice.

On the negative-minimum branch, $x=0.47$ remains before the first minimum passage identified by the separate gravity event calculation. The additional late diagnostic is:

| Seed | Kinetic choice | $N_2$ | $N_4$ | $N_8$ | Off-diagonal $\beta$ sum |
|---:|---|---:|---:|---:|---:|
| $+0.001$ | $B=1/f$ | 0.00653243 | 0.00438613 | 0.000999778 | 0.000960156 |
| $+0.001$ | $Z=1$ | 0.00777468 | 0.00356165 | 0.000530275 | 0.0000285430 |
| $-0.001$ | $B=1/f$ | 0.00904796 | 0.00308946 | 0.000600283 | 0.00291464 |
| $-0.001$ | $Z=1$ | 0.00777760 | 0.00355259 | 0.000521364 | 0.0000523926 |

The round de Sitter control already has nonzero diagonal finite-time occupations under this convention. At $x=0.12$, $(N_2,N_4,N_8)=(0.0000502944,0.0000477847,0.0000398664)$, with zero off-diagonal mixing to numerical precision. Consequently, the total finite-time occupation cannot be identified entirely with the inhomogeneous quadrupole or with an oscillating parent.

## Numerical controls and their scope

The zero-duration calculation gives a largest occupation $1.74\times10^{-31}$. Canonical commutators and both Bogoliubov column identities are preserved to $7.11\times10^{-15}$.

For the homogeneous zero-seed case we independently propagate the exact minimally coupled closed-de-Sitter solution. With $\eta=\arctan\sinh(hx)$, its mode solutions are proportional to
$(\sin\eta-i(n+1)\cos\eta)e^{-i(n+1)\eta}$.
The relative transfer-matrix disagreement is below $8.42\times10^{-14}$, for both kinetic choices. Their off-diagonal squared sums are below $1.42\times10^{-28}$.

The early results compare daughter cutoffs $n_{\max}=12,18,24$; gravity degrees 32 and 48; daughter quadrature orders 128 and 192; and DOP853 maximum time steps 0.006 and 0.003 with tightened tolerance. The late $x=0.47$ calculation uses gravity degree 48 against 80, daughter quadrature 160 against 224, and daughter bands 12,18,24,30. Across all reported low modes, time, spatial quadrature and geometry refinement change occupations by less than $4.00\times10^{-16}$ in absolute units. At the late endpoint, the largest low-mode change from band 12 to 24 is $2.28\times10^{-8}$; from 18 to 24 it is $2.10\times10^{-13}$; from 24 to 30 it is $5.12\times10^{-17}$.

These are convergence statements for the listed low modes within one angular sector. The calculation does not sum all angular sectors or establish convergence of the full quantum stress tensor. It shows that dynamical finite-band mode mixing can be followed on the actual evolved geometry under an explicit endpoint convention.
