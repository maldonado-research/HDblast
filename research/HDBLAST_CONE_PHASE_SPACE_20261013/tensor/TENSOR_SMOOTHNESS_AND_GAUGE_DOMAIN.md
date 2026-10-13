# Tensor regularity, energy completion, and cone gauge charges

13 October 2026. Independent continuation calculation on the frozen regular
shell spectral checkpoint. This note modifies no repository or frozen file.
It recovers a known tensor horizon diagnostic and applies it to the
nonconstant Einstein–scalar background. It is not a new stability theorem.

## 1. Precise conclusions

1. A nontrivial four-dimensional Einstein graviton lifted with the even
   tensor radial zero profile `h(y)=1` has gauge-invariant five-dimensional
   Weyl components `rho^-2 delta C[gamma]`. It cannot extend as a C2 metric
   perturbation through the cone vertex whenever its slice Weyl tensor is
   nonzero. Radial normalizability therefore does not imply smooth extension.
2. This does **not** remove the massless mode from the radial energy
   completion. Cutting the profile off near the vertex changes both its
   squared radial norm and radial quadratic form by `O(epsilon^3)`.
   At a regular inertial initial slice the corresponding local H1 × L2
   data can also be approximated, with energy cost `O(epsilon^2)`.
   Smooth configurations and their energy completion are distinct objects.
3. The smooth translation residual in the Hessian-degenerate scalar sector
   is a bona fide smooth Cartesian vector. The other Gaussian-normal
   integration constant generally is not C1 at the vertex. This classifies
   the smooth local gauge generators more precisely than their formal
   expressions in singular radial coordinates.
4. Smooth gauge generators and perturbations have vanishing charge on
   shrinking fixed-de-Sitter-time vertex cuts. This is a small-area estimate,
   not a no-charge or no-flux result at a finite cut of the null cone.

These statements preserve the conditional spectral theorem. They disallow
equating its Hilbert space with a pointwise smooth configuration space, while
leaving equality of suitable completed physical spaces unresolved.

## 2. Cone geometry: vertex and finite null cuts

Assume, more explicitly than `rho ~ y`, that the regular background has

    rho(y) = y + O(y^3),       phi(y) = u + O(y^2),

with smooth even coefficient functions in `y^2`. The polynomial potential
and regular analytic background used by the checkpoint satisfy this local
form. Define a conformal radius `r` by

    dr/r = dy/rho(y),        r/y -> 1.

Then `rho = r Omega(r^2)` and `dy = Omega(r^2) dr`, with `Omega(0)=1`.
For global unit dS4 coordinates,

    gamma = -dt^2 + cosh^2(t) dOmega_3^2,
    X^0 = r sinh(t),         X^i = r cosh(t) n^i,

the metric becomes exactly

    g = Omega(s)^2 eta,      s = eta_AB X^A X^B = r^2.

The scalar is a smooth function of `s`. A local smooth continuation to
negative `s` exists; its global extent is not fixed by this local statement.
The radial patch is the spacelike exterior of the null cone.

Set `u_c = r exp(-t)` and `v_c = r exp(t)`. Then

    X^0=(v_c-u_c)/2,   R=(v_c+u_c)/2,   s=u_c v_c,
    eta = du_c dv_c + ((u_c+v_c)/2)^2 dOmega_3^2.

At fixed finite `t`, `y -> 0` reaches only `X^A=0`. A finite future null
cut has `u_c -> 0`, `v_c > 0` fixed, requiring `t -> +infinity`.
The past cone requires the opposite correlated limit. No fixed-time
radial trace estimate establishes uniform control in these limits.

## 3. Exact tensor Weyl identity on the nonconstant scalar background

Let `q_mu nu` be a transverse-traceless solution of the four-dimensional
linearized Einstein equation

    delta R_mu nu[gamma; q] = 3 q_mu nu.

This is the physical zero-mass slice equation in the checkpoint's tensor
reduction. Consider

    delta g_yy = delta g_y mu = 0,
    delta g_mu nu = rho(y)^2 q_mu nu(x),      delta phi = 0.

Its radial profile is exactly `h=1`: it solves `h''+4 a h'=0`, has the
even Neumann junction, and has squared radial norm `integral rho^2 dy <
infinity`. It is not admitted by the odd Dirichlet junction.

The following identity does not require `p=phi'` to vanish. For any
four-dimensional Einstein metric `gamma_hat` satisfying
`Ric[gamma_hat]=3 gamma_hat`, the warped curvature formulas are

    R^5_mu nu alpha beta
      = rho^2 R^4_mu nu alpha beta
        - rho^2 rho'^2 (gamma_mu alpha gamma_nu beta
                       - gamma_mu beta gamma_nu alpha),
    R^5_y mu y nu = -rho rho'' gamma_mu nu,
    Ric^5_yy = -4 rho''/rho,
    Ric^5_mu nu = (3-rho rho''-3 rho'^2) gamma_mu nu.

Substitution in the five-dimensional Weyl decomposition gives

    C^5_mu nu alpha beta = rho^2 C^4_mu nu alpha beta,
    C^5_y mu y nu = C^5_y mu nu alpha = 0.

Linearizing this identity at dS4 is legitimate for `q` obeying the
linearized Einstein equation; no exact nonlinear slice family is needed.
For an orthonormal tangential frame it gives

    delta C^5_hat a hat b hat c hat d
         = rho^-2 delta C^4_hat a hat b hat c hat d.             (1)

The background Weyl tensor vanishes, so `delta C^5` is invariant under
linearized diffeomorphisms: `delta_xi C^5 = L_xi C_background = 0`.
Along a fixed finite-dS-coordinate ray the displayed orthonormal frame
has a finite Cartesian limit. If any component of `delta C^4` at that
slice point is nonzero, (1) diverges like `y^-2` in a bounded frame.
A C2 metric perturbation of the smooth extended background has continuous,
bounded linearized curvature. Therefore this perturbation has no C2
extension at the vertex, in any smooth gauge.

This uses the Weyl tensor, not the scalar `C_ABCD C^ABCD`. The latter is
quadratic around a conformally flat background and can vanish for nonzero
Lorentzian radiative Weyl tensors. A vanishing invariant scalar would not
evade the componentwise obstruction.

An elementary local example uses the flat dS chart and a constant spatial
trace-free polarization `e_ij`:

    gamma = -dt^2 + exp(2t) delta_ij dx^i dx^j,
    q_ij = exp(2t) f(t) e_ij,       f(t)=exp(-3t).

The tensor equation is `f''+3f'=0`, and the linearized orthonormal Weyl
component for a diagonal polarization is

    delta C^4_hat 0 hat i hat 0 hat i = -(3/2) f(t) e_ii.

This example is a local curvature witness, not a globally finite-energy
plane wave. Globally smooth finite-energy slice examples are obtained from
TT tensor harmonics on the compact S3 sections of global dS4. For a harmonic
of degree `ell >= 2`, its time amplitude obeys

    f'' + 3 tanh(t) f' + ell(ell+2) sech^2(t) f = 0.

Smooth non-gauge initial data have nonzero Weyl tensor and finite energy on
each finite-time S3 slice. Their warped `h=1` lifts satisfy the same
obstruction.

The case `delta C^4=0` is deliberately excluded from this proposition.
Locally on de Sitter, linearized Einstein deformations with zero linearized
Weyl are locally gauge (up to global qualifications). The argument cannot
turn such directions into physical singular gravitons.

## 4. Why this obstruction does not exclude the energy-completed zero mode

The tensor radial Hilbert norm and factorized form, written in `h`, are

    N[h] = integral_0^T rho^2 |h|^2 dy,
    Q[h] = integral_0^T rho^4 |h'|^2 dy.                         (2)

Choose a fixed smooth cutoff `chi(s)` equal to zero for `s <= 1` and one
for `s >= 2`, and set `h_epsilon(y)=chi(y/epsilon)`, with `2epsilon<T`.
It equals one near the shell, so it retains the even Neumann condition.
From `rho/y -> 1`,

    N[h_epsilon-1] <= C integral_0^(2epsilon) y^2 dy
                      = O(epsilon^3),
    Q[h_epsilon-1] <= C epsilon^-2
                       integral_epsilon^(2epsilon) y^4 dy
                      = O(epsilon^3).                         (3)

The profiles have smooth initial data at the vertex because they vanish in
a neighborhood of it. They are initial-data approximants, not separated
zero-mass eigenfunctions: their subsequent evolution excites other radial
masses. A smooth core approximating an eigenvector need not consist of
eigenvectors at the same eigenvalue.

For a direct inertial check, take the regular slice `X^0=0`, corresponding
to `t=0`. A smooth angular TT perturbation with `h=1` has bounded Cartesian
spatial metric components away from `r=0`, homogeneous of degree zero to
leading order. Its spatial derivatives and its inertial velocity are
`O(r^-1)`, since `partial_X0 = r^-1 partial_t` there. In four spatial
dimensions the local H1 × L2 energy is bounded by

    integral_0^epsilon r^3 r^-2 dr = O(epsilon^2).

The derivative of the cutoff adds `O(epsilon^-1)` on its annulus, at energy
cost `O(epsilon^2)` as well. Multiplication by smooth positive `Omega`
does not change these local estimates. For global S3 tensor harmonics the
spatial metric perturbation and its momentum have only angular components,
with zero angular divergence and trace. Multiplication by a radial cutoff
preserves these properties: the possible radial divergence term is
proportional to the angular trace. On the time-symmetric background slice
`X^0=0`, the spatial background Ricci tensor is isotropic in the angular
directions and the scalar gradient is radial. Consequently these data
satisfy the linearized Hamiltonian and momentum constraints, including
after the radial cutoff. A detailed global completed phase-space and
evolution theorem with the shell is still separate.

More explicitly, `delta K` is angular TT and `delta phi=delta pi_phi=0`.
The momentum constraint is its spatial divergence, which vanishes. The
Hamiltonian curvature variation
`delta R=-Ric_ij delta h^ij + div div(delta h)-Delta tr(delta h)` vanishes
because the last two terms are zero and the angular Ricci contraction is
proportional to the zero angular trace. The radial scalar-gradient term is
unchanged. These are linear constraints, not nonlinear gluing claims.

Thus the curvature obstruction distinguishes classical smooth configurations
from weak finite-energy states. It does not prove that the radial completion
is too large as a finite-energy phase space. It supplies no reason to delete
the normalizable massless tensor eigenvalue, and proves no equality of
radial and covariant completions either.

## 5. Hessian-degenerate harmonics: the regular gauge vector

Only the actual five-dimensional harmonic subspace

    D_mu D_nu Y = -gamma_mu nu Y

is at issue here. Write `Y=c_A n^A`, where `X^A=r n^A` embeds the unit
de Sitter hyperboloid. These harmonics have `Box Y=-4Y`; other solutions
at the same mass value are not included.

The checkpoint's complete local Gaussian-normal residual has

    xi^y=e0 Y,       xi^mu=L(y) D^mu Y,
    L=ell0-e0 F,     F'=rho^-2,
    delta phi=p e0 Y,
    bar psi=(a+F)e0-ell0.

Fix the additive constant of `F` by its odd Laurent expansion at the
vertex. Let `A(s)=-r F` and `B(s)=[Omega(s)^-1-A(s)]/s`. The equation for
`F` becomes

    A-2s A_s=Omega^-1,       A(0)=1.

An exact regular choice, beyond a formal expansion, is

    F(r)=-1/r + integral_0^r [Omega(v^2)^-1-1]/v^2 dv.

The integrand is smooth and even, so `-r F(r)` is smooth in `r^2`; its
constant term is one. Thus `B` is smooth as well.

With `Omega=1+omega s+nu s^2+...`,

    A=1+omega s+(nu-omega^2)s^2/3+...,
    B=-2omega+4(omega^2-nu)s/3+... .

The translation part is therefore the smooth Cartesian vector

    xi^A_e0 = e0 [A(s)c^A+B(s)(c_B X^B)X^A].                   (4)

In particular `xi^A_e0 -> e0 c^A`. What appears as `L~e0/y` in radial
coordinates is a smooth translation-type vector, with smooth curvature
corrections. The scalar perturbation `p e0Y` is also smooth because
`p=O(y)`.

The remaining independent vector is

    xi^A_ell0 = ell0 [r c^A-(c_B X^B)X^A/r].                  (5)

It need not even be continuous in the full ambient approach to the vertex.
For example take `c=(0,1,0,...)`, `X^1=epsilon`,
`X^0=sqrt(epsilon^2-epsilon^6)` and the other components zero, with
`0<epsilon<1`. Then `s=epsilon^6`, `r=epsilon^3`, and
`xi^1_ell0=ell0(epsilon^3-epsilon^-1)`, which diverges if `ell0!=0`.
Along bounded finite-de-Sitter-time rays the vector is O(r), but its
first derivatives already have direction dependence. Either test excludes
a C1 Cartesian vector. A smooth gauge convention therefore fixes `ell0=0`
in this normalization. The
background equation gives the useful additional cancellation

    (a+F)'=-p^2/3,
    a+F=-(1/3) integral_0^y p(v)^2 dv = O(y^3),

where the zero integration constant follows from the regular expansion.

Equation (5) is **not** declared a charged physical degree of freedom.
A weaker completed gauge group might contain it or its cutoff limits.
Its charge and admissibility would have to be established in that group's
specified phase space. The smooth local classification alone cannot decide
this. Likewise the constant scalar harmonic has its separate static Jacobi
and shell determinant argument; it is not covered by (4), and the entire
mass values zero or minus four must not be discarded.

If smooth bulk diffeomorphisms are gauge and transform the shell embedding
at the same time, (4) is an admissible local bulk gauge generator. The
checkpoint's two-copy scalar matching still excludes the invariant relative
shell displacements through `p != 0` and `B_boundary > 0`. This note adds
no proof that every smooth perturbation admits the generic longitudinal
representative.

## 6. What can be proved about gauge charges at the vertex

For the action density `R/2-(nabla phi)^2/2-U`, use the Lee–Wald potential
and the field-independent-vector Noether form

    theta^a = sqrt(-g) [ (nabla_b h^ab - nabla^a h)/2
                         - nabla^a phi chi ],
    Q_xi = -(1/2) *d(xi_flat),
    k_xi(delta)=delta Q_xi - i_xi theta(delta).

Signs depend on the common orientation, without affecting the vanishing
estimate. Suppose `g,phi,delta g,delta phi,xi` and all the derivatives
entering `k_xi` are bounded in a smooth Cartesian chart near the vertex.
For the cuts `S_epsilon={y=epsilon,t=t0}`, at a fixed finite `t0`,

    area(S_epsilon)=O(epsilon^3),
    |integral_S_epsilon k_xi(delta)| <= C(t0) epsilon^3 -> 0.   (6)

This includes the smooth translation generator (4). No apex charge
survives under these bounded-jet hypotheses. On a complete closed Cauchy
surface, ordinary smooth diffeomorphisms compatible with the doubled shell
and its embedding are gauge; the shell contributions must be included in
the full presymplectic identity. If a region is cut at the cone instead,
boundary conditions and possible surface degrees of freedom become part of
the definition of that different phase space.

At fixed nonzero `v_c`, the null cut has finite radius `v_c/2` and finite
area. The constant `C(t0)` in (6) was never uniform as `t0 -> infinity`.
The small-area argument is inapplicable. A vector nonzero on that cut is
not automatically a zero-charge gauge direction for the isolated exterior
region. A sufficient genuinely zero-charge class at a physical boundary
consists of generators whose value and first jet vanish there, or whose
support avoids it. These statements do not establish a finite-cone charge
classification for arbitrary radiative data.

## 7. Relation to prior work and remaining gate

The singular zero-mode diagnostic is already explicit in J. Garriga and
M. Sasaki, *Brane-world creation and black holes*, Phys. Rev. D 62, 043523
(2000), arXiv:hep-th/9912118v2, section III, p. 9. They give
`C5^2=C4^2/a^4` and explain that the zero mode is singular on the Rindler
horizon, while a mode expansion is a tool for constructing Green functions.
Their caution is central here. The present calculation recovers that known
mechanism and verifies its componentwise applicability to the registered
nonconstant-scalar warp; it makes no priority claim. The cutoff estimates
make explicit why smoothness failure cannot by itself remove the mode from
an energy completion. Lee–Wald and Iyer–Wald supply the standard covariant
phase-space and charge framework used in section 6.

The physical interpretation gate now has a sharper form: specify whether
the desired space is smooth extendible solutions, their finite-energy
completion on a regular Cauchy surface, or an exterior phase space with
independently prescribed characteristic data; then prove the reduction and
gauge quotient identify that space. The reduced self-adjoint theorem solves
none of these identifications by its limit-point statement alone.

The companion script checks curvature coefficients, a nonzero slice Weyl
witness, the exact power laws in the cutoff estimates, conformal null
coordinates, and the regular-vector expansion. It does not mechanically
prove PDE well-posedness, global gauge completeness, or functional analysis.
