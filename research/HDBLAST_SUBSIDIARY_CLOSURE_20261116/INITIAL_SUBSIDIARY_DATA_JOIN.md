# Exact initial subsidiary data for the constructed reduced solution

16 November 2026 UTC. New mathematical join for the selected
canonical family. This proves zero initial relative-harmonic subsidiary
value and first normal/time derivative at the actual finite-regularity
reduced solution. It uses the exact physical reconstruction/reduction
interface specified below. It proves no bulk subsidiary propagation or
physical boundary uniqueness; those are separate companion arguments.
All preceding source files remain immutable.

## 1. Precise inputs and conclusion

Fix one admitted original parameter epsilon in the selected canonical
family. Use its unchanged lapse binding, density normalization, Gamma
action and the same-state elapsed past reflection. On
`Sigma=[0,6ell] x S3`, `ell=10^-6`, the already constructed reduced
solution U has

```
partial_s^j U in C([0,T];H^(14-j)), j=0,...,14, T=13ell/8,
U_j(0)=J_j,           j=0,...,14.
```

The interval, point-state patch and high-norm bounds are those of the
accepted reduced theorem. No interval is changed in this join. Let V be
the original-t state, so

```
U(s)=J_E V(-s),
V(0)=J_E U(0), V_t(0)=-J_E U_s(0), V_tt(0)=J_E U_ss(0).
```

J_E is the fixed parallel tensor time-slot reflection, with square
identity. It is not a change of the selected canonical parameter.

The physical reconstruction companion supplies the following exact
interface. The literal inverse point chart, with all background,
proper-y root, shell scales and moving-frame factors, gives
`(g,phi)=Psi(t,x,V)` in the retained original physical coordinates.
In particular g,phi have strong initial time jets in H14,H13,H12,
the first two jets are precisely the actual canonical harmonic Cauchy
completion in Section2 below, and the reconstructed equations are

```
E_AB=Ric_AB-phi_A phi_B-(2/3)mathcalU(phi)g_AB
    =nabla_(A C_B),
S=Box_g phi-mathcalU'(phi)=0,
C^A=g^BC(Gamma[g]^A_BC-Gamma[g0]^A_BC).                 (1)
```

Parentheses include one half. The potential is the retained
`mathcalU(phi)=Uhat(phi-1)`. All reference-connection terms are
included. Equation (1) holds continuously in H12, hence also at t=0;
its initial acceleration is the reconstruction of the actual strong
J2, not separately chosen physical data. This physical interface is
the exact inverse-chart/reduced-equation join delegated to its companion,
not an assumption that the unreduced Einstein equations already hold.

The future unit normal to the initial slice is
`u=alpha^(-1)(partial_t-beta)`. It is distinct from the outward
spacelike shell normal. Our conclusion, in physical covector components,
is

```
C|_(t=0)=0 in H13,
(nabla_u C)|_(t=0)=0 in H12,
partial_t C|_(t=0)=0 in H12.                            (2)
```

The same zero statements hold after lowering/raising indices and after
the declared elapsed tensor pullback; in particular the elapsed
coordinate-covector subsidiary has zero value and zero s derivative.
They apply per copy and on the retained identical-two-copy union.

The constraints required here belong to the exact selected canonical
datum. They are not asserted for generic small analytic Cauchy pairs
or for its unconstrained smooth approximants.

## 2. The actual half-action canonical conversion and first Cauchy jet

Write `s_gamma=sqrt(det gamma)`, `rho0=sqrt(det gamma0)`,
`j_det=det gamma/det gamma0`, `lambda=j_det^(-1/2)=rho0/s_gamma`.
The exact canonical family has

```
gamma=gamma0+h, phi=phi0+psi,
pi_act=rho0 tau, varpi_act=rho0 omega,
alpha=rho(1+epsilon f)>0,       beta^i=0,
f=u_profile(y)(n1^2-n2^2),      dy=-rho dz.
```

The action is `R/2-|dphi|^2/2-mathcalU(phi)`. Its canonical
definitions and inverse in four spatial dimensions are

```
pi_act^ij=(s_gamma/2)(K^ij-k gamma^ij),
varpi_act=s_gamma Pi,
K=(1/2)Lie_u gamma, k=tr_gamma K, Pi=u(phi),
Ttau=gamma_ij tau^ij,
K_ij=2lambda[gamma_ik gamma_jl tau^kl-Ttau gamma_ij/3],
k=-(2/3)lambda Ttau,               Pi=lambda omega.     (3)
```

The doubled-action inverse printed in the older generic gauge note
cannot be substituted for (3). The factors two in the metric momentum
and its inverse, and the scalar normalization in (3), are literal.

The actual physical initial values are

```
g_tt=-alpha^2, g_ti=0, g_ij=gamma_ij, phi=phi0+psi.       (4)
```

Put `R^A=g^BC Gamma[g0]^A_BC`, evaluated at this perturbed initial
metric, and `Gamma_gamma^i=gamma^jk Gamma[gamma]^i_jk`.
The completed velocities are

```
gamma_t=2alpha K,             phi_t=alpha Pi,
alpha_t=alpha^2 k,
beta_t^i=alpha^2(Gamma_gamma^i-R^i)-alpha gamma^ij d_j alpha. (5)
```

Here R^0=0 on the whole selected initial family: at t=0 the reference
coefficients `GammaBar^0_ij,GammaBar^0_00` vanish, and every initial
`g^0i` is zero. The other R^i are evaluated at the perturbed metric;
they have not been frozen to their background contractions.

Equations (3)--(5) give the source's literal component velocities

```
g_ij,t=4rho(1+epsilon f)lambda
                  [gamma_ik gamma_jl tau^kl-Ttau gamma_ij/3],
g_tt,t=(4/3)alpha^3 lambda Ttau,
phi_t=alpha lambda omega,
g_ti,t=gamma_ij beta_t^j
      =alpha^2 Ccov_i(gamma)+gamma_ij rho(D0)^j rho-alpha d_i alpha,
Ccov_i(gamma)=gamma^jk[(D0)_j h_ki-(1/2)(D0)_i h_jk].    (6)
```

D0 is the spatial gamma0 connection. In particular the four mixed
velocities in (6) are present and need not vanish. The equality follows
from the full reference contraction and
`GammaBar^i_tt=rho(D0)^i rho`. The lapse alpha is temporal; it is
not the radial component root `sqrt(g_yy)` in the inverse point chart.

The nonlinear initial-chart construction and its inverse are mutually
inverse on the actual patch. Their chain rule identifies (4),(6) with
the reconstructed value and velocity of the reduced solution. In point
map notation, the full identities are

```
Xphys(0)=Psi(0,V(0)),
Xphys_t(0)=Psi_t(0,V(0))+Psi_V(0,V(0))V_t(0),
Xphys_tt(0)=Psi_tt+2Psi_tV[V_t]+Psi_VV[V_t,V_t]+Psi_V V_tt,
```

with all second-line/third-line map coefficients evaluated at `(0,V(0))`.
The actual even-time frame has vanishing explicit first derivative at
zero for every state; the general chain formula keeps every term
regardless. The last line is determined by J2 and the reduced PDE.
We neither prescribe nor modify it in order to make (2) hold.

## 3. All five initial relative-harmonic equations

The contracted-Christoffel identity
`Gamma[g]^A=-|g|^(-1/2)partial_B(|g|^(1/2)g^AB)` gives in ADM variables

```
Gamma[g]^0=k/alpha
       -(partial_t-beta^j partial_j)alpha/alpha^3,
Gamma[g]^i=-beta^i Gamma[g]^0
       -alpha^(-2)(partial_t-beta^j partial_j)beta^i
       +Gamma_gamma^i-alpha^(-1)gamma^ij d_j alpha.      (7)
```

These use K with the plus Lie-normal sign in (3). On the actual initial
slice beta=0, so substituting (5) gives separately

```
Gamma[g]^0=k/alpha-alpha^2 k/alpha^3=0=R^0,
Gamma[g]^i=-[Gamma_gamma^i-R^i-alpha^(-1)gamma^ij d_j alpha]
             +Gamma_gamma^i-alpha^(-1)gamma^ij d_j alpha
           =R^i,                 i=1,...,4.           (8)
```

Thus all five contravariant components C^A are zero. Multiplication by
the actual metric gives C_A=0 as well. These equalities hold in H13:
every connection term contains one metric derivative, and inverse-metric
and lapse operations are valid on the positive initial patch.

They depend only on the strong original physical value and first
velocity just identified. They remain equalities for the actual finite
datum of the nonlinear limit. No condition on the gauge residual of an
individual analytic approximant has been used. Spatial covariant
tangential derivatives of the zero covector trace therefore vanish in
H12, including their connection terms, which multiply C itself.

## 4. Exact canonical constraints are the required geometric projections

Define the trace-adjusted residual

```
H_AB=E_AB-(1/2)g_AB tr_g E
    =G_AB-[phi_A phi_B-(1/2)|dphi|_g^2 g_AB-mathcalU(phi)g_AB]. (9)
```

This is the Einstein-minus-scalar-stress residual for the actual
half-action convention. The Gauss--Codazzi identities with
`K=+Lie_u gamma/2`, the retained curvature convention and u timelike
unit are

```
G(u,u)=[R(gamma)+k^2-|K|_gamma^2]/2,
Ric(u,partial_i)=nabla^gamma_j K^j_i-partial_i k.         (10)
```

For example in Gaussian normal coordinates at the slice,
`Gamma^0_ij=K_ij`, `Gamma^i_0j=K^i_j`. In Ric(u,i), their
contraction gives exactly the second identity of (10). The contraction
forming G(u,u) cancels the normal derivative of k and gives the first.
These are tensorial slice identities; this verification changes neither
the actual chosen lapse/shift nor its completed gauge velocities.

Since `|dphi|_g^2=-Pi^2+|Dphi|_gamma^2`, (9),(10) give

```
H(u,u)=[R(gamma)+k^2-|K|^2-Pi^2-|Dphi|^2-2mathcalU]/2,
H(u,partial_i)=nabla^gamma_j K^j_i-partial_i k-Pi partial_i phi. (11)
```

From (3), direct full-tensor contraction gives

```
(k^2-|K|^2)/2
 =2lambda^2[(Ttau)^2/3-tr(gamma tau gamma tau)],
Pi^2=lambda^2 omega^2.                                 (12)
```

Substitution into (11) is precisely the source's full Hamiltonian
constraint, with `Lambda=lambda^2`. Consequently
`H_constraint=H(u,u)`, with no extra factor two and no discarded metric
or momentum component.

For clarity the density divergence and its sign are also explicit.
The physical covariant derivative on a weight-one symmetric density
satisfies

```
nabla^(gamma,dens)_j pi_act^jk
   =s_gamma nabla^gamma_j(pi_act/s_gamma)^jk,
```

because the physical volume density is parallel. Therefore the actual
canonical momentum density obeys

```
(M_density)_i/s_gamma
 =-2gamma_ik nabla^gamma_j(pi_act/s_gamma)^jk+Pi partial_i phi
 =-nabla^gamma_j K^j_i+partial_i k+Pi partial_i phi
 =-H(u,partial_i).                                     (13)
```

This also reproduces the exact fixed-rho0 source rather than replacing
rho0 by s_gamma. Set
`A^k_ij=Gamma[gamma]^k_ij-Gamma[gamma0]^k_ij`. The two density trace
actions cancel in the contracted divergence, so

```
nabla^(gamma,dens)_j(rho0 tau^jk)
 =rho0[(D0)_j tau^jk+A^k_jl tau^jl],
-2gamma_ik A^k_jl tau^jl
 =-2tau^jl(D0)_j h_li+tau^jl(D0)_i h_jl.
```

Thus

```
(M_density/rho0)_i
 =-2gamma_ik(D0)_j tau^jk-2tau^jk(D0)_j h_ki
    +tau^jk(D0)_i h_jk+omega partial_i phi,
lambda(M_density/rho0)_i=-H(u,partial_i),               (14)
```

exactly the current nonlinear constraint map. The source proves its
Hamiltonian and all four momentum components zero for the selected
canonical family. Equations (11)--(14) therefore give

```
H(u,u)=0,             H(u,partial_i)=0,   i=1,...,4.     (15)
```

The identities lie in H12 with these finite inputs; spatial derivatives
of K and the spatial Ricci tensor have that regularity. The induced
metric, scalar, K and Pi in (11) are the ones of the constructed metric
by its strong physical Cauchy value and velocity. The projections in
(15) do not depend on an independently assigned second normal metric
jet. Gauss--Codazzi in (10) explains their cancellation of that jet.

## 5. The reduced initial acceleration forces the subsidiary velocity to zero

The strong reconstructed J2 and (1) make every term of
`E_AB=nabla_(A C_B)` continuous in H12 up to t=0. We can thus evaluate
this equation at that endpoint using its actual one-sided time trace.
Let `X_A=(nabla_u C)_A|0`. The trace C|0 is zero by (8), and
`nabla_(e_i) C|0=0` for every spatial tangent by Section3. In a
g-orthonormal frame `(u,e_i)` this gives

```
E(u,u)=X_u,   E(u,e_i)=X_i/2,   E(e_i,e_j)=0,
tr_g E=-X_u,
H(u,u)=X_u/2,                 H(u,e_i)=X_i/2.           (16)
```

The minus sign in the trace is the timelike unit sign. The one half in
the mixed component is the declared symmetrization. Orthogonalizing
the spatial frame is an invertible algebraic operation on the positive
initial metric; (15) therefore supplies exactly the zero projections
in (16). It follows that X_u=0 and all X_i=0, so

```
(nabla_u C)|0=0 in H12.                                (17)
```

This is the point at which the actual reduced acceleration is used.
Canonical constraint zero alone would not identify a freely specified
normal derivative of an unrelated covector field. Here C is the gauge
residual of the same reconstructed metric, and (1) ties its derivative
to the actual reduced metric equation. No extra condition beyond the
completed velocities (5),(6) is imposed on g_tt,t, and no independent
g_tt,tt is prescribed.

Finally `partial_t=alpha u+beta^i partial_i`. The covariant derivative
along partial_t of C is zero by (17) and its zero spatial tangential
derivatives. The difference between this covariant derivative and the
ordinary derivative of a coordinate covector is a connection times C,
which is zero initially. Hence

```
partial_t C_A|0=0 in H12.                               (18)
```

Raising the index differentiates the metric, but its extra differentiated
factor also multiplies C=0. Thus the contravariant gauge vector has the
same zero time derivative. This proves all original-t data in (2).

## 6. Finite regularity, elapsed zero data and exact scope

The operations above need only the actual physical jets through order
two. The metric/scalar values H14, first velocities H13 and second
velocities H12 give C|0 in H13 and its derivative in H12. Smooth fixed
coordinate/frame conversions and uniformly invertible metric/lapse
factors preserve those orders. The physical equations and Gauss--Codazzi
involve at most two spacetime derivatives; these high Sobolev orders on
the four-dimensional spatial slice justify their classical identities
and their strong H12 endpoint traces. Equivalently finite-order smooth
approximation proves these local identities, followed by the stated
Sobolev product limits. Such an algebraic identity check does not assign
physical constraints to the approximants used to construct the solution.

Let iota(s,y,Omega)=(-s,y,Omega). The elapsed subsidiary covector is
the natural tensor pullback `C_hat=iota^* C`, with the reflected physical
and reference metrics. The spatial components are C_i(-s), and the
time-covector component is `-C_t(-s)`. Equation (2) gives

```
C_hat|s0=0 in H13,        partial_s C_hat|s0=0 in H12.    (19)
```

The time-slot reflection is constant, so there is no additional frame
derivative. The elapsed time unit normal has the corresponding opposite
original future orientation; its derivative of this zero covector is
zero either way. In the subsidiary energy companion's fixed product
covector control bundle one may use C_y,C_a. The alternative normal
component C_z=-rho C_y has equivalent zero initial value and derivative
because rho is time-independent. Any bounded Robin straightening of
that control field preserves (19); no conjugation is needed to state
the zero-data interface.

All constructions are natural under the retained Gamma action. The
current canonical family, lapse profile and initial completion are the
same invariant data used by the reduced theorem. The proof is per-copy
and hence gives zero subsidiary data on the specified identical copies
with their global Hilbert norms, without a new copy factor. The future
slice normal used throughout is not either copy's outward shell normal.

The conclusion is exact zero subsidiary initial value and first
normal/time derivative for the constructed reduced solution at each
admitted canonical parameter. No boundary Robin identity, subsidiary
bulk uniqueness, later-time harmonic propagation, Israel/force claim,
evolved parameter regularity or observer result is inferred here. Those
must use their separate exact reconstructed equations and boundaries.

## Read-only source basis and companion join

At immutable base `6ee27c36dac67497a2454f7c1283fd7177a91888`:

- `research/HDBLAST_UNIFORM_CAUCHY_DATA_20261113/UNIFORM_CAUCHY_DATA.md`,
  Sections1--3, SHA256
  `785b20a06ecac643d635a89f70a800a765d20164d22da7e2df3c777a31b3bfa9`:
  actual half-action conversion, selected lapse/shift and all sixteen
  value/velocity fields.
- `research/HDBLAST_NONLINEAR_CONSTRAINT_BALL_20261113/NONLINEAR_CONSTRAINT_BALL.md`,
  Section3 equation(12) and exact selected-family conclusion, SHA256
  `04f35658bd66d2c3cfb9d1b11cbf942d349595e0231dd11b3a8d847241ab1f54`.
- `research/HDBLAST_QUASILINEAR_SHELL_EVOLUTION_20261015/gauge-parameter/`
  `INITIAL_GAUGE_SUBSIDIARY_AND_PARAMETER.md`, Sections1--2, SHA256
  `74bce09080d4015b2eccf8c183b33f47b80516216a9547837e3e0e34251d7b03`:
  geometric gauge identities and the trace-adjusted projection argument;
  its doubled-action canonical inverse is replaced by the actual (3).
- `research/HDBLAST_EVOLUTION_BOUNDS_20261106/COEFFICIENTS.md`, exact
  reconstruction/physical reduction sections, SHA256
  `7b28775e6c44b8ef1213e4789c4efb3c6b82727226cd18c5953c0d20038f9f39`.
- `research/HDBLAST_NONLINEAR_PICARD_LIMIT_20261114/`
  `REDUCED_EVOLUTION_AND_CAUSAL_SUPPORT.md`, SHA256
  `94a860dc021dde62e3110ea8cee1688e3e43ad48a6d04585c2686849dc45637c`:
  actual solution, unchanged admission and strong J0/J1/J2 endpoint traces.

The accepted new companion
`PHYSICAL_RECONSTRUCTION_AND_BOUNDARY_JOIN.md`, accepted original-proof SHA256
`c9b6bde860cca55f581bbf58b29e26e97f2e057f6fe38ad6ce1093e33719d5db`,
applies the literal inverse chart and complete physical reduction to
that solution. Its Sections1,5 give the finite physical regularity and
reduced equations; Section7 gives the exact value/velocity inverse
identity and full J2 chain used in Section2 here. Authentication and
finite replay do not replace these analytic hypotheses or this written
initial-data argument.
