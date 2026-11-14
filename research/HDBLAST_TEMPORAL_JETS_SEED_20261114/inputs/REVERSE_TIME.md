# Quantitative same-state past evolution with the full tensor reflection

13 November 2026 UTC. Reviewed conditional derivation. Earlier frozen source
is preserved. This closes the orientation join between the positive-time
first/second wave estimates and the original state's past proper-observer
problem. All conditional analytic existence, state-flow, support, weighted
inverse, gauge and uniqueness premises of the retained theorems remain
in force. No numerical nonlinear amplitude radius is supplied.

Let `ell=10^-6`, `T=13ell/8`, `Z=6ell`. For the original selected
initial state, including its nonzero velocity, the result is

```
sup_(-T<=t<=0) sum_(j=0)^14
       ||partial_t^j V1(t)||H^(14-j)
       <B1=10^4000000 exp(10^152017),

sup_(-T<=t<=0) sum_(j=0)^13
       ||partial_t^j V2(t)||H^(13-j)
       <B2=10^10200000 exp(10^152019).                    (1)
```

The parameter convention is unchanged:
`V_e=e V1+e^2 V2/2+o(e^2)`. The parameter e is not reversed.
The actual past first wave is not identified with the actual future first
wave. A reflected elapsed-time problem with different velocity data is
used to prove (1), then pulled back to the original coordinates.

## 1. Exact source and domain joins

All following paths are relative to the research directory of the checkout
at merged PR137, commit `7fb54e385e5250deae408910da55bf000a53f5a2`.
The needed sources are:

* `HDBLAST_HIGH_ORDER_FIRST_WAVE_20261112/HIGH_ORDER.md` and
  its pinned `inputs/COEFFICIENTS_06G.md`,
  `inputs/BOUNDARY_FLOW_PROOF_06F.md`,
  `inputs/INITIAL_TRANSPORT_06H.md`, and
  `inputs/EXACT_FINITE_REGULARITY_SHELL_EVOLUTION.md`;
* `HDBLAST_SECOND_WAVE_EVOLUTION_20261112/MIXED_SOURCE.md` and
  `SECOND_WAVE.md`, including the actual inhomogeneous initial-jet and
  normal-recovery recurrences;
* `HDBLAST_NONLINEAR_OBSERVABLE_BUDGETS_20261105/inputs/NONLINEAR_BRIDGE.md`,
  sections 1–3, for the precise conditional full-family past application;
* `HDBLAST_FIRST_PROPER_GEOMETRY_20261106/FIRST_CUT_BOUND.md`,
  section 1, for the original nonzero-velocity state and its warning
  against reversing every coordinate metric velocity indiscriminately.

The raw coefficient theorem states the real physical domain explicitly:

```
0<=z<=6*10^-6,  |t|<=13*10^-6/8,  ||V||star<=10^-3.       (2)
```

It therefore supplies the negative-time coefficient values and all the
required derivatives. This is not inferred from a forward numerical
trajectory. The background is

```
g0=rho(z)^2[dz^2-dt^2+cosh(t)^2 kappa],
phi0=phi0(z),  129<=rho<=131.
```

The originally registered exact-family binding stated a slab ending at
`3ell/2`. It is not silently enlarged. NONLINEAR_BRIDGE, section 2,
explicitly applies the retained prescribed-slab theorem to the reflected
exact initial curve on the longer slab with

```
Lsupport=5ell/2, c*=2, Z=6ell, T=13ell/8,
Lsupport+c*T=23ell/4<6ell,
T+Z=61ell/8<8ell.                                       (3)
```

The causal margin is ell/4. This full-family application retains the
same state-flow, inverse, support, regularity and subsidiary-system
hypotheses on that slab. The arithmetic is not itself a proof of those
premises or an evaluated parameter radius. The new contribution here is
the complete quantitative transport of the already proved jet estimates
to that actual past development.

## 2. Tensor time reflection on the fixed sixteen-field bundle

Let `chi(s,z,Omega)=(-s,z,Omega)`, with `0<=s<=T`. On the
tangential coframe `(dt,cosh(t)theta^A)`, its pullback is represented by

```
S=diag(-1,1,1,1),     S^2=I,
S J4 S=J4,            J4=diag(-1,1,1,1).
```

The equality `cosh(-s)=cosh(s)` is essential. On

```
E=R_v + Sym^2(R dt+T*S3)_H + R_w + (R dt+T*S3)_D
```

define the constant fiber reflection

```
J(v,H,w,D)=(v, S H S, w, S D).                            (4)
```

Thus H_tt and H_AB are even, H_tA is odd, D_t is odd and D_A
is even. The root and scalar slots are even. There are four odd
independent fields, not sixteen odd fields. J is an involution and an
isometry of the positive direct-sum/Frobenius norm and of the star
norm. It preserves the shell Dirichlet/Neumann summands separately.

The fixed product connection acts only on the angular slots, so J is
parallel. It commutes with spatial covariant derivatives, radial
derivatives, natural angular rotation lifts and both shell projectors.
The base pullback supplies the additional sign on each time derivative
slot. These signs are orthogonal transformations of all positive tensor
norms used in the estimates.

For any physical metric and scalar, the actual tensor pullback is

```
g^-_ss(s)=g_tt(-s),
g^-_sI(s)=-g_tI(-s),
g^-_IJ(s)=g_IJ(-s),
phi^-(s)=phi(-s),                                        (5)
```

where I,J are spatial directions, including the radial direction.
The background is fixed by (5). The outward spatial shell normal and
the inward coordinate z are not reversed.

## 3. The exact boundary chart commutes with this reflection

This must be checked for the actual point map rather than assumed from
the tensor principal symbol. In the pinned proper-y normalization let
`H=J4+Htilde`, `n=1+v`, `eta=eta_b+w`. The exact frozen generator is

```
B_v=-rho_b{(2/3)[sigma(eta)n^2-sigma_b]
                  -a_b[n^3 tr(H^-1 J4)-4]},
B_H=-rho_b[sigma(eta)n H-sigma_b J4]/3,
B_w= rho_b[n sigma'(eta)-sigma'_b]/2,
B_D=0.                                                  (6)
```

Since `(S H S)^-1=S H^-1 S` and `S J4 S=J4`,

```
tr[(S H S)^-1 J4]=tr(H^-1 J4).
```

The scalar rows in (6) are unchanged and the H row is conjugated by S.
The zero D row is compatible with D -> S D. Therefore

```
B(JU)=J B(U),
F_z(JU)=J F_z(U),    F_-z(JU)=J F_-z(U),                  (7)
```

on the admitted flow domain, by uniqueness of the finite-dimensional
flow. Its radial parameter z is not reflected.

The reconstruction has
`G_zz=rho^2 n^2`,
`G_ab=rho^2 J4_ab+rho_b^2 Htilde_ab`,
`G_za=-rho rho_b D_a`, and `phi=phi0+w`.
Together with (7) and the even angular coframe, these identities show
that the transformed unknown of the pulled-back physical fields is
exactly

```
V^-(s,z,Omega)=J V(-s,z,Omega).                           (8)
```

The physical normalized metric patch is also preserved: the full
five-dimensional reflection is an orthogonal diagonal sign matrix
preserving J5, so `||Gamma-J5||op` is unchanged. No additional
state-domain constant appears.

## 4. Reflected Cauchy data and the actual conditional past family

At the canonical initial section the background momenta vanish.
Time reflection of exact data is

```
(gamma,pi,phi,varpi) -> (gamma,-pi,phi,-varpi).             (9)
```

The Hamiltonian constraint is unchanged and the momentum constraint
changes its overall sign. Exact constraint satisfaction, regularity,
reflection symmetry, radial support and the initial open background
collars are preserved. This is the reflection of the entire selected
curve, including its second and higher completion, not only of its
first tangent. In the radial weighted construction the metric correction
potential is retained and the momentum correction potential changes
sign. Thus the reflected curve remains in the same weighted class.
The retained localized theorem permits the reflected velocity profile;
its hypotheses do not require a favorable sign of that profile.

In ADM language the reflected initial lapse is unchanged and the shift
is negated. The selected shift is zero. The exact relative-harmonic
completion is natural under the tensor pullback, because the full
reference g0 is invariant. Equivalently one can pull back its completed
metric velocities directly. Equation (8) gives, at every amplitude,

```
V^-(0)=J V(0),
partial_s V^-(0)=-J partial_t V(0).                       (10)
```

For the amplitude derivatives this becomes

```
V1^-(0)=J V1(0),  (V1^-)_s(0)=-J (V1)_t(0),
V2^-(0)=J V2(0),  (V2^-)_s(0)=-J (V2)_t(0).              (11)
```

Even-time-slot velocities reverse; odd-time-slot velocities do not.
For example the actual first completed mixed velocity is

```
(g_yt,1)_t=-6rho^2(u_y+4a u)Y,
```

where u in this line is the registered radial profile. It is generally
nonzero. Since g_yt has one time slot, (5) implies
`(g^-_ys,1)_s(0)=(g_yt,1)_t(0)`, not its negative. A rule negating
every metric component's time derivative would fail this join.

Apply the retained prescribed-slab theorem to this exact reflected curve
and its complete reflected gauge data, on `0<=s<=T` with (3). It
provides a conditional exact reflected family and C2 parameter dependence
into X10 on some common, unevaluated parameter neighborhood. Pulling it
back by the inverse reflection produces an exact original-coordinate
family on `-T<=t<=0` whose position and velocity at t=0 are precisely
the original selected data:

```
J V^-(0)=V(0),
-J partial_s V^-(0)=partial_t V(0).
```

Uniqueness identifies it with the original state's past development.
It may be joined to the original future branch after taking the smaller
of their two qualitative parameter neighborhoods. This gives the actual
conditional full-family past slab needed by a proper-observer argument,
not merely an independent formal linear wave. Its numerical amplitude
radius is not inferred from the at-zero bounds below.

## 5. Reflected background operator and unchanged coefficient budgets

Write the actual background energy-sign operator in fixed product
directions as

```
Lw=beta^2[w_tt-w_zz-sech(t)^2 Delta_E w]
                       +B^t w_t+B^i D_i w+A w,
beta=1/rho.
```

Here i is spatial and B is the complete fixed-bundle first-jet map.
For `Rw(s)=Jw(-s)`, direct substitution gives

```
L^- Rw=R(Lw),

L^-U=beta^2[U_ss-U_zz-sech(s)^2 Delta_E U]
                       +B_-^s U_s+B_-^i D_i U+A_-U,
A_-(s)=J A(-s)J,
B_-^s(s)=-J B^t(-s)J,
B_-^i(s)= J B^i(-s)J.                                   (12)
```

This identity already suffices for all estimates. Geometric covariance
also identifies these coefficients as the same natural reduced equation
written in the reflected tensor coordinates; no freezing of the varying
reference Christoffel symbols is involved.

Every mixed coefficient derivative in (12) is the corresponding old
derivative at -s, conjugated by constant isometries and multiplied by
signs. The state-derivative input slots transform in the same way. The
old coefficient domain (2) includes all these points and is invariant
under J. Hence all entry-sum and full tensor operator bounds through the
orders used in the first/second wave proofs are unchanged. In particular

```
C_raw=10^150000, Dstar=10^152000,
Hcoef=10^152300, R_init=10^152310, M_norm=10^152320.         (13)
```

The normalized lower maps are still multiplied by the same rho^2,
whose time derivatives vanish. The natural angular generators commute
with J, so their exact operator commutators, rank estimates and Sobolev
comparisons are unchanged.

## 6. Forced quadratic source and initial temporal recurrences

Let `u=V1`, `v=V2`, and `u^-=Ru`, `v^-=Rv`. The original energy-sign
second force is

```
F2[u]=2Q1[u]:H u-C2(u,u)-2E2[u]Du-2D2(Du,Du).            (14)
```

On a first-jet direction i, time reflection contributes a sign r_i,
with r_t=-1 and r_i=1 spatially. On a Hessian there is one such sign
for each base derivative slot. Thus, for instance,

```
C2^-[U,U]=J C2(-s)[JU,JU],
Q1^-[U]^{ij}=r_i r_j Q1(-s)[JU]^{ij}.
```

The analogous first-jet coefficient maps contain the corresponding r_i
and fiber J in each input slot. Substitution into (14), or differentiating
the full reflected nonlinear equation twice in the unchanged amplitude,
gives the exact identity

```
F2^-[u^-](s)=J F2[u](-s),
L^-v^-=F2^-[u^-].                                       (15)
```

There is no overall minus sign on the source and no change of any
factor two. The Hessian's two time signs and the output fiber reflection
are already included. Since all these maps are isometries, the complete
mixed source product bounds retain exactly their old constants.

For the initial temporal recurrences let `G=rho^2 F2` and
`G^-=rho^2 F2^-`. The source jets satisfy

```
partial_s^j G^-(0)=(-1)^j J partial_t^j G(0).              (16)
```

The forced time equation retains the positive source sign,

```
(v^-)_ss=(v^-)_zz+sech(s)^2 Delta_E v^-
           -P_-^s(v^-)_s-P_-^i D_i v^- -Q_-v^-+G^-.
```

Its normal recovery retains the negative source sign,

```
(v^-)_zz=(v^-)_ss-sech(s)^2 Delta_E v^-
           +P_-^s(v^-)_s+P_-^i D_i v^- +Q_-v^- -G^-.
                                                               (17)
```

To verify every initial step explicitly, set
`J_j=partial_t^j v(0)`, `J_j^-=partial_s^j v^-(0)`.
The claimed relation is `J_j^-=(-1)^j J J_j`.
The k-th time derivative of the reflected temporal lower coefficient
P_-^s has sign `(-1)^(k+1)`, while the spatial and zeroth-order
ones have sign `(-1)^k`. The corresponding factors in the initial
recurrence are J_(j-k+1) or J_(j-k), respectively. Each product
therefore has the common sign `(-1)^j` and fiber J. The same is true
of the principal terms and, critically, the nonzero source jet (16).
This proves the relation inductively, starting from (11). The homogeneous
first-wave recursion has the same argument with G=0.

All corner jets still vanish on their initial open collars. The complete
source initial hierarchy is reflected, not omitted. In particular
the norms of the first/second initial data, the recursively recovered
initial temporal jets, and the source jets at zero are identical to
their original-coordinate counterparts.

## 7. Reapply the energy argument in elapsed time

One must not reverse the direction of a one-sided differential inequality
and keep its sign. Instead apply the proved energy identity to L^- in
the increasing variable s. Its scalar principal coefficients are

```
a=beta^2, b_principal=0, c=beta^2 I.
```

The moving angular frame is `X_A^-=sech(s)E_A`, and its exact identity
is

```
partial_s D_(X_A^-) U=D_(X_A^-)U_s-tanh(s)D_(X_A^-)U.
```

Thus the signed moving-frame term in the elapsed-time energy derivation
has exactly the old form. The inverse metric positivity, radial divergence
budget, full lower-map norms and absolute time/coefficient derivative
budgets are unchanged. The same evaluated constants apply:

```
Kstar=10^150006,  C_comm=10^152020,
Lambda=Kstar/2+95*C_comm<10^152022,
Lambda*T<10^152017.                                     (18)
```

The shell has the same four Dirichlet and twelve Neumann rows, because
J preserves that split and commutes with D_z. The outer Dirichlet wall
is also preserved. Consequently the boundary flux cancellations remain
exact. No boundary row is imposed on a recovered normal derivative.

At each fixed s, all spatial Sobolev norms are isometric under J. For
each time derivative,

```
partial_s^j (Rw)(s)=(-1)^j J partial_t^j w(-s),
||partial_s^j (Rw)(s)||H^m=||partial_t^j w(-s)||H^m.       (19)
```

The same equality holds for all mixed radial/time/natural-rotation word
norms. Therefore the full initial, energy, source and normal ledgers of
the predecessors apply without an extra component, reflection or frame
factor. Their relevant numerical interfaces are

| Interface | Unchanged ceiling |
| --- | --- |
| First initial H14/H13 sum | `10^5200` |
| Second initial H14/H13 sum | `10^10000` |
| First full mixed X14 sum | `B1=10^4000000 exp(10^152017)` |
| Evolving second-force mixed sum S12 | `10^8200000 exp(10^152018)` |
| Initial second-force sum I11 | `10^4200000` |
| Normalized radial/time/rotation source Z11 | `10^250 S12` |
| Second full mixed X13 sum | `B2=10^10200000 exp(10^152019)` |

For the first full mixed sum, use the explicitly proved W14 comparison
in MIXED_SOURCE, rather than infer it from the displayed N14 norm alone.
All these estimates are uniform over compatible data with the stated
initial budgets; they do not require the original positive velocity sign.
Applying them to the elapsed-time solutions, then using (19), proves (1).
There is no claim that the two opposite sides of the original solution
have equal field values, only that both admit the same upper bounds.

## 8. Actual parameter derivatives and proper-observer output

The reflected full family is C2 into the retained lower X10 topology.
Its first derivative solves the reflected homogeneous linear problem;
its second derivative solves (15) with the actual reflected second
data. The stronger first X14 and second X13 linear solutions are
identified with those actual derivatives by the same energy uniqueness
argument as in PR136/137. Pullback then identifies the bounded fields
with the actual original-coordinate derivatives on the negative slab.
This does not assert C2 dependence of the entire nonlinear curve into
X13 or X14.

For the original past proper observers, the background path is still

```
q=-tau/rho_b,  X0(q,n)=(-q,n).
```

The quantitative inputs in (1) are therefore available on that original
path, without substituting a different future solution. In particular
`q*=25ell/16<T` has the strict coordinate margin

```
T-q*=ell/16>0.
```

The retained qualitative C2 family and its regular background permit
an unevaluated common smaller parameter neighborhood for a compact
proper-observer tube inside this margin. Proving numerical observer
variational constants at zero can now use the actual past shell traces.
The margin alone does not give a numerical amplitude, uniform metric
patch, global finite-amplitude observer domain or normalized-observable
remainder. Those claims are not part of this reverse-time join.

The proof uses the actual full tensor time pullback, the exact boundary
chart, the two-sided coefficient certificate, the correctly reflected
canonical curve and harmonic velocities, and the inhomogeneous source
recurrence. It preserves the original nonzero-velocity state and all
retained analytic qualifications.
