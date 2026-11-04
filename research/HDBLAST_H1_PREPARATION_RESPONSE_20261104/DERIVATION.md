# A direct H1 preparation-to-response estimate

4 November 2026. Conditional coefficient theorem. The exact profile norm radius
`R<=2e31`, tilt, boundary weight and numerical guards were declared before the
final numerical gates. No prior capsule is modified.

## 1. Retained model and the weaker geometric datum

Use the same fixed background, first packet, complete second constraints,
conditional exact-family development, harmonic conditioning, mixed linear
uniqueness, physical shell/anchor and proper-time target as the published
preparation-dependence and intrinsic-preparation-robustness results. The
parameter convention is actual second derivatives: no extra Taylor factor 1/2
is inserted in the response coefficients. All background, operator and nonlinear
analytic premises remain conditional.

The initial area profile difference D is measured between two families sharing
that first tangent, at fixed background scalar label on the specified physical
Cauchy slice and continuing the selected scalar-level S3 component. Its label
derivative is `D_r=(-ell rho p) partial_v D`, with `ell=10^-6` and the fixed
background calibration. The slice, component and parameter normalization remain
part of the datum. Define the complete-increment norm

```
I=[1,T], T=5/4, Delta=T-1=1/4,
R^2=integral_I (D^2+D_r^2) dr.                              (1)
```

For `base+lambda*h`, use `|lambda| R(h)`. This norm is not an identified physical
energy. Nor is its radius selected by a physical action or state.

The retained compact metric-mean h is smooth, flat outside the packet and has
zero weighted moment. With `chi_y=rho^2 h`, `z=chi/(ell rho^3)`,
`b=ell rho a`, `g=ell rho p^2/(3a)`, it obeys

```
z_r=3b z-h,
h_r=3ell^2 z+(-b+g)h-(2/3)gD, z(1)=h(1)=0.               (2)
```

Here `p=phi0_y!=0`, `a=rho_y/rho>0`. The same background enclosures used in the
preceding capsule give, writing `Rinf=sup_I |D|`,

```
|h|<=H Rinf, |h_r|<=Hr Rinf, |V_L|<=V Rinf,
H=131/73426500000000,
Hr=3ell^2 H+(bM+gM)H+(2/3)gM, V=Hr+bM H,
bM=917/62500000, gM=131/16317000000000.                     (3)
```

These bounds need only Rinf, not the supremum of D_r: the forcing in (2)
contains D alone. They follow from the same continuous Gronwall calculation,
whose matrix exponential is bounded by 4/3 on this initial interval. The new
arithmetic reconstructs (3) exactly from the authenticated prior caps; it does
not evaluate a new background or evolved-field oracle.

## 2. Exact operator, normalization, domain and retained caps

In mean order `U=(K,L,N,J,w)`, with the constant positive inner product
`Hgram=diag(1,2,1,3,1)`, the complete homogeneous preparation difference solves

```
U_xx-U_rr=A(x,r)U_x+B(x,r)U_r+Dmat(x,r)U,
L(x,0)=0, PN U_r(x,0)=M(x) PN U(x,0).                    (4)
```

The normal metric L is the Dirichlet slot; the other four slots form PN. These
projections are orthogonal for Hgram. The boundary derivative in (4) is the
increasing-r derivative already converted from `U_r=-ell rho U_y`; changing to
an outward-normal sign without changing M would be an error. The initial
scalar velocity is zero. Every metric/scalar, connection and Robin coupling is
retained in (4); M is not required to be symmetric, positive or constant.

The selected source file
`research/HDBLAST_TOTAL_RESPONSE_PREFLIGHT_20261023/analytic/EXPLICIT_PREMISES.json`
explicitly declares the *conditional continuous* bounds in this Hgram norm:

```
||A||,||B||<=10^-4, ||Dmat||<=10^-6, ||M||<=1/2000.        (5)
```

Its status is `EXPLICIT_CONDITIONAL_PREMISES_NOT_REPROVED_BY_PORTABLE_CHECKER`.
It declares the domain

```
Dplus={0<=x<=3/2, max(0,x-5/4)<=r<=7/4-x}.
```

The new backward cone `D_T={0<=x<=T,0<=r<=T-x}` is contained in Dplus: `T=5/4`
makes the lower radial bound zero and `T-x<=7/4-x`. The same Hgram, scale and
caps were bound by the prepared-family response authentication. The exact
five-mean projection is a retained analytic premise. Thus (5) restricts to
this cone and block. Byte authentication of these assertions is not a fresh
all-point proof of the coefficient bounds, not a source re-evaluation and not
a PDE proof. This distinction is present in the report and source-pins record.

The compatible solution is zero on `x+r<=1` by the retained finite propagation
and initial collar. This is a coefficient-level cone statement, not a common
finite-parameter nonlinear cone.

## 3. Tilted analytic energy and the actual Robin boundary

Write `V=U_x`, `Q=U_r`, and use the Hgram inner product below. Let

```
k=1/2, alpha=1/2, F=A V+B Q+Dmat U,
e=(|V|^2+|Q|^2+2k<V,Q>+|U|^2)/2,
u(x)=PN U(x,0),
E(x)=integral_0^(T-x) e(x,r)dr + alpha |u(x)|^2/2.          (6)
```

The tilt is a mathematical estimating device. It is not the gravitational
Hamiltonian, canonical Noether energy or a new physical positivity law.
Because `0<k<1`,

```
e >= (1-k)(|V|^2+|Q|^2)/2+|U|^2/2.
```

Direct differentiation, using the constant principal wave operator and the
constant Hgram, gives

```
partial_x e=partial_r[<V,Q>+k(|V|^2+|Q|^2)/2]
              +<V+kQ,F>+<U,V>.                            (7)
```

No derivative of A,B,Dmat or M is used. At the moving outer endpoint `r=T-x`,
the flux minus the moving energy density is exactly

```
-(1-k)|V-Q|^2/2-|U|^2/2 <=0.                              (8)
```

At r=0 the Dirichlet slot has `V_L=0` and contributes `-k|Q_L|^2/2<=0`.
On the PN slots put `P=PN(V+Q)`. Since `PN Q=M u`,

```
u_x=P-Mu,
B_N=-k|P|^2/2-(1-k)<P,Mu>+(1-k)|Mu|^2.                  (9)
```

Adding the derivative of `alpha|u|^2/2`, and setting
`m=1/2000`, `beta=alpha+(1-k)m`, gives

```
B_N+alpha<u,u_x>
 <= -k|P|^2/2+beta|P||u|+[(1-k)m^2+alpha m]|u|^2
 <= -k|P|^2/4+[beta^2/k+(1-k)m^2+alpha m]|u|^2.          (10)
```

The last line is Young's inequality. It remains valid for nonsymmetric M,
because only its induced Hgram norm was used. In terms of the boundary
summand of E, the growth coefficient is

```
Cboundary=(2/alpha)[beta^2/k+(1-k)m^2+alpha m]
          =2003001/1000000=2.003001.                      (11)
```

For the bulk source use `a0=b0=10^-4`, `d0=10^-6`. Expanding absolute values in
`<V+kQ,F>+<U,V>` and bounding each product by half the sum of squares gives
coefficients of `|V|^2,|Q|^2,|U|^2` bounded respectively by

```
cV=a0+(b0+k a0)/2+d0/2+1/2,
cQ=k b0+(b0+k a0)/2+k d0/2,
cU=d0/2+k d0/2+1/2.
```

Consequently the bulk growth coefficient is

```
Cbulk=max(2cV/(1-k),2cQ/(1-k),2cU)=1000351/500000=2.000702.
```

Combining (8)--(11) proves the continuous differential inequality
`E_x<=C E`, for the predeclared `C=201/100`. Gronwall therefore gives

```
E(x)<=exp(Cx)E(0)<=13E(0), 0<=x<=T.                        (12)
```

The scalar exponential bound is strictly below13; multiplication by E(0)=0
may give equality. The strict exponential bound is proved by a rational Taylor sum through
order24 at `CT=201/80`, with tail bounded by the first omitted term divided
by `1-(CT)/26`. No exponential oracle or numerical mesh is used. The identity
holds first for `x<T`; the regular solution and vanishing interval length
allow passage to x=T. The positive boundary-state summand remains at that
limit, so the scalar trace is still controlled there.

## 4. H1 profile norm controls all initial energy channels

For an H1 function on length1/4, the elementary one-dimensional inequality
applied to D gives

```
sup|D| <=2||D||_2+(1/2)||D_r||_2
        <=(sqrt(17)/2)R <=(21/10)R.                         (13)
```

The actual retained smooth directions obey this inequality; no higher profile
regularity is newly assumed. In contrast, applying such an inequality to D_r
would require H2, which is not needed here.

After pure second harmonic conditioning the initial difference is exactly

```
U0=(h,0,0,h,w), U0_x=(0,V_L,0,0,0), U0_r=(h_r,0,0,h_r,w_r),
w=kappa(h-2D/3), kappa=p/(2a), V_L=h_r+b h,
w_r=kappa_r(h-2D/3)+kappa(h_r-2D_r/3).                   (14)
```

Here `kappa` is distinct from the energy tilt k. Use the old exact absolute
bounds `kappaM=1/1554`, `kappa_rM=509197/3449880000000`. Since the initial
support in the cone is I and its length is1/4, equations (3),(13),(14) give
L2 caps per R

```
h2=(21/20)H, hr2=(21/20)Hr, v2L=(21/20)V,
w2=kappaM(h2+2/3),
wr2=kappa_rM(h2+2/3)+kappaM(hr2+2/3),
D0=2h2+w2, Dr=2hr2+wr2, Dv=2v2L.                         (15)
```

The 2 factors safely bound the metric Gram norm and the L velocity weight.
Every one of the five channels, the metric mixed velocity, the scalar
pullback and the derivative of kappa is retained. The bound on wr uses
`||D_r||_2<=R`, not its pointwise supremum.

The initial tilted cross term vanishes **exactly**: velocity is in the L
slot, while the radial derivative has zero L slot, and Hgram is diagonal.
The initial boundary value also vanishes by the common collar. Thus

```
E(0)<=C0 R^2, C0=(D0^2+Dr^2+Dv^2)/2 <1/5000000.           (16)
```

This last strict rational inequality is a new exact-arithmetic gate. Combining
(6),(12),(16), for every scalar boundary value `m(x)=w(x,0)` and R>0,

```
|m(x)|^2 <= (2/alpha)E(x) <[52/5000000]R^2 <R^2/90000,
|m(x)| <= R/300.                                         (17)
```

The R=0 case follows directly from E(0)=0 and Gronwall. The weak energy
estimate is used on the existing regular solution; no new weak-solution
existence or uniqueness theorem is asserted.

## 5. Complete proper observable differences

The two families share their first solution and hence every first-only
quadratic bulk, shell, clock, measure and Gauss correction. Those terms
cancel only in their **difference**; the absolute base coefficients retain
them. The exact proper-time Green formulas from the preserved preparation
comparison therefore apply, with `tau*=rho_b ell T` and `F'_b/2`:

```
Delta v2=(rho_b ell)^2(F'_b/2)integral_1^T kernel(T,x)m(x)dx,
Delta Theta2=rho_b ell(F'_b/2)integral_1^T kernel_T(T,x)m(x)dx.
```

The retained exact kernels satisfy `0<=kernel<=T-x`, `0<=kernel_T<=1`.
The zero wedge gives m=0 before x=1. With `rho_b<=131`, `|F'_b|<=1/1000`,
(17) implies

```
|Delta v2|<=Cv R, Cv=(131ell)^2(1/1000)Delta^2/(4*300),
|Delta Theta2|<=Ct R, Ct=131ell(1/1000)Delta/(2*300).        (18)
```

The denominators include the physical forcing factor1/2 and the exact
integrals of the nonnegative kernel upper bounds. No Euclidean-time slope
is substituted for proper-time expansion.

At the predeclared closed radius `R<=2e31`, adding both signs of (18) to the
same authenticated exact negative base intervals gives safe outward gates

```
v2 in [-7.2400e16,-8.6400e15],
Theta2 in [-1.0113e22,-3.3540e21].                         (19)
```

Both upper endpoints are strictly negative. As before, C2 dependence implies
that each admitted exact family separately has smaller positive volume and
slower positive expansion than background at the registered proper time for
sufficiently small nonzero parameter of either sign. No common effective
nonlinear parameter window follows from this coefficient norm bound.

This result bypasses the failed implication from H1 control into the earlier
C1-radius class. It neither shows set containment between the two specified
balls nor claims that any H1-bounded high-frequency example has an unstable
observable response. A physical preparation rule, calibrated scale and a
physical bound implying (1) remain open. No force, radiation generation,
observational conclusion or higher-dimensional Big Bang origin is established.

## 6. Actual admissible preparations beyond the former pointwise criterion

The [actual-image construction](image/PROOF.md) does not choose arbitrary profiles
that might fail the constraints. On the fixed smooth nonsingular background, let

```
psi(u)=u^3 exp(1-1/(1-u^2)) for |u|<1, zero otherwise,
r0=9/8, z_j(r)=j^-5 psi(j^2(r-r0)), integer j>=4,
chi_j=ell rho^3 z_j, h_j=3b z_j-(z_j)_r.
```

Both chi_j and h_j are compact and flat inside the packet. The exact identity
`(chi_j)_r=-ell rho^3 h_j` gives `chi_y=rho^2 h` and the weighted zero moment.
The same retained constraint construction adds only the raw metric mean and fixes
the other canonical sectors. Smoothness of this particular sequence uses the
retained smooth background, separately from the lower C4 parameter topology.

Writing `Aprofile=3/(2g)`, its actual geometric profile is

```
D_j=Aprofile[(z_j)_rr+(-2b-g)(z_j)_r
                    +3(ell^2-b^2+bg-b_r)z_j].
```

For k=0,1,2,3, the exact localization identity is
`||z_j^(k)||_2=j^(2k-6)||psi^(k)||_2`. The finite background coefficient bounds
on the fixed interior interval therefore give a finite C with `||D_j||H1<=C`
for every j. The bump jets at its center are `psi=psi'=psi''=0`, `psi'''=6`.
Consequently, with no asymptotic remainder,

```
(D_j)_r(r0)=9j/g(r0),
(w_j)_r(r0)=-3p(r0)j/[a(r0)g(r0)].
```

The former pointwise radius and its initial derivative-supremum quantity N
therefore diverge. For any budget B>0, fixed multiplication of the entire datum
by beta=B/(2C) places the full sequence inside `R_H1<B`, while those suprema remain
unbounded. Taking B=2e31 gives admitted data controlled by the new theorem that
violate every finite old pointwise cutoff. Each fixed direction has its own
conditional exact-family parameter neighborhood; their higher norms are not
uniformly bounded.

This does not claim containment of either entire specified norm ball in the other.
It proves no unstable or divergent observable: the new theorem supplies the stated
negative coefficients for the scaled sequence. Its analytic bump certificate
checks exact polynomial jets and continuous bounds without numerical bump queries.

## 7. Compatibility and physical selection

The earlier opposite-sign preparation lies outside the new class. The authenticated
old witness at the initial bulk scalar label r=T has
`|lambda D_h(T)|>3.7621e33`. Equation (13) gives

```
R_H1 >= (10/21)|lambda D_h(T)| >1.7914e33 >2e31.
```

This is an initial level label, not the future physical shell event. No new
background or bump evaluation is needed for this comparison.

The [supporting selection note](support/PREPARATION_SELECTION.md) distinguishes
this preparation-sensitive mathematical energy from quadratic Noether tests of
the fixed first tangent. The latter tests, with their structures fixed, are shared
across the old preparation fiber and cannot select a finite preparation budget.
The specific old matched-domain Noether generator has additional hypotheses not
implied by the local five-mean problem. These standard observations do not rule
out nonlinear physical energies: their second variations can contain a linear
term in the second preparation. An explicit physical selection rule remains open.

The [primary-source comparison](REFERENCES.md) identifies the established energy
method and its different boundary hypotheses. The finite checks establish rational
inequalities under retained premises, not a new continuous coefficient proof,
nonlinear PDE existence result, external peer review or observational discovery.
