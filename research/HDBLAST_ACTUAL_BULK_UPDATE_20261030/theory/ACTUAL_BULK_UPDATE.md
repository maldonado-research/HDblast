# A boundary-preserving update using the actual reduced bulk operator

30 October 2026. This note advances the preceding auxiliary reconstruction
from its frozen bulk forcing to one application of the retained actual bulk
operator and actual affine bulk source. The result is an integral-defined
iterate, with a conditional C1 majorant. It is not a completed second-order
Einstein/constraint solve, numerical C2 certificate, nonlinear family or
physical observable. The iteration and contraction argument are classical;
no new iteration method or originality claim is made.

## Fixed data, operators and domain

All fields use the same registered positive Gram norm H and retained cone

```
D+ = {0<=x<=3/2, r>=max(0,x-5/4), x+r<=7/4},
X(V)=max(sup_D+ ||V||_H/3, sup_D+ ||V_x||_H, sup_D+ ||V_r||_H).
```

Let `PN` and `PD` be the Robin and Dirichlet projectors. The auxiliary
principal-wave propagator uses the mixed reflection `Pi=PN-PD`. The map `K`
is the zero-initial-data wave response with homogeneous Dirichlet/Neumann
conditions, and `L` is the inward-sign characteristic boundary lift,

```
(Lg)(x,r)=-integral_0^max(x-r,0) g(s) ds.
```

The physical bulk source is completed on `r>=0` before its **output** is
reflected with `Pi`. No bulk coefficient matrix is commuted through `Pi`,
`K`, or a boundary resolvent.

Let `Bs` be the principal free wave with the existing smoothed nominal initial
data. Let `Ba` be the principal free wave for the retained actual initial data.
The target U is the same retained formal response as in the earlier
certificates. Its actual coefficient functions and sources are fixed:

```
La Y = A(x,r) Y_x+B(x,r) Y_r+D(x,r)Y,
M=PN M PN,  PD G=0,
U = Ba+K(La U+F)+L(M PN U+G).
```

The actual shell matrix M is constant in shell time, lies within its inherited
coefficient boxes, and need not equal the earlier frozen matrix M0. The same
fixed actual M, G, A, B, D and F are used in every iterate and in U. Interval
enclosures may discard correlations for bounds, but independent coefficient
choices at every occurrence would not define this iteration. The lower-order
operator is linear in the second-order response Y; the affine source F is
fixed by the retained first-field problem. This proof does not apply unchanged
to a nonlinear operator whose coefficients or sources depend on Y.

Compatible fields, sources, traces and errors retain the inherited zero-wedge
and corner premises. In particular fields vanish on `x+r<=1`, and the active
shell history has length at most `1/4`. The already established estimates are

```
sup_D+ ||V||_H <= (3/4) X(V),
sup_shell ||PN V||_H <= (1/4) X(V),
X(Kf) <= (3/2) sup_D+ ||f||_H,
X(Lg) <= sup_shell ||g||_H.
```

These estimates and the retained first-field/background/source interpretations
remain conditional premises. A source magnitude cap alone does not establish
zero support, corner compatibility, differentiability, or its physical origin.

## Close the actual boundary law before updating the bulk

Define operators on the compatible C1 space by

```
QY=K(La Y),
BY=L(M PN Y|shell),
RB=(I-B)^(-1),
Jactual=RB Q.
```

For an arbitrary base `S` having the mixed homogeneous auxiliary boundary
conditions, `RB(S+LG)` is exactly the Robin reconstruction from the preceding
checkpoint. Its trace solves `Z=PN S|shell-integral(MZ+G)`. Thus it preserves
homogeneous Dirichlet data and satisfies the actual affine law

```
PN V_r=M PN V+G.
```

The old candidate W0 has frozen forcing `f0=L0 Bs`, where
`L0Y=A0Y_x+B0matY_r+D0Y`. Its exact integral identity is

```
W0=RB[Bs+K f0+LG].
```

Here `B0mat` is a matrix, distinct from the free wave Bs. The new affine map
and first update are

```
Psi(Y)=RB[Bs+K(La Y+F)+LG],
W1=Psi(W0).
```

Consequently W1 uses the actual nonconstant bulk coefficients and actual
affine bulk source, evaluated on W0. Its Cauchy data are still those of Bs,
and its shell law is exactly the same actual affine Robin/Dirichlet law.

Set

```
r0=La W0+F-f0
  =(La-L0)Bs+La(W0-Bs)+F.
```

Subtracting the two exact integral identities gives

```
C0=W1-W0=RB K r0,                         (1)
C0=K r0+L(M PN C0|shell).                 (2)
```

The affine shell source G cancels in (2) only because both fields use the
same actual source and same actual M. The correction has zero initial data,
homogeneous Dirichlet trace and homogeneous actual Robin law

```
PD C0=0,    PN (C0)_r=M PN C0.
```

With sufficient inherited regularity the principal-wave identities are

```
(W1)_xx-(W1)_rr=La W0+F,
(C0)_xx-(C0)_rr=r0.
```

These put the actual reduced bulk operator into the forcing once. They do not
give `(W1)_xx-(W1)_rr=La W1+F`. With the same forcing-minus-principal convention
as r0, the next registered forcing defect is

```
r1=La W1+F-[La W0+F]=La(W1-W0)=La C0.    (3)
```

The opposite principal-minus-actual-RHS convention would instead give
`-r1=La(W0-W1)`. The fixture below verifies that this defect need not vanish.
No Einstein, gauge-subsidiary or junction-constraint residual is declared zero
by replacing the principal source once.

## Boundary-preserving contraction and correction bound

Let `na,nb,nd,nm` be uniform H operator-norm upper bounds for the actual
coefficients, including their difference from the frozen coefficients. Set

```
k=(3/2)(na+nb+(3/4)nd),
b=nm/4,
q=k+b,
beta=k/(1-b).
```

Then `X(QV)<=k X(V)` and `X(BV)<=b X(V)`. If `q<1`, then `b<1` and

```
RB=sum_(j>=0) B^j,
||RB||_X <= 1/(1-b),
||Jactual||_X <= beta < 1.
```

The boundary Volterra reconstruction itself exists more generally than this
small-norm bound; the present verified smallness is used to obtain simple
rational estimates. No inverse of the shell matrix M is used. In particular,
`RB K La` is not replaced by `K La RB`.

Let `R_X` be a bound on `X(K r0)`. It is a response-norm bound, not necessarily
a pointwise or derivative bound on r0. Equation (1) gives

```
X(C0)<=Dcap=R_X/(1-b).                     (4)
```

The existing direct residual ledger supplies such an `R_X`: its lower-order
residual H bound is multiplied by `3/2`, and the affine bulk response is
bounded as `X(KF)<=SOURCE+Gcap`, because SOURCE bounds `X(KF+LG)`.
This response number must not be relabeled as `sup||F||_H` or multiplied by an
extra `3/2` as though it were a raw source cap.

All later updates are `W_(n+1)=Psi(W_n)`. Linearity of the fixed actual bulk
operator gives

```
W_(n+2)-W_(n+1)=Jactual[W_(n+1)-W_n],
X(W_(n+1)-W_n)<=beta^n Dcap.               (5)
```

For instance, the component correction at `(9/8,0)` has magnitude at most
`Dcap/(8 sqrt(H_i))`, while generic first-derivative corrections are at most
`Dcap/sqrt(H_i)`. Robin rows have the sharper boundary relation

```
|(C0)_r,i| <= (Dcap/8)
              sqrt(sum_j max|M_ij|^2/H_j).
```

These are norm-derived enclosures of an integral-defined correction. They do
not evaluate its pointwise sign or establish that a norm bound is attained.

Equation (3) supplies a useful further corollary without requiring individual
second derivatives of W1 to be bounded. The registered forcing difference is
well-defined from the C1 field and its specified principal forcing. On the
whole retained cone,

```
sup ||r1||_H <= (na+nb+(3/4)nd) Dcap.       (6)
```

When interpreted as a classical principal-wave defect, the regularity
conditions below are still required. It is not the complete pointwise
Einstein tensor or subsidiary residual.

For the authenticated model caps in the companion bulk-bound calculation,
the right side of (6) is approximately `3.5739437385060643e22`, strictly below
`3.574e22`. A separate complete old-sample residual enclosure gives
`||r0(9/8,0)||_H>9e22`. Consequently the new **global** forcing-defect upper
bound is below the old sample's lower bound divided by 2.5. Under the same
retained premises, this establishes a reduction of the registered reduced
forcing defect, rather than only an improvement in a response-error upper
certificate. It does not establish a reduction of the unknown true field
error, a numerical C2 error, or a physical effect. The source and norm
conventions must match the companion reports for this comparison to apply.

## A continuous local triangle majorant

For the selected shell event set `a=9/8`, `d=a-1`, and use the complete past
physical triangle

```
Omega_a={(t,r):t>=0,r>=0,t+r<=a}.
```

This triangle is inside the retained cone. Suppose a physical source f
vanishes on `t+r<=1` and has uniform H bound `Fcap` on Omega_a. At an arbitrary
output `(t,r)` in that triangle, the original wave diamond has source time
`0<=tau<=t` and spatial range
`r-t+tau<=y<=r+t-tau`. After folding negative y to `v=|y|`, the direct and
reflected physical intervals are respectively

```
[max(0,r-t+tau), r+t-tau],
[0,max(t-r-tau,0)].
```

An interval with nonpositive length contributes nothing. Both lie in
`tau+v<=t+r<=a`. Any physical point `(tau,v)` has at most two preimages,
the original y values `+v` and `-v`. The Duhamel factor is one half, and
reflection by Pi is an H isometry. Therefore its absolute integration weight
after the triangle inequality is at most one. The active part of the common
physical triangle has area

```
Area=integral_1^a alpha d alpha=(a^2-1)/2,
```

where `alpha=tau+v`. Consequently, **throughout** Omega_a,

```
||Kf||_H <= Area Fcap,
||(Kf)_t||_H,||(Kf)_r||_H <= a Fcap.        (7)
```

The derivative estimates follow from the two characteristic edge integrals:
each has time length at most `t<=a` and carries factor one half. The upper
edge has `tau+v=t+r`, so when the output is active it can remain active for
the full time history. Replacing a by d in the derivative bound, or replacing
Area by `d^2/2`, would be invalid. At the even shell output `(a,0)`, the
comparison integrand `1_{tau+v>1}` attains Area and has outgoing derivative
a; with odd output reflection its normal derivative is a. This comparison
integrand is a bound fixture, not an asserted discontinuous model source.

For example, suppose local initial-profile bounds (including convolution
support) are `U0,U1,V` for displacement, its derivative and velocity. Positivity
of the free-wave coefficients and the zero collar give

```
||Bs||_H <= U0+d V,
||(Bs)_t||_H,||(Bs)_r||_H <= U1+V.
```

The velocity integral has at most two folded copies of `[1,a]` with factor
one half. The same distinction from the full characteristic edge length
applies to this value estimate. A uniform frozen-source cap is then
`f0cap=(na0+nb0)(U1+V)+nd0(U0+d V)`.

Let `hcap` bound `PN(Bs+Kf0)` over the shell history; the above estimates give
`hcap=U0+d V+Area f0cap`. The total boundary load of W0 is
`g=M PN W0+G`. Its trace equation implies
`g=Mh-MJg+G`. Thus, for `m>=||M||` and `m*d<1`,

```
gcap=(m hcap+Gcap)/(1-m d),
||W0-Bs||_H <= Area f0cap+d gcap,
||(W0-Bs)_t||_H,||(W0-Bs)_r||_H <= a f0cap+gcap.             (8)
```

All bounds are uniform on Omega_a, not just evaluations at its shell vertex.
With actual coefficient-error caps and an authenticated **pointwise** local
bound on the genuine affine F, these estimates give a uniform local residual
cap Rcap for r0. In that case the correction identity (2) and (7) yield

```
Cvalue=Area Rcap/(1-m d),
sup_Omega_a ||C0||_H <= Cvalue,
sup_Omega_a ||(C0)_t||_H,||(C0)_r||_H <= a Rcap+m Cvalue.
```

At the shell vertex its more precise component structure is

```
PD C0=PD (C0)_t=0,
||PN (C0)_t||_H <= a Rcap+m Cvalue,
PN (C0)_r=M PN C0,
||PD (C0)_r||_H <= a Rcap.                 (9)
```

If a pointwise F bound is unavailable, the already established response bound
`X(KF)<=SOURCE+Gcap` may instead be kept as a **separate** contribution:
`||KF||_H<=d(SOURCE+Gcap)` on Omega_a and its first derivatives are bounded by
`SOURCE+Gcap`. One may then apply (7) only to the local lower-order residual.
The response bound cannot be inserted into Rcap as a pointwise F bound.

Equations (7)–(9) give conservative local integral majorants. They are not a
validated quadrature of the actual bulk correction or an independent
evaluation of F. They can improve sample enclosures when the local residual
caps are appreciably smaller than their whole-cone counterparts.

## Initial-data bias and the meaning of the limit

The map Psi is a contraction on the compatible affine C1 space. Its limit
`Us` solves the retained reduced integral evolution with the **smoothed
nominal** initial data, actual coefficient functions, actual affine F, actual
shell M and actual G. Its initial data differ from the target U's actual
initial prescription. Smoothing need not preserve the original Einstein
initial constraints. Thus Us is an auxiliary reduced-evolution limit, not a
newly established physical second-order or nonlinear spacetime solution.

The target identity can be written

```
U=RB[Ba+KF+LG]+Jactual U.
```

Let `einitial>=X(Ba-Bs)`, the existing actual-minus-smoothed-free-wave initial
bound. Subtracting the iterate identity gives the exact biased error relation

```
U-W_(n+1)=RB(Ba-Bs)+Jactual(U-W_n).        (10)
```

Therefore a valid recurrence for upper certificates is

```
E_(n+1)=alpha+beta E_n,
alpha=einitial/(1-b).                     (11)
```

The initial term cannot be dropped merely because W1 has zero correction
Cauchy data. W1 and W0 share initial data; U and W0 generally do not.

The distinct limiting comparison and the geometric tail satisfy

```
X(U-Us)<=einitial/(1-q),
X(Us-W_n)<=beta^n Dcap/(1-beta),
X(U-W_n)<=einitial/(1-q)+beta^n Dcap/(1-beta).
```

The first term is a limiting **upper majorant**, not a lower bound on the
true error and not a proven unavoidable physical error floor. Both true
differences can be smaller than the available bounds.

When the initial certificate is computed consistently as
`E0=(einitial+R_X)/(1-q)`, these formulas are algebraically equivalent to

```
E_n=(einitial+beta^n R_X)/(1-q).           (12)
```

In particular, one actual-bulk update contracts only the bulk residual part;
it retains the actual initial-data uncertainty. Alternatively, if using the
looser inherited total cap `qret>=q`, use `kret=qret-b>=k` and
`betaret=kret/(1-b)` consistently throughout (11)–(12). This yields the
conservative family
`E_n,ret=(einitial+betaret^n R_X)/(1-qret)`. An arbitrary readable E0 can also
be propagated by (11), but is not asserted equal to the exact closed family.

## Regularity and finite verification

The contraction proof is a C1 integral-map statement. It needs compatible
zero-wedge/corner hypotheses, continuous bounded actual coefficients and
forcing in the stated function spaces. The reflected source is formed after
the complete physical coupling. Physical-side C2 principal-wave identities
for a finite update additionally need sufficient initial regularity, actual
coefficients and F of class C1, the preceding W0 of class C2, and G of class
C1 with the matching corner jets. These are separate inherited premises when
used. The C1 contraction alone neither proves C2 regularity of its limit nor
provides a quantitative C2 error or differentiated constraint estimate.

`check_bulk_update_identities.py` uses exact rational polynomial algebra and
imports no numerical HDBLAST evaluator. It constructs a known full
variable-coefficient polynomial target with a distinct actual initial datum,
then derives nonzero affine bulk and shell sources from that target. A
nonzero nilpotent Robin matrix makes its boundary resolvent finite without a
matrix inverse. The actual pointwise bulk matrices are dense and do not
commute with parity or the Robin matrix. The forcing fixture is globally
parity-compatible so a separate exact polynomial principal-wave recurrence
can be compared with physical reflected triangle integration. This is a
controlled finite fixture, not an assumption of commuting coefficients in
the model proof.

The checks verify the actual affine boundary law, homogeneous correction law,
principal source update, exact residual decomposition, nonzero remaining
bulk defect, contact jets through order two, and the initial-biased error
identity. Controls reject wrong reflection order, frozen-M substitution,
omission of F/G, a reversed lift sign, freezing the actual bulk coefficients,
and omission of the initial bias. Separate rational fixtures verify the
majorant recurrence and closed form through six iterates. Exact polygon/edge
integration of the active-source comparison envelope checks the local
geometry at ninety physical output/sign fixtures and rejects the tempting
short-history value and derivative caps. All 39 named groups give identical
normal and optimized Python output.

These finite algebra checks support the displayed analytic identities. They
do not numerically resolve the HDBLAST affine bulk source, re-prove its
first-field certificates, or supply external specialist/proof-assistant
verification. Publication-level physical and nonlinear claims remain outside
this result.
