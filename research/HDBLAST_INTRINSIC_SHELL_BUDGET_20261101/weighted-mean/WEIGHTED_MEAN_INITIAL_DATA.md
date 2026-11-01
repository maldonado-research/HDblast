# Selected weighted family: explicit monopole second Cauchy problem

1 November 2026. Analytic continuation, source commit
`b51908d7e62e7ec80cef82438edd00923ff94723`. Public repository source hashes
and Git objects are recorded in `SOURCE_PINS.json`. This note uses the same actual first packet,
canonical curve and harmonic completion as the exact-family TT transfer.
It does not import the older explicit physical second jet. The analytic
weighted inverse and nonlinear evolution theorems remain conditional inputs.

The scalar initial field and its coordinate velocity are determined:
`phi2=0` and `partial_t phi2=0`. The mean metric and its velocity require
the radial weighted problems below. The cited family specification leaves
the defining function and weight parameters uninstantiated; consequently
these equations do not specify a unique numerical scalar history. Once those
choices are fixed, the theorem's dotted-space inverse is unique; there is
no additional arbitrary inverse to select within that construction.

Nevertheless, the actual current packet forces a strictly negative initial
monopole moment for every such weighted choice. Section 6 proves that sign
using the retained coefficient boxes and an exact integral identity. This is
a gauge-fixed initial constraint result, not a shell energy or expansion
measurement.

## 1. Conventions and authenticated inputs

All coefficients are parameter derivatives, so `X_e=X0+e X1+e^2 X2/2`.
On the initial section, use

```
gamma0=dy^2+rho^2 kappa,  a=rho'/rho,  p=phi0',
a'=-rho^-2-p^2/3,  p'=V'-4ap,
V=p^2/2-6a^2+6rho^-2.
```

Primes here mean outward proper-radial `partial_y`, not scaled inward `r`.
Thus `partial_y=-rho^-1 partial_z`; all integrals below run with increasing
`y`. Write `s=ut/rho`, `j=vt/rho`, with `u,ut,v,vt` precisely as in source
S1 section 2. The normalized first canonical data in the fixed background
orthonormal frame are

```
h1=diag(-4u,2u,2u,2u)Y,
tau1=pi1/rho0=diag(-3s/2,0,0,0)Y,
phi1=vY, varpi1/rho0=jY,
p j=-3(s'+3as).
```

`rho0` is the background volume density, distinct from the radius `rho`.
The last equality is the first momentum constraint; the scalar reconstruction
also gives `pv=-3(u'+2au)`. Sphere averages use probability measure:

```
<Y^2>=1/6,  <|DY|^2>=4/3,  Delta_kappa Y=-8Y.
```

The canonical scalar constraint and normalized momentum constraint are

```
H=R/2+2[(tr pi)^2/3-|pi|^2]/det(gamma)
  -varpi^2/[2det(gamma)]-|dphi|^2/2-V,
Mhat_i=M_i/rho0.
```

This `H` is half the physical constraint in S3. The normalized momentum has
the opposite sign to the physical `div K-d tr K-Pi dphi`, multiplied by
`sqrt(det gamma)/rho0`. These factors matter for the second source.

## 2. Explicit second sources for this canonical affine tangent

Let `QH=<D^2 H[v_seed,v_seed]>` and
`QM=<D^2 Mhat_y[v_seed,v_seed]>`. Then

```
QH = -2u'^2-24a u u'-2u u''
     -(16a'+32a^2+8p^2/3)u^2-4u^2/rho^2
     -s^2-(j^2+v'^2+V''v^2)/6
     -(4/3)p u v'-(4/3)v^2/rho^2,

QM = -4u s'-12a u s-2s u'+j v'/3
   = -2s u'+j(v'+4pu)/3.
```

The full angular canonical momentum Hessian is

```
(D^2 Mhat)_y=(-24u s'-72a u s-12s u'+2j v')Y^2,
(D^2 Mhat)_A=(6us+jv) D_A(Y^2).
```

Its angular mean is zero. To derive the radial expression directly, expand
the exact density identity in S2 with background-covariant derivatives:

```
Mhat_i=-2 gamma_ik D0_j tau^jk-2 tau^jk D0_j gamma_ki
       +tau^jk D0_i gamma_jk+(varpi/rho0) partial_i phi.
```

The quadratic radial terms are respectively
`-24u(s'+3as)Y^2`, `-24su'Y^2`, `+12su'Y^2`, and `+2jv'Y^2`.
For the Hamiltonian, the source S3's physical Hessian can be divided by two
and averaged: the nonlinear canonical conversion of `K,Pi` does not change
that Hessian because the Hamiltonian is quadratic in them at `K0=Pi0=0`.
This reuse concerns its source formula only, not S3's chosen second solution.

For an independent geometric route, set
`L^2=1-4e uY`, `B=rho sqrt(1+2e uY)`, `k=B'/(LB)` and expand

```
R4=B^-2[6-4 Delta log B-2|D log B|^2]
   -12k^2-6L^-1 k'-2L^-1 B^-2[Delta L+D log B.D L].
```

Its second derivative has coefficients
`[-24u'^2-288auu'-24uu''-(192a'+384a^2)u^2-336u^2/rho^2]Y^2`
and `36u^2|DY|^2/rho^2`. Adding the kinetic, gradient and potential terms
gives `QH` above. No second-order scalar profile is inserted.

## 3. Radial weighted adjoint problem

Write the monopole corrections in the fixed frame as

```
<gamma2>=diag(A,B,B,B),
<pi2/rho0>=diag(T,U,U,U),
phi2=varpi2=0.
```

For a radial covariant metric tensor `(A,B)`, the metric-only Hamiltonian
linearization is

```
Pg(A,B)=-3B''/2+(3a/2)A'-6aB'
        +(3a'+6a^2+p^2/2)A-3B/rho^2.
```

For a radial scalar adjoint potential `N`, define

```
Jr[N]=[-3aN'+(3a'+3a^2+p^2)N]/2
     =[-3aN'+3(a^2-rho^-2)N]/2,
Jt[N]=[-N''-2aN'+(a'+3a^2-2rho^-2)N]/2.
```

These are exactly the radial and tangential components of S1's `Pg* N`.
Use `alpha_w` for the integer weight parameter, avoiding confusion with
the background `a`, and put

```
qweight=d^2, w=d^(2alpha_w-4) exp(-beta/d),
WH=w^2 qweight^4=d^(4alpha_w) exp(-2beta/d),
WM=w^2 qweight^2=d^(4alpha_w-4) exp(-2beta/d).
```

For the radial one-form potential `Z=xi dy`, the selected correction obeys

```
A=WH Jr[N], B=WH Jt[N],
Pg(WH Jr[N], WH Jt[N])=-QH,                         (H)

T=2WM xi', U=2WM a xi,
-4[(WM xi')'+3a WM(xi'-a xi)]=-QM.                 (M)
```

Equation (H) is fourth order; equation (M) is second order. Their function
spaces are the radial parts of S2's dotted weighted closure, with `k=12`.
The physical corrections have zero Sobolev extension across the endpoints
of `A'={ell/2<z<5ell/2}`. These are not arbitrary Dirichlet/Neumann boundary
conditions on the potentials. Exponentially singular potential behavior
compatible with that closure must be treated through its weighted weak
formulation. Naive regular-endpoint shooting would need justification.

For compactly supported radial tests `n,chi`, the weak forms are

```
integral rho^3 WH [Jr[N]Jr[n]+3Jt[N]Jt[n]]dy
    =-integral rho^3 QH n dy,

4 integral rho^3 WM [xi'chi'+3a^2 xi chi]dy
    =-integral rho^3 QM chi dy.
```

The full linearized operator and the radial weights commute with rotations.
Although the first quadrupole is not radial, its quadratic source can be
averaged before inversion. There is no radial invariant adjoint kernel by
the retained source classification. Hence these mean equations determine
the mean of the same weighted solution once its weights are fixed.

## 4. Physical Cauchy and initial harmonic conversion

The canonical-to-physical formulas in S1 imply pointwise

```
Pi2=-2uj Y^2,
K2=2[tau2-I tr(tau2)/3]+diag(36us,-6us,-6us,-6us)Y^2.
```

Therefore

```
<Pi2>=-uj/3,
<K2>rr=4T/3-2U+6us,
<K2>angular=-2T/3-us.
```

The selected lapse is `N_e=rho(1+e uY)` and the initial shift vanishes
for every parameter. Since `(sqrt(det gamma))1/rho0=uY`,

```
partial_t phi_e=N_e varpi_e/sqrt(det gamma_e)
              =e rho jY+O(e^3),
phi2=0, partial_t phi2=0.
```

This cancellation holds in every angular sector, not only the mean.
The metric's monopole coordinate velocity is

```
<partial_t gamma2>rr=2rho[4T/3-2U+16us/3],
<partial_t gamma2>angular=2rho[-2T/3-2us/3].
```

Also `<g2_tt>=-rho^2 u^2/3` and `g2_ti=0` initially, because
`g_tt=-N_e^2` and the shift is exactly zero. Lapse/shift time derivatives
must still be obtained from the same relative-harmonic completion in S1;
they are not free second data. These formulas supply the non-gauge spatial
Cauchy part and do not by themselves propagate the monopole through contact.

## 5. What is specified, and exact dependence on the missing choices

S1 section 3 says to select radial defining function and weights as in S2.
S2 permits any admissible radial `d`, a sufficiently small defining-function
rescaling, `alpha_w in N`, and `beta>0`. Neither source gives explicit
values or a formula for those choices. The theorem selects a curve only
after they have been fixed. The weighted complement/inverse is then fixed
by the stated pairing and closure. A numerical mean computation requires
instantiating `d,alpha_w,beta`, the background evaluation and a justified
weighted solver. Replacing those data by convenient choices would define
a newly instantiated family; it cannot be described as recovering a unique
unrecorded old scalar history.

The exact choice dependence is sharper than saying merely “nonunique.” For
either block write `R_W=W P* (P W P*)^-1`, with the appropriate tensor
pairing, and let `b=-Q`. Then `h_W=R_W b`. For a differentiable variation
of admissible weights on a common domain where the inverse differentiates,

```
delta h=(I-R_W P)[(delta log W)h_W],  P(delta h)=0.
```

The differential formula is conditional on that common-domain regularity;
the finite difference identity `P(h_W1-h_W2)=0` needs no such assumption.
Multiplying `W` by a positive constant cancels from `R_W`; only its shape
can affect the physical correction. Thus rescaling `d` while rescaling
`beta` by the same constant changes each displayed multiplier by only a
constant and leaves that second correction unchanged. Varying beta alone
changes `delta log WH=delta log WM=-2 delta beta/d`. These formulas quantify
the possible change but do not prove that every permitted weight variation
changes this particular packet's correction.

## 6. A fixed monopole moment and explicit admissible scalar freedom

There is a useful source-only result even before choosing weights. Direct
substitution gives

```
Pg*(rho)=diag(-3/(2rho),-1/rho,-1/rho,-1/rho).
```

Integrating `Pg(A,B)=-QH` by parts with the zero extension yields

```
integral rho^2(A+2B)dy=(2/3)integral rho^4 QH dy.     (*).
```

This is the areal-radius adjoint moment already present in S4, expressed
in the present coordinates and canonical normalization. It is not a new
gauge-invariant energy. It fixes one moment, not the entire radial profile
or shell history; every pair of weighted choices has zero difference in
this moment.

### Strict sign for the current packet, independent of the missing weights

Let `I=integral rho^4 QH dy`. Compact support and the background equation
give the exact identity

```
I=integral rho^4 [(4/rho^2+8p^2/3)u^2-s^2-j^2/6
                 -(v'+4pu)^2/6-(V''/6+4/(3rho^2))v^2]dy.
```

At the integrand level, the difference from the displayed `QH` is the
derivative of `-2rho^4 u u'-8rho^4 a u^2`, whose endpoint values vanish.
Set `F=rho^2 u`. The first scalar reconstruction gives
`F'=-rho^2 p v/3`. The proper-radial length of the first packet's support
interval is at most `Lmax=131 ell`, because `dy=-rho dz` and its z length
is `ell`. The elementary estimate
`F(y)=integral_(y_-)^y F'`, Cauchy-Schwarz and integration imply

```
integral rho^4 u^2 dy = integral F^2 dy
 <= Lmax^2/2 integral (F')^2 dy
 <= Lmax^2 pmax^2/18 integral rho^4 v^2 dy,
pmax=7/50000, ell=1/1000000.
```

Source S6 supplies `rho in [129,131]`, `|p|<=pmax` and
`eta=phi0-1 in [-9/100000,-1/12500]` throughout this support. The exact
potential has

```
V''=28/9+(100/9)eta-2eta^2-(80/9)eta^3-(20/9)eta^4>3.
```

For `E=9/100000`, subtracting the absolute upper values of all nonconstant
terms already gives a rational lower bound greater than 3. In addition

```
epsilon_pos=(4/129^2+8pmax^2/3) Lmax^2 pmax^2/18 < 1/100.
```

Dropping the negative squares therefore proves the deliberately coarse
but exact bound

```
I <= -(49/100) integral rho^4 v^2 dy < 0,
integral rho^2(A+2B)dy <= -(49/150) integral rho^4 v^2 dy < 0.
```

Strictness follows because the actual `u` is nonzero: if `v` vanished
identically then `F'=0`, and compact support would force `u=0`. This extends
the useful retained moment calculation to the actual current displacement
and velocity packet, with affine canonical scalar sectors. It does not
evaluate the full magnitude, identify a unique `A,B`, or imply a sign for
the anchored shell-volume response after contact.

To exhibit the remaining radial homogeneous data explicitly, take any
smooth compactly supported `f` in the correction annulus satisfying
`integral rho^2 f dy=0`. Set

```
C[f]=Pg(f gamma0)=-3(f''+3af')/2-Vf,
xi_f(y)=-(rho^4 p^2)^-1 integral_(y_-)^y rho^4 C[f] ds,
h_f=f gamma0+Lie_(xi_f partial_y) gamma0,
(h_f)rr=f+2xi_f', (h_f)angular=f+2a xi_f.
```

Indeed

```
Pg(Lie_(xi partial_y)gamma0)=rho^-4 (rho^4 p^2 xi)',
integral rho^4 C[f]dy=-(9/2)integral rho^2 fdy=0.
```

Consequently `xi_f` is compactly supported, `Pg h_f=0`, and the correction
with only this metric component solves all linearized initial constraints.
It has zero angular TT part and preserves both affine scalar sectors. A
nonzero compact radial metric-only perturbation with the scalar held fixed
cannot be a pure radial diffeomorphism of the full data because `p!=0`:
such a diffeomorphism also changes `phi` by `xi p`. This is actual freedom
in completing the first seed, not a relabeling of the specified weighted
range. This infinite-dimensional class includes sign-changing radial
bumps with their one weighted integral canceled.

The same IFT as S2 applied to
`x0+e v_seed+e^2(h_f,0,0,0)/2+A u_tilde(e)` has unchanged first derivative
and unchanged `u_tilde''(0)` because `P h_f=0`. Its second metric datum
differs exactly by `h_f`. These are alternative exact initial families,
conditional on the same theorem, and are not the selected weighted curve.
Their zero first-scalar difference and nonzero second metric demonstrate
why first-pulse/reflection/zero-TT data alone cannot select a scalar history.

## 7. Checks and scope

`checks/check_weighted_mean.py` verifies radial linearization, adjoints,
canonical conversion, the source reductions, moment and kernel identities
symbolically. It includes deliberate wrong-sign/factor controls and runs
under normal and optimized Python. The checks do not solve the singular
weighted boundary-value problem, establish the retained analytic inverse,
or provide a shell response or hot-matter result. The source model remains
the same pure-tension Einstein-scalar model. This is internal AI-assisted
research and review, not external specialist peer review.

Source key: S1 = exact-family TT `DERIVATION.md` (20261031); S2 = exact
localized initial-data theorem (20261014); S3 = explicit physical initial-jet
derivation (20261015), used only for the Hessian; S4 = localized constraint
manuscript (20261014); S5 = retained shell proper-time observable; S6 =
actual potential/background coefficient certificate retained by S3. Exact
paths and hashes are in `SOURCE_PINS.json`.

Public source links at the pinned commit:

- [S1: exact-family TT transfer](https://github.com/maldonado-research/HDblast/blob/b51908d7e62e7ec80cef82438edd00923ff94723/research/HDBLAST_EXACT_FAMILY_TT_TRANSFER_20261031/DERIVATION.md).
- [S2: exact localized initial theorem](https://github.com/maldonado-research/HDblast/blob/b51908d7e62e7ec80cef82438edd00923ff94723/research/HDBLAST_EXACT_LOCALIZED_PRECONTACT_20261014/EXACT_LOCALIZED_INITIAL_DATA_AND_PRECONTACT.md).
- [S3: physical initial-jet Hessian](https://github.com/maldonado-research/HDblast/blob/b51908d7e62e7ec80cef82438edd00923ff94723/research/HDBLAST_EXPLICIT_SECOND_INITIAL_JET_20261015/EXPLICIT_COMPACT_SECOND_INITIAL_JET.md).
- [S4: localized constraint and moment manuscript](https://github.com/maldonado-research/HDblast/blob/b51908d7e62e7ec80cef82438edd00923ff94723/research/HDBLAST_SECOND_ORDER_LOCALIZED_20261014/manuscript/SECOND_ORDER_LOCALIZED.md).
- [S5: shell proper-time observable](https://github.com/maldonado-research/HDblast/blob/b51908d7e62e7ec80cef82438edd00923ff94723/research/HDBLAST_TOTAL_RESPONSE_PREFLIGHT_20261023/observable/SHELL_PROPER_TIME_OBSERVABLE.md).
- [S6: retained background coefficient certificate](https://github.com/maldonado-research/HDblast/blob/b51908d7e62e7ec80cef82438edd00923ff94723/research/HDBLAST_EXPLICIT_SECOND_INITIAL_JET_20261015/inputs/ACTUAL_POTENTIAL_JET_BOUNDS.json).
