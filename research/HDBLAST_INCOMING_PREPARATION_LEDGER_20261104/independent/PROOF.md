# Independent source ledger and preparation checks

4 November 2026. Conditional mathematical reconstruction of retained inputs.
The companion executable uses Python's standard library and exact rational
arithmetic, without importing the project's scientific implementations. Its
45 named gates authenticate the selected source bytes, reconstruct four source
budgets, and compare the shared consequences of the primary reports. These
finite checks supplement the analytic arguments; they do not reprove the
retained continuum PDE, exact-family or coefficient-enclosure premises.

## Forced energy with the actual shell rows

The retained five-mean problem has H=diag(1,2,1,3,1),

```
U_xx-U_rr=A U_x+B U_r+D U+S,
L(x,0)=0, PN U_r=M PN U+g.
```

The L Dirichlet slot is homogeneous. The quadratic scalar/metric/gauge terms
remain in S and g; they cannot be dropped when the first perturbation is
nonzero. In increasing-r coordinates,

```
S=P0(-2 ell^2 rho^2 QEred, ell^2 rho^2 Qphi),
g=(-ell rho GN,0,2ell rho QR00,2ell rho QRJ,ell rho QP),
GN=2QR00-6QRJ-2(QC)_y.
```

Here (QC)_y is the covector component, not a new y derivative. These are second
parameter derivatives, with no extra division by two. The old Oct31 native
velocity correction changed a second-response approximation, while the first
state and first quadratic sources stayed fixed; the selected Oct31 exact-family
and Nov1 prepared-family derivations supply that retained binding.

Put V=U_x,Q=U_r,k=alpha=1/2 and u=PN U(x,0). On the shrinking future cone
0<=x<=T=5/4, 0<=r<=T-x define

```
e=(|V|_H^2+|Q|_H^2+2k<V,Q>_H+|U|_H^2)/2,
E=integral e dr+alpha|u|_H^2/2.
```

The retained conditional coefficient caps are
||A||,||B||<=1/10000, ||D||<=1/1000000, ||M||<=1/2000.
They apply on the declared future domain. No past extension is inferred.
The exact bulk identity

```
2e-|V+kQ|_H^2=(1-k^2)|Q|_H^2+|U|_H^2
```

gives <V+kQ,S>_H<=e+|S|_H^2/2. At the shell set P=PN(V+Q), so
u_x=P-Mu-g. Expanding the full nonsymmetric Robin expression and subtracting
its homogeneous counterpart gives exactly

```
-(1/2)<P,g>+<Mu,g>+(1/2)|g|^2-(1/2)<u,g>.
```

The homogeneous boundary proof leaves -|P|^2/8. With m=1/2000,

```
-|P|^2/8+(1/2)|P||g| <= |g|^2/2,
(m+1/2)|u||g| <= |u|^2/4+(m+1/2)^2|g|^2.
```

The bulk and boundary homogeneous rates are respectively 1000351/500000 and
2003001/1000000. Adding one energy rate to each remains below 301/100. Thus
every existing compatible regular solution satisfies

```
E' <= (301/100) E +(1/2)||S||_L2_r,H^2
                     +(5002001/4000000)|g|_H^2.
```

The outer moving characteristic flux and homogeneous Dirichlet flux are
nonpositive. No outer boundary condition or symmetry of M is required. An
independent order-24 rational Taylor upper bound gives exp(301/80)<44. Hence

```
B=E(0)+(1/2) integral |S|_H^2 +(5002001/4000000) integral |g|_H^2,
E(x)<=44B, |w(x,0)|^2<=176B.
```

This is an estimating norm. It is not the physical Noether energy or a measured
stress flux. An auxiliary zero-initial reduced response is not by itself a
constraint-compatible physical second state.

## Independent source integration

For a retained rectangular cell [a,b]x[c,d], define z_+=max(z,0). Its area below
x+r<=s is exactly

```
((s-a-c)_+^2-(s-b-c)_+^2-(s-a-d)_+^2+(s-b-d)_+^2)/2.
```

This translated-triangle inclusion/exclusion route is independent of the primary
implementation's polygon clipping. Subtracting the corresponding s=1 area
from the s=T area implements the retained exact support wedge x+r<=1. The
whole cone has area25/32, the active part9/32. Five active half-unit cells have
areas3/32,3/32,1/32,1/32,1/32; one additional cone cell has zero active area.

The direct mean square cap is formed from the interval means of
(metric yy,y0,00,ab; scalar), with weights(1,2,1,3,1), then bounded also by the
retained full-source global cap FH^2. Integrating each cap over its exact cell
area gives the same four full/active, global/cellwise budgets as the primary
report. The full angular norm is positive and its invariant mean is orthogonal;
the mean projection is contractive. The independent direct-mean route gives
the same chosen caps for these cells.

Boundary caps retain the actual GN combination above, including its QC term.
The normalized source square is bounded by

```
(ell*131)^2 [GN^2+4QR00^2+12QRJ^2+QP^2],
```

using interval upper magnitudes, and also by GH^2. The active boundary length
is1/4, and the global GH cap is the smaller active bound. All endpoints and
source budgets are exact fractions reconstructed from old continuous boxes;
there are no new background, bump, evolved-field or source-oracle queries.
The resulting auxiliary trace upper bound remains looser than the retained
Duhamel bound, so this ledger does not improve the numerical prediction.

## Current packet diagnostics and the stronger two-endpoint inequality

The registered first packet has W(0,r)=Q(r), W_x(0,r)=ell exp(-ell r)Q(r),
ell=10^-6, with the flat bump supported on(1,2). At r=3/2, Q=1,Q_r=0;
both principal components J+=W_x+W_r and J-=W_x-W_r are strictly positive.
For the integrals over the support, endpoint flatness gives

```
integral(J+^2-J-^2)=2ell^2 integral exp(-ell r)Q^2 >0,
integral(J+^2+J-^2)=2 integral(ell^2 exp(-2ell r)Q^2+Q_r^2).
```

On a unit interval, write
Q(x)=(1-x) integral_0^x Q' -x integral_x^1 Q'. Weighted Cauchy implies
Q(x)^2<=x(1-x) integral_0^1 Q'^2. Integrating x(1-x) gives1/6, so the
relative principal imbalance is less than ell^2/6. This implies the common
publication guard ell^2/4=1/4000000000000. The stronger route need not share
the primary proof's intermediate constants.

An entire-current-slice one-way condition would impose Q'=+/-ell exp(-ell r)Q
with Q(1)=0, forcing Q=0 by linear ODE uniqueness. The registered packet is
therefore mixed at the current slice. This does not exclude an incoming-only
past state: potential scattering and the physical shell can produce both
current components. Principal combinations are not a universal exact radiation
condition for a wave equation with a potential.

## Affine null and normalized scalar-history scope

For twist-free affine null cuts with three-dimensional screens, the second
Raychaudhuri variation about shear zero is

```
theta2'+(2theta0/3)theta2
 =-2theta1^2/3-2|sigma1|^2-2(Phi1')^2-2Phi0'Phi2'.
```

Phi denotes the scalar pulled back to the physical affine cuts, including the
embedding. Under the newly declared cut condition Phi2=0 and zero corner
theta2, the integrating factor is R0^2 and theta2 is a nonpositive weighted
integral. The first constraint must also hold. For log density L, the actual
area second derivative contains average(L2+L1^2), so its sign does not follow
from L2<=0. These are partial null-constraint identities, not a completed
characteristic solution or a proper-shell expansion/H1 sign bound.

For the retained admissible h-only family, let w=-p zeta initially. Substituting
D=3h/2-3aw/p gives z_r=3bz-h and h_r=3ell^2z-bh+cw,
c=2ell rho p/3, with z(1)=h(1)=0. Retained initial bounds give

```
cM=131/10500000000, Hw=cM/3,
Hrw=(3ell^2+bM)Hw+cM, Vw=3ell^2Hw+cM,
N<=max((1+2Hw+2T Vw)/3,1+2Hrw+2Vw)Y(w)<=Nguard Y(w),
Nguard=20000001/20000000, q=1201/3200000,
kappa=q Nguard/(1-q)=24020001201/63975980000000<47/125000<1.
```

The numerical majorant multiplying Y(w) is strictly below Nguard; the last
inequality is strict when Y(w)>0 and includes equality for the zero profile.
The reflected free scalar trace is w(x); the full coupled retained resolvent
and exact zero wedge imply Y(m-w)<=kappa Y(w). Therefore
(1-kappa)Y(w)<=Y(m)<=(1+kappa)Y(w). This proves causal injectivity only on
the retained admissible smooth image, not surjectivity onto arbitrary histories.

The observable m is the difference of second derivatives of the normalized
proper-clock shell scalar mean (integral phi dV)/V. Formally differentiating
this quotient cancels phi_b times the second volume variation. Common first
measure/scalar/clock products cancel between the two families, background
phi_b is constant along the shell, and the retained pure second normal gauge
component is zero there. An unnormalized integral does not have this property.
No past state, independent physical preparation law, unrestricted inverse or
cosmological origin follows from the estimate.

## Selected retained inputs

- [Five-mean physical and gauge rows](../../HDBLAST_PREPARED_SHELL_VOLUME_20261101/gauge/MEAN_INITIAL_GAUGE.md).
- [Conditional future operator premises](../../HDBLAST_TOTAL_RESPONSE_PREFLIGHT_20261023/analytic/EXPLICIT_PREMISES.json).
- [Old continuous actual-source boxes](../../HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/numeric/ACTUAL_SOURCES_NORMAL.json).
- [First packet and higher jets](../../HDBLAST_COMMON_GAUGE_FORCING_20261016/inputs/HIGHER_JETS.md).
- [Oct31 exact-family first-jet identification](../../HDBLAST_EXACT_FAMILY_TT_TRANSFER_20261031/DERIVATION.md).
- [Retained H1 energy argument](../../HDBLAST_H1_PREPARATION_RESPONSE_20261104/energy/PROOF.md).
- [Retained initial constants](../../HDBLAST_INTRINSIC_PREPARATION_ROBUSTNESS_20261104/bounds/references/NORMAL.json).
