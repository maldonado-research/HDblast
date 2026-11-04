# Actual retained first-packet forcing: a continuous future-cone ledger

4 November 2026. Portable conditional source derivation and exact postprocessing. No
retained source was edited and no background, bump, master, gauge, source,
or second-evolution oracle was called. This note uses the existing actual
solution's continuous enclosures; it does not establish a past-retarded
preparation law.

## Registered field identity

The native update in `HDBLAST_EXACT_CAUCHY_PACKET_20261031/README.md` changes
the initial velocity of a **second-response approximation**. Its actual
second target, first packet, forcing and operators remain fixed. In
particular it does not replace the actual first velocity by zero.

The current exact-family seed is explicitly given in
`HDBLAST_EXACT_FAMILY_TT_TRANSFER_20261031/DERIVATION.md`, section 2:

\[
u=\alpha Q/\rho^2,\quad u_t=E u,\qquad
v=B(Q_z+kQ),\quad v_t=BE[Q_z+(k-1)Q],
\quad B=3\alpha/(\rho^3p),\quad E=e^{-z}.
\]

The prepared initial executable `initial/scientific.py:first` uses exactly
`u=C*q`, `v=B*(q.dz()+k*q)`, `ut=E*u`, and
`vt=E*(v-3*u/(rho*p))`. The identity `B Q=3u/(rho p)` proves the displayed
velocity equality. Its replay wrapper imports the geometric source
checkpoint's same `numeric/interval_core.py`, with pinned actual shell,
background ODE, bump and derivative data. No nominal midpoint or completed
second-candidate exterior is substituted for an actual first field.

The old master certificate uses `W(0,r)=Q(r)` and
`W_x(0,r)=ell*exp(-ell*r)*Q(r)`. Its correlated reconstruction is
`u=(alpha/rho^2)cosh(t)^(-3/2)W` and
`v=B cosh(t)^(-3/2)(W_r/ell+kW)`. Since the time derivative of the cosh
factor vanishes initially, these are the same four physical first profiles.

Values alone are insufficient to bind the sources. The actual common
harmonic first gauge has `f=b=d=f_t=b_t=d_t=0` initially, but
`f_tt=-6rho^2(u_y+4au)`, giving the nonzero mixed metric velocity
`q_ty,t=-6rho^2(u_y+4au)Y`. The Oct31 exact-family derivation, section 4,
matches this full one-jet and invokes the retained mixed linear uniqueness
to identify the full evolved master/gauge pair. The Nov1 prepared gauge
derivation, section 7, applies the same identification with a changed
second initial curve. This is an inherited conditional analytic identity;
the present arithmetic does not reprove the PDE or uniqueness theorem.

The prepared published `gauge/SOURCE_PINS.json` binds that Oct31 derivation
and the old harmonic completion and full response proof. Separately, the
old numeric predecessor pins bind exactly the same common first-fields
proof bytes. `source_budget.py` checks these selected joins
before and after its arithmetic. The accompanying `REGISTERED_PACKET_AUDIT.md` records the wider
code/input identity trace. A genuinely changed
first packet, including a new choice `W_x(0)=0`, would fail this identity
and needs another source certificate.

## Source adapter, signs, normalizations and regularity

Put `ell=10^-6`, `x=t/ell`, `r=z/ell`, `dy=-rho dz`.
Thus `partial_x=ell partial_t` and `partial_r=-ell rho partial_y`.
Parameter expansions use `X_e=X0+e X1+e^2 X2/2+o(e^2)` throughout.
The quadratic sources are second parameter Hessians; their coefficients
in an ordinary `e^2` polynomial are half these values.

In the parallel background frame the five means are
`U=(K0,L0,N0,J0,w0)`, with positive probability-sphere Gram
`diag(1,2,1,3,1)`. The closed subsystem is

\[
U_{xx}-U_{rr}=A_0U_x+B_0U_r+D_0U+S_0,\qquad
L_0(x,0)=0,\qquad P_NU_r=M_0P_NU+G_0.
\]

The actual bulk source, from the independent-jet geometric Laurent map,
is exactly

\[
S_0=P_0(-2\ell^2\rho^2Q^{\rm red},\;\ell^2\rho^2Q_\phi).
\]

Explicitly its five entries are the constant scalar/tensor coefficients
of `(-2ell^2rho^2 QEred_yy, -2ell^2rho^2 QEred_y0,
-2ell^2rho^2 QEred_00, -2ell^2rho^2 QEred_ab,
+ell^2rho^2 Qphi)`. `QEred` includes the geometric Ricci and gauge terms,
not just the non-Ricci matter adapter. The exact source map is the retained
`projection/GEOMETRIC_COEFFICIENT_MAPS.json`, whose SHA is recorded in
the accepted actual report. Its polynomial/Laurent inputs are the seven
first frame amplitudes and their ordinary physical two-jets, the first
scalar two-jet, and actual background jets. These are evaluated as
derivatives of the same actual master and gauge solutions. Normalized
Taylor coefficients are converted to ordinary derivatives with their
factorials before source substitution. A derivative of an independently
selected interval error is never used.

For the physical shell sources use the one-copy outward affine Hessians
`QR`, `QP`, and the actual-metric-lowered quadratic gauge covector `QC`.
The outward shell coordinate rows before conversion to inward r are

\[
K_y:\ GN,\qquad N_y:\ -2QR_{00},\qquad
J_y:\ -2QR_J,\qquad w_y:\ -QP,
\quad GN=-2\operatorname{tr}_\gamma QR-2(QC)_y.
\]

Consequently the actual boundary forcing in mean order is

\[
G_0=\ell\rho_b(-GN_0,\;0,\;2QR_{00,0},\;2QR_{J,0},\;QP_0),
\quad GN_0=2QR_{00,0}-6QR_{J,0}-2QC_{y,0}.
\]

This is the source of the **boundary derivative**, not its time integral.
The reflected lift is `Lg=-integral_0^max(x-r,0)g(s)ds`, whose r derivative
equals g at the shell; it introduces no additional source normalization.
The L Dirichlet condition is homogeneous. No derivative of G is needed
by the accompanying forced-energy estimate.

There is also a concrete closed shell functional of the actual first
traces. In the numeric amplitude convention
`H=A`, `m00=D Y`, `mab=L Y*kappa+M Hess(Y)`, `eta=sY`, set

\[
c=aH^2/4-\sigma'Hs/6-\sigma''s^2/6,\quad
d=aH+\sigma's/3,\quad
q=pH^2/4+\sigma''Hs/2+2s^2.
\]

Then `QR00,0=(-c-dD)/6`,
`QRJ,0=(c-dL)/6+4dM/9`, and `QP0=q/6`.
These use the actual on-shell first physical rows only; none is
differentiated normally. `QC_y,0` is supplied by the same geometric
source map, so G is a fully specified quadratic functional of the actual
first field, not a newly named unspecified norm.

The existing actual master/gauge C3 derivative certificates give all
ordinary first metric/scalar two-jets used here. They retain smooth
background/master/gauge existence and first-equation premises. The later
exact-family C4 spacetime control suffices for the differentiated equations,
sphere projection, traces and the continuous energy ledger. None of these
analytic regularity premises is established by the present finite joins.

## Continuous history already available on the future cone

The retained actual source report covers
`0<=x<=3, r>=0, x+r<=8` by 81 rectangles of side 1/2 clipped at the last
diagonal. Each reported interval encloses the actual source throughout its
whole clipped cell, including shell and reflected cells. It supplies all
bulk angular coefficients and six whole-time physical shell cells on
`0<=x<=3`. Its actual field errors and background uncertainty remain.

The global positive bundle bounds are
`F_H<=17387513156573` and `G_H<=3851400`. The bulk number is the sum of
probability-sphere angular Gram block bounds, with mixed metric blocks
counted twice. The shell join includes physical metric/scalar forcing
and `GN`, using the actual harmonic gauge source. Orthogonal invariant
projection is contractive in this same positive norm.

Restrict to `D_T={x,r>=0,x+r<=T}`, `T=5/4`. Both the retained actual
quadratic sources and all first jets vanish exactly for `x+r<=1` by the
flat initial collar and finite propagation. Hence the uniform safe bounds
are

\[
\int_{D_T}|S_0|_H^2\,dx\,dr\le\frac9{32}F_H^2,
\qquad \int_0^T|G_0|_H^2\,dx\le\frac14G_H^2.
\]

Without using the exact zero wedge the equally safe geometric factors
are `25/32` and `5/4`. The postprocessor records both.

The new arithmetic clips each existing rectangle with the two exact
half-planes `x+r<=5/4` and `x+r>=1`, using rational polygon intersections
and the shoelace area formula. The old grid has pairwise disjoint
interiors, so shared edges have measure zero and are not double counted.
The five positive-area active pieces have total area exactly `9/32`.
For each piece it independently forms the direct mean Gram square from
the stored source intervals; it also forms the full angular Gram square
using the retained exact Gram certificates. The smaller of each valid
cap and the global cap is used. There is no new actual source call.

On the shell the old physical QR/QP intervals are combined with QC_y from
the matching old r=0 boundary rectangle through the GN identity above.
This is a continuous conservative shell enclosure. The global G bound
wins on the sole active time interval `[1,5/4]`; no sharper G certificate
is asserted. In particular shell QR alone is not substituted for all G.

The resulting exact rational bounds are in `references/NORMAL.json`:

| Continuous source integral | Global active-cone cap | Existing-cell cap |
| --- | ---: | ---: |
| Five-mean bulk squared L2 | about 8.5029078873e25 | about 8.5264181293e24 |
| Full-angular bulk squared L2 | about 8.5029078873e25 | about 2.3021328949e25 |
| Five-mean boundary squared L2 | 3708320490000 | 3708320490000 |

The decimals are explanatory approximations; the JSON fractions are the
certificates. This is an integral upper bound for the actual continuous
source, not numerical integration of a retarded response or an exact
value of the integral.

## Forced energy consequence and limits

Conditional on the accompanying independently reviewed forced-energy
proof, let

\[
\mathcal B_S=\frac12\int_{D_T}|S_0|_H^2+
\frac{5002001}{4000000}\int_0^T|G_0|_H^2.
\]

That proof yields `E'<=3.01E+0.5||S0||^2+Cg|G0|^2`,
`exp(3.01T)<44`, so

\[
E(x)\le44(E(0)+\mathcal B_S),\qquad
|P_NU(x,0)|_H^2\le176(E(0)+\mathcal B_S).
\]

The cellwise ledger has `B_S` about `4.2632090646e24`. For the auxiliary
zero-initial reduced response only, this gives a trace upper estimate
about `2.7392057159e13`, compared with about `8.6501785767e13` from the
uniform energy cap. It remains looser than the retained sup-norm
Duhamel bound of about `5.4e12`. No improved response or observable bound
is claimed. The useful new result is a justified actual integrated source
ledger and a route compatible with the weaker energy framework.

The actual absolute second initial datum cannot be set to zero: its
constraint is `P x2=-D^2 C[x1,x1]`, with a retained nonzero source witness.
The present zero-initial columns are explicitly auxiliary reduced-equation
calculations, not an exact-family preparation. The exact open zero shell
collar supplies G(0)=G_x(0)=0 and the retained corner identities, but does
not remove the bulk initial constraint away from that collar.

A physically selected past-retarded preparation still needs a specified
compatible past initial or incoming state and its constraint/gauge
completion, the corresponding actual first solution on that past domain,
background/operator coverage there, continuous quadratic bulk and shell
source boxes there, and a causal identification of the resulting second
state at x=0. The present source certificate is future-only (`x>=0`);
time symmetry of the background does not certify a changed or time-reversed
first packet. A finite past slab estimate must also retain its genuine
initial energy rather than silently let it vanish.

Finally, an absolute proper-volume or expansion coefficient includes the
retained first-only clock, measure, and Gauss terms as well as U. This
source-energy ledger alone proves no physical sign, Noether law,
radiation flux, matter production, SI scale, common nonlinear parameter
window, or origin model.

## Replay

From any complete checkout with the pinned old sources, use a fresh
absolute output outside both that checkout and this component:

```
python -I -S -B source_budget.py --repo-root REPO --output OUT/NORMAL.json
python -I -S -B -O source_budget.py --repo-root REPO --output OUT/OPTIMIZED.json
```

Both reports have identical bytes. These two processes perform exact
postprocessing only: zero project imports, source/background/bump queries,
or PDE solves. The strict wrapper retains the two canonical children and separately
records the eight actual rejection subprocesses.
