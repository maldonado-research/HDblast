# First-master interior derivatives through order three

This is a conditional exact rational certificate for the retained smooth actual
first master. It extends the separately retained first-jet work. It does not
evaluate a second-order solution or the full quadratic bulk source. The original
coordinate-source checkpoint is unchanged.

Put `ell=10^-6`, `x=t/ell`, `r=z/ell`, and `U(x,r)=w(ell*x,ell*r)`. Every bound
below holds on the backward-closed triangle

\[
 D=\{0\le x\le3,\ r\ge0,\ x+r\le8\}.
\]

The exact equation and initial data are

\[
 U_{xx}-U_{rr}+f=0,\quad f=\ell^2\mathcal W(\ell x,\ell r)U,
 \quad U(0,r)=Q(r),\quad U_x(0,r)=\ell V(r),
 \quad V=e^{-\ell r}Q,
\]

where `Q=exp(4-1/((r-1)(2-r)))` on `(1,2)` and is zero elsewhere.
The boundary value is the retained actual `y(x)=U(x,0)`. Its trace is
`H=U_r(x,0)-ell*gamma*y`. The retained existence/compatibility argument is
an analytic hypothesis of this certificate; the interval cells enclose that
actual function. Independent choices inside different cells need not define a
continuous solution and cannot be differentiated as though they did.

## Exact result

Each row bounds **every mixed partial** of the indicated total order. Values
are strict rational upper bounds chosen for display; the JSON receipts retain
the sharper exact fractions used by later source calculations.

| Scaled total order | Uniform bound on D |
| --- | ---: |
| 0 | 1.000193 |
| 1 | 16.000488 |
| 2 | 256.002081 |
| 3 | 6665.267480 |

For physical `t,z` derivatives of total order `k`, divide the scaled bound by
`ell^k`. Additional boundary bounds are

\[
 \|y\|\le0.000188884282937508,\qquad
 \|y'\|\le0.000486339292639931,
\]
\[
 \|y''\|<0.002059031,\quad
 \|H'\|<256.002081,\quad
 \|y'''\|<0.124286658.
\]

These are global enclosures, not pointwise derivative tables or sign claims.

## Reflected higher derivatives and the corner

For physical forcing `f(x,r)` on `r>=0`, let `fo` be its odd extension and

\[
 Kf(x,r)=\frac12\int_0^x\int_{r-x+s}^{r+x-s}f_o(s,q)\,dq\,ds,
 \qquad U=U_0-Kf.
\]

The forcing need not vanish at the shell at positive time. Consequently its
odd extension must not be differentiated as a smooth function:
`d_r fo = f_r(s,abs(q)) + 2 f(s,0) delta_0` in distributions. Set

\[
 a_\pm=r\pm(x-s),\quad g(s,q)=f_r(s,|q|),\quad
 h(s,q)=\operatorname{sgn}(q)f_{rr}(s,|q|),
\]
\[
 B=\mathbf1_{x>r}f(x-r,0),\qquad
 B_1=\mathbf1_{x>r}f_x(x-r,0).
\]

For `r>0` away from `x=r`, the exact identities are

\[
 K_{rr}=\tfrac12\int_0^x[g(s,a_+)-g(s,a_-)]ds-B,
 \quad K_{xr}=\tfrac12\int_0^x[g(s,a_+)+g(s,a_-)]ds+B,
 \quad K_{xx}=f+K_{rr},
\]
\[
 K_{rrr}=\tfrac12\int_0^x[h(s,a_+)-h(s,a_-)]ds+B_1,
 \quad K_{xrr}=\tfrac12\int_0^x[h(s,a_+)+h(s,a_-)]ds-B_1,
\]
\[
 K_{xxr}=f_r+K_{rrr},\qquad K_{xxx}=f_x+K_{xrr}.
\]

One obtains the first two formulas by differentiating the endpoint-ray
integrals and retaining the crossing of the odd forcing's shell jump. The
function `g` is continuous across zero, so its next derivative is `h` without
another delta term. Differentiating `B` gives the displayed `B1` signs;
it would additionally give a characteristic delta if `f(0,0)` were nonzero.

The retained initial collar is flat near `r=0`: `Q` and `V` vanish there.
Finite propagation and the zero dynamic-boundary data imply `y=0` for `x<1`.
Thus `y^(j)(0)=0`, and `f(0,0)=f_x(0,0)=0`. These are sufficient for `K`
and the causal lift separately to be C3 across `x=r`. Explicitly, the inner
minus outer second-derivative jumps are `(-f00,+f00,-f00)` in the order
`rr,xr,xx`; with `f00=0`, the third jumps are
`(+fx00,-fx00,+fx00,-fx00)` in the order `rrr,xrr,xxr,xxx`.
They all vanish here. General nonflat compatible data can instead use
cancellation between the lift and `K`; no necessity of separate flatness
is claimed for those different problems.

At `r=0`, derivatives are one-sided and the joins give

\[
 K_{rr}=-f,\quad K_{xx}=0,\quad K_{xrr}=-f_x,\quad K_{xxx}=0,
\]
\[
 K_{rrr}=\int_0^x f_{rr}(s,x-s)ds+f_x(x,0),
 \qquad K_{xxr}=f_r+K_{rrr}.
\]

Every reflected point satisfies `abs(q)+s<=x+r<=8`. No physical potential at
negative depth or beyond the audited collar enters these identities.
Smoothness of the actual physical forcing makes the ray formulas continuous
within each branch; the explicit zero jumps and one-sided traces extend the
displayed derivatives continuously to the whole triangle.

## The incoming route to Hprime

Set `J=U_x+U_r`. The equation gives `(d_x-d_r)J=-f`. On its incoming ray,

\[
 H=J_0-y'-\ell\gamma y-\int_0^x f(s,x-s)ds,
 \quad J_0=Q'+\ell e^{-\ell x}Q.
\]

Leibniz differentiation of this actual smooth integral gives

\[
 H'=Q''+\ell e^{-\ell x}(Q'-\ell Q)-y''-\ell\gamma y'
       -f(x,0)-\int_0^x f_r(s,x-s)ds.
\]

The endpoint `f(x,0)` is essential. The dynamic boundary equation supplies

\[
 y''=(\ell/\mu)H-Ay,\quad
 y'''=(\ell/\mu)H'-A'y-Ay',
\]
\[
 A=\ell^2[5/2-(35/4)\tanh^2(\ell x)],\quad
 |A|\le(5/2)\ell^2,\quad
 |A'|\le(35/2)\ell^3(3\ell).
\]

The checker reads all 4,353 exact retained H cells. It bounds the first and
second primitives of `H` by integrating their whole-cell endpoints, including
each quadratic primitive's possible interior vertex. The positive interval
`ell/mu` then gives `F=(ell/mu)int_0^x(x-s)H(s)ds` and `F'`. With
`q_b=(5/2)ell^2*9/2`, it computes

\[
 B_y=\|F\|/(1-q_b),\quad
 B_{y'}=\|F'\|+3(5/2)\ell^2 B_y,\quad
 B_{y''}=(\ell/\mu)_{\max}B_H+(5/2)\ell^2 B_y.
\]

No derivative of an H interval selection is taken.

## Bump envelopes and closure of the rational estimates

Write `s=1/((r-1)(2-r))>=4` and `v=3-2r`, with `v^2=1-4/s` and `|v|<=1`.
Direct differentiation gives

\[
 Q'=vs^2e^{4-s},\quad
 Q''=(s^4-6s^3+6s^2)e^{4-s},
\]
\[
 Q'''=v(s^6-10s^5+24s^4-12s^3)e^{4-s}.
\]

For `s>=4`, `|s^2-6s+6|<=s^2` and
`|s^3-10s^2+24s-12|<=s^3`. The checker verifies these inequalities by shifting
`s=u+4` and checking nonnegative coefficients of the relevant sums and
differences. The maxima of `s^k e^(4-s)` occur at `max(4,k)`. Since
`exp(2)>sum_(j=0)^4 2^j/j!=7`, this proves

\[
 B_Q=1,\quad B_{Q'}=16,\quad B_{Q''}=256,\quad B_{Q'''}=46656/7.
\]

Flatness of the bump at its endpoints validates its smooth zero extension.
Also `B_V<=1`, `B_V'<=16+ell` and
`B_V''<=256+32ell+ell^2`. The free field

\[
 U_0=\tfrac12[Q_o(r+x)+Q_o(r-x)]
 +\tfrac\ell2\int_{r-x}^{r+x}V_o(q)dq
 +\mathbf1_{x\ge r}y(x-r)
\]

therefore has every order-k mixed derivative bounded by
`B_Q^(k)+ell*B_V^(k-1)+B_y^(k)` for `k=1,2,3`. A radial derivative of the
lift has its factor `-1`; taking absolute values does not remove the need for
that sign in the reflected identities or corner joins.

The retained actual background and new independent potential certificate give

\[
 M=102000,\quad M_z=30000000,\quad M_t=21/400000,
 \quad M_{zz}=30000000000.
\]

These are physical derivatives of this actual potential, derived from its
registered background ODE and polynomial potential. The higher-jet checker
authenticates that report and its matching actual-background input; it does
not claim to rederive the Wzz interval. The capsule's potential checker does.

With `q=M(3ell)^2/2=459/10^9`, the C0 and first-jet formulas give

\[
 B_U=(1+3\ell+B_y)/(1-q),\quad
 B_1=16+\ell+B_{y'}+3\ell^2MB_U.
\]

Define `F=ell^2 M B_U`, and use subscripts `0` for shell bounds:

\[
 F_r=\ell^3M_zB_U+\ell^2MB_1,\quad
 F_x=\ell^3M_tB_U+\ell^2MB_1,
\]
\[
 F_0=\ell^2MB_y,\quad
 F_{x0}=\ell^3M_tB_y+\ell^2MB_{y'}.
\]

The reflected second errors for `U-U0` are
`E_rr=E_xr=3Fr+F0`, `E_xx=F+3Fr+F0`. Add the common free second bound to
obtain `B2_rr,B2_xr,B2_xx`. The incoming formula bounds `Hprime`, then
the boundary ODE bounds `ythird`, using those already certified first jets.
There is no circular dependence. Finally

\[
 F_{rr}=\ell^4M_{zz}B_U+2\ell^3M_zB_1+\ell^2MB_{2,rr}.
\]

The third errors are `E_rrr=E_xrr=3Frr+Fx0`,
`E_xxr=Fr+3Frr+Fx0`, `E_xxx=Fx+3Frr+Fx0`. Adding the common free third
bound yields the table. This route needs no Wtt or Wtz bound.

## Reproduction and finite-check limits

From the standalone capsule root:

```sh
python -I -B jets/check_higher_jets.py --root . --output /tmp/HIGHER_JETS_NORMAL.json
python -I -O -B jets/check_higher_jets.py --root . --output /tmp/HIGHER_JETS_OPTIMIZED.json
python -I -B jets/check_higher_jets.py --run-controls --output /tmp/HIGHER_JETS_CONTROLS.json
```

The program imports only the Python standard library. It resolves all retained
inputs relative to `--root`, pins five minimal retained files plus the new
potential report, and writes only the requested output. Its independent
polynomial checks use exact tensor Lagrange differentiation of the integrated
kernel, separately from the reflected derivative formulas. Twelve physical
monomials cover reflected/unreflected points, shell traces, both characteristic
traces and compatible joins. They support the analytic argument above; finite
witnesses do not replace the all-point proof. Eleven intentional corruptions
must each fail by an explicit exception in both normal and optimized Python.

All first-master derivatives needed before a direct longitudinal bulk Hessian
evaluation are now bounded. Metric/scalar reconstruction coefficients, a
consistent first gauge bridge, full quadratic-source evaluation, compatible
second initial data, coupled second evolution, and a nonlinear remainder are
separate tasks. This certificate alone establishes none of those later steps.
