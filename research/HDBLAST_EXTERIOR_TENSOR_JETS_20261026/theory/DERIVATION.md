# Exterior TT coefficient two-jet: differentiated scalar resolvent

This note derives an implementable interval ledger for the actual total,
registered-gauge TT coefficient. It retains the smooth formal background,
master, gauge, source identities, and second-response existence premises of
`HDBLAST_FULL_ANGULAR_RESPONSE_20261026`. It neither replaces the actual
first solution by its comparator nor asserts gauge-invariant curvature.
The two TT equations, sources, and zero initial data coincide, so the same
coefficient and all its derivatives occur in both TT slots by uniqueness.

## 1. Domain and equation

Set `ell=1/1000000`, `T=1/4`, `r0=5/4`, and

```
C = {(x,r): 0 <= x <= T, |r-r0| <= T-x}.
u_xx-u_rr = A u_x + B u_r + D u + F,
u(0,r)=u_x(0,r)=0.
```

Every past characteristic cone of a point of C is contained in C. The
domain has `r-x>=1`; neither the shell nor any reflected contribution is
in its domain of dependence. Thus the exterior solution of this scalar
problem is the actual total TT coefficient, with no free or shell term.
The complete 78-term authenticated `metric_forcing/ab/5` row is F; row 6
is identical. The zero-data premise holds identically on an initial
interval, not only at finitely many points.

The scalar coefficients inherited from the full operator are

```
A = -3 ell H,                         H=tanh(ell*x), c=cosh(ell*x),
B = -3 ell rho a,
D = ell^2[-6/c^2 + rho^2(-2a^2+2h^2-2b-8d+4V/3)],
h = H/rho,  a_y=-rho^-2-p^2/3,  b=-(a_y+a^2),  d=rho^-2-a^2.
```

Here physical radial coordinate y obeys `dy=-ell*rho dr`. All quantities
below are evaluated on C or its explicitly certified enclosing boxes.

## 2. Wave map and boundary formulas

For continuous f let

```
(Kf)(x,r) = 1/2 int_0^x int_[r-(x-s)]^[r+(x-s)] f(s,z) dz ds.
```

Write `f_+(s)=f(s,r+x-s)` and `f_-(s)=f(s,r-x+s)` at a fixed `(x,r)`.
Direct differentiation gives

```
(Kf)_x = 1/2 int_0^x [f_+(s)+f_-(s)] ds,
(Kf)_r = 1/2 int_0^x [f_+(s)-f_-(s)] ds.
```

If the actual spatial derivative f_r is continuous, a second radial or
mixed derivative follows without using f_x:

```
(Kf)_xr = 1/2 int_0^x [(f_r)_+(s)+(f_r)_-(s)] ds,
(Kf)_rr = 1/2 int_0^x [(f_r)_+(s)-(f_r)_-(s)] ds,
(Kf)_xx = f(x,r)+(Kf)_rr.
```

The endpoint contribution in the last identity is exactly f(x,r).
Differentiating in x first would introduce a nonzero initial acceleration
for u_x and unnecessary F_x input. The spatial derivative avoids both.

For the whole-cone norm

```
X(w)=max(2 sup_C|w|/T, sup_C|w_x|, sup_C|w_r|),
```

the formulas imply `X(Kf)<=T sup_C|f|`: the value bound is
`|Kf(x,r)|<=x^2 sup|f|/2`, and both derivative bounds are at most
`x sup|f|`. The same formulas give the sharper apex value factor T/2
relative to the common first-derivative bound.

## 3. Spatially differentiated equation and resolvent

Let `v=u_r`, and introduce

```
P w = A w_x+B w_r+D w,
Q w = A_r w_x+B_r w_r+D_r w.
```

Differentiation of the retained smooth actual PDE, with commuting ordinary
coordinate partial derivatives, gives

```
u = KF + KP u,
v = KF_r + KP v + KQ u,
v(0,r)=v_x(0,r)=0.
```

There is no new initial datum: both zero statements follow by radial
differentiation of the two identically zero initial functions. In
particular, the generally nonzero `u_xx(0,r)=F(0,r)` is irrelevant to the
initial data for v. The construction only uses continuity of F and F_r
and the required coefficient derivatives in addition to the inherited
smoothness premise; no quantitative F_x bound is used.

Let nonnegative exact rational upper caps be

```
a>=sup|A|, b0>=sup|B|, d0>=sup|D|,
ar>=sup|A_r|, br>=sup|B_r|, dr>=sup|D_r|,
f>=sup|F|, g>=sup|F_r|.
sigma=a+b0+(T/2)d0,         q=T*sigma,
tau=ar+br+(T/2)dr,          s=T*tau.
```

The symbol b0 here is the coefficient cap, distinct from background b.
Require `0<=q<1`. Then `X(KP w)<=q X(w)` and
`X(KQ w)<=s X(w)`. Sequential scalar inversion gives

```
M_u = T f/(1-q),
M_v = (T g+s M_u)/(1-q),
X(u)<=M_u,                  X(v)<=M_v.
```

No additional assumption `q+s<1` is necessary. The block operator is
triangular; its nonnegative norm majorant is `[[q,0],[s,q]]`. Its
n-th power for n>=1 is `[[q^n,0],[n*s*q^(n-1),q^n]]`, and its
resolvent is `[[1/(1-q),0],[s/(1-q)^2,1/(1-q)]]`. This proves the
two displayed bounds component by component even if s exceeds one.

These inequalities apply to the actual complete scalar solution. They are
not estimates for an auxiliary first-iterate model. They retain every
number of variable-coefficient feedback interactions and all differentiated
coefficient terms.

## 4. Apex enclosure ledger

Define free-source apex intervals by rigorous quadrature:

```
J0  encloses (KF)(T,r0),
Jx  encloses 1/2 int_0^T [F(s,r0+T-s)+F(s,r0-T+s)] ds,
Jr  encloses 1/2 int_0^T [F(s,r0+T-s)-F(s,r0-T+s)] ds,
Jxr encloses 1/2 int_0^T [F_r(s,r0+T-s)+F_r(s,r0-T+s)] ds,
Jrr encloses 1/2 int_0^T [F_r(s,r0+T-s)-F_r(s,r0-T+s)] ds.
```

The complete omitted-feedback allowances are

```
e1 = q M_u,
e0 = (T/2)e1,
e2 = q M_v+s M_u
   = q*T*g/(1-q)+s*T*f/(1-q)^2.
```

Consequently the actual apex jet obeys

```
u    in J0  +[-e0,e0],
u_x  in Jx  +[-e1,e1],
u_r  in Jr  +[-e1,e1],
u_xr in Jxr +[-e2,e2],
u_rr in Jrr +[-e2,e2],
u_xx in I_rr + I_A I_x + I_B I_r + I_D I_0 + I_Fapex.
```

`I_0,I_x,I_r,I_rr` in the last row are already widened total-response
intervals, and `I_A,I_B,I_D,I_Fapex` are certified apex intervals. The
last row uses the exact actual scalar PDE. It is a sound interval for
u_xx; it is not an independent computation of a pointwise PDE residual.
No extra geometric tail is added after e0/e1/e2: those terms already sum
the complete response feedback.

The global value and first-derivative caps also imply
`|u|<=T M_u/2`, `|u_x|,|u_r|<=M_u`, and
`|u_xr|,|u_rr|<=M_v`; intersecting corresponding apex intervals with
these independent bounds is sound. A further optional consistency
intersection comes for free if the cone cells already enclose F_r:

```
u_r=v in (KF_r)(T,r0)+[-T e2/2,T e2/2].
```

This is the same derivative represented by an area integral instead of
the characteristic boundary integral. Empty intersections must fail the
certificate, not silently select one computation.

### Quadrature with no midpoint or smoothness assumption

For every clipped source cell C_j, use its exact positive area alpha_j
and a rectangular enclosure containing C_j. If `F(C_j)` is enclosed in
`[l_j,h_j]`, its J0 contribution is
`[alpha_j*l_j/2,alpha_j*h_j/2]`. Exact cell intersections cover C with
disjoint interiors and total area T^2. Taking the maximum of all cell
absolute caps supplies f and, separately, g for the whole cone.

For each time interval `[s_i,s_(i+1)]`, enclose the whole corresponding
segment of each characteristic in a valid source-oracle box. Multiply
each boundary enclosure by the exact positive weight
`(s_(i+1)-s_i)/2`. Use interval addition for the plus combination and
interval subtraction for the minus combination; subtraction reverses
the second interval's endpoints. The segments have total length T.
These range-enclosure quadratures need no source Hessian, endpoint
interpolation estimate, or numerical derivative of interval values.

## 5. Explicit coefficient derivatives

The background identities give

```
rho_r = -ell rho^2 a,
a_r   = -ell rho a_y,
p_r   = -ell rho p_y,
V_r   = -ell rho p Vp,
A_r   = 0,
B_r   = 3 ell^2 rho^2(a^2+a_y).
```

For a compact D_r expression set

```
S = 8a^2-(2/3)p^2+(4/3)V.
D/ell^2 = -6/c^2-10+2H^2+rho^2 S,
D_r = -ell^3 rho^3[2a S+16a a_y-(4/3)p p_y+(4/3)p Vp].
```

The rewrite uses only the displayed definitions of b,d,h and the exact
background equation for a_y. It does not assume a numerical Hamiltonian
cancellation. Direct differentiation of the unsimplified D expression is
an independent algebraic check. Interval evaluation of either form is
valid; intersecting two certified forms can sharpen the result.

## 6. Why the actual TT forcing derivative needs only inherited C3

A generic derivative of a quadratic metric source involving physical
two-jets of first-order fields can require four master/gauge derivatives,
because some first-order metric amplitudes already contain a derivative.
That generic count is unnecessarily strong for this exact TT source.

An independent scan of all 78 encoded monomials in the frozen source map
finds that the *only* second-derivative field symbols are `Myy`, `Myt`,
and `Mtt`, in zero-based rows 7,22,50 respectively:

```
row  7: -(4/3) A Myy ell^2 rho^2,
row 22:  (8/3) B Myt ell^2 rho,
row 50: -(4/3) D Mtt ell^2.
```

Here uppercase D in row 50 is a first-order metric amplitude, not the
TT operator potential D. The registered construction has `M=2d_gauge`,
with no derivative of the master or gauge in M itself. All other field
symbols in the source have differential order at most one; after one
radial derivative they need only the already certified physical two-jets.
The scalar occurs only as s and needs only its first radial derivative.
Thus the three displayed M terms need gauge d through total scaled order
three, which is included in the inherited actual differentiated C3 input
certificate. No C4 certificate is needed for this particular source row.

Explicitly, writing `d_gauge=z` to avoid all other uses of d,

```
partial_r Mtt = 2 z_xxr/ell^2,
partial_r Myt = -2 z_xrr/(ell^2 rho) - 2a z_xr/ell,
partial_r Myy = 2 z_rrr/(ell^2 rho^2)
                +6a z_rr/(ell rho)+2(a^2-a_y)z_r.
```

For any ordinary field symbol h,h_y,h_t already supplied by the oracle,
`h_r=-ell*rho*h_y`, `(h_y)_r=-ell*rho*h_yy`, and
`(h_t)_r=-ell*rho*h_yt`. Background radial derivatives follow from their
ODEs, including `(Vp)_r=-ell*rho*p*Vpp` and
`(a_y)_r=-2ell*a/rho+(2/3)ell*rho*p*p_y`. The temporal symbols c,theta
and constant ell have zero radial derivative. Differentiate each complete
encoded monomial with these certified actual derivatives. Negative
exponents require validated nonzero denominator intervals. This is
differentiation of exact identities with actual derivative enclosures;
it is never differentiation of freely selected interval endpoints.

The accompanying independent script audits the symbol union, the three
monomial positions, and the registered M construction without importing
the producer or reevaluating the actual source oracle. Full frozen-tree
authentication remains the evaluator's separate prerequisite.

## 7. Natural-coordinate coefficient two-jet and PDE closure

At the apex the ordinary coefficient partial derivatives satisfy

```
u_t  = u_x/ell,
u_y  = -u_r/(ell rho),
u_tt = u_xx/ell^2,
u_ty = -u_xr/(ell^2 rho),
u_yy = u_rr/(ell^2 rho^2)+a u_r/(ell rho).
```

These must use actual apex background enclosures with strictly positive
rho. The `a u_r/(ell rho)` correction and the mixed-derivative sign
are essential. Substitution into the scaled scalar equation gives

```
rho^-2 u_tt-u_yy
  = -3H rho^-2 u_t +4a u_y
    +[D/(ell^2 rho^2)]u+F/(ell^2 rho^2).
```

This is exact closure of the registered scalar TT coefficient equation
in the natural y,t coordinates. An interval residual assembled from
these dependent enclosures must contain zero, but that containment by
itself is not an accuracy estimate or independent field-equation test.
Because u_xx is reconstructed from the PDE, an algebraic closure check
must be described as such.

These ordinary coefficient partials are not, on their own, the physical
covariant two-jet of the full moving-frame metric, the Einstein tensor,
or a gauge-invariant curvature observable. Such an observable requires
its separately derived projection, frame and connection terms, gauge
transformation law, and all relevant non-TT sectors.

## 8. Optional sharper feedback and declared acceptance

If an actual feedback iterate is evaluated, the same triangular bound
gives a rigorous remainder after n feedback powers:

```
X(u-u_[n]) <= q^(n+1) M_u,
X(v-v_[n]) <= q^(n+1) M_v+(n+1)s q^n M_u,
```

where `(u_[n],v_[n])=sum_(j=0)^n Lblock^j(KF,KF_r)`.
The second formula includes the differentiated coefficient insertion;
replacing it by a scalar q-only tail would omit a required term.

A second optional bound uses `v=u_r` to absorb `B_r v` into v's potential:
`Pv+Qu = A v_x+B v_r+(D+B_r)v+A_r u_x+D_r u`.
With certified cap `d_v>=sup|D+B_r|`, put
`q_v=T[a+b0+(T/2)d_v]`, `s0=T[ar+(T/2)dr]`.
If `q_v<1`, then `X(v)<=(Tg+s0 M_u)/(1-q_v)` and its free-source
tail is at most `q_v M_v+s0 M_u`. Either proven ledger, or the
intersection of both, is valid; this variant does not alter the PDE.

Before numeric evaluation the acceptance target should be: finite,
authenticated exact-rational enclosures for all six actual apex jet
entries; complete source differentiation and feedback accounting; exact
coordinate conversion; successful independent geometry/arithmetic and
manufactured identity controls. The finite partition/refinement budget
is declared separately by the evaluator. No derivative sign, narrow
relative error, or new physical invariant is promised. Intervals which
contain zero remain a successful bounded-two-jet result. If the declared
budget yields broad bounds, report their widths and the limiting input
uncertainties rather than changing the target after looking at signs.

## 9. Independent exact controls

`check_theory_controls.py` uses only Python's standard library and exact
rational polynomial algebra. For a polynomial forcing `x^i r^j`, it
computes the retarded map independently from the beta-integral formula

```
K[x^i r^j] = sum_(0<=k<=j, k even)
             binom(j,k) i! k!/(i+k+2)! x^(i+k+2) r^(j-k).
```

It checks the two characteristic derivative formulas, both second-jet
formulas, zero initial traces, a nonconstant-coefficient manufactured
u/v pair with all differentiated coefficients present, a nonzero initial
acceleration, the block-majorant power/resolvent identities, and the
natural-coordinate coefficient factors. Controls deliberately dropping
a B_r term, changing the characteristic minus sign, omitting the retarded
one-half factor, or omitting the radial conversion term must be rejected.
These checks test the identities and error-accounting structure without
mirroring the actual source numerical producer. They do not substitute
for authentication and validation of the actual derivative oracle.
