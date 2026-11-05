# A conditional metric-to-area-to-moment estimate

5 November 2026 (UTC). Working note.
This is standard determinant, pullback, and quotient algebra applied to the
specified observable. It is not a nonlinear field-equation estimate, an
evaluated physical amplitude radius, or a claim of optimal moment constants.

## 1. The full cut metric and the required domain

Retain the selected completion, intrinsic shell normal-geodesic observers,
transported initial labels, and common noncrossing proper-time cuts of the
frozen nonlinear bridge. On a certified compact domain
`|e|<=a<a_adm`, `tau in I`, `Omega in S3`, assume all pullbacks below are
twice continuously differentiable in `e`, uniformly. Parameter derivatives
hold **actual proper time and transported labels fixed**. Let `k_e` be the
full positive-definite three-dimensional metric of these cuts, pulled back
to the fixed label sphere, and let `h0=k_0`. Define

```
G_e=h0^(-1/2) k_e h0^(-1/2),   G_0=I,
j_e=sqrt(det G_e),            dA_e=j_e dA_0,
dA_0=c(tau)dOmega,            c(tau)>0.
```

This notation means a background-orthonormal frame, or equivalently the
corresponding endomorphisms on each tangent space. No global sphere frame
is required. Suppose, uniformly on this same domain,

```
lambda I <= G_e <= Lambda I,   lambda>0,
||G'_e||op <= P,              ||G''_e||op <= Q,
p=P/lambda,                  q=Q/lambda.
```

These are **fixed-background** operator budgets. `P,Q` include the
observer, clock, angular transport, and shell embedding. A sufficient way
to obtain the spectral bounds from anchoring alone is

```
a P<1,   lambda=1-aP,   Lambda=1+aP.                 (1)
```

Indeed `||G_e-I||<=|e|P` by integration. Alternatively use any independently
certified positive spectral gap and upper spectral bound; (1) is then not
required. No numerical values of these budgets are supplied here.

## 2. Exact determinant identities and a sharp coefficient

Put `A=k_e^(-1)k'_e`, `B=k_e^(-1)k''_e`. Their traces agree with the
background-matrix expressions. Their relative operator norms mean those
of the symmetric matrices `G^(-1/2)G'G^(-1/2)` and
`G^(-1/2)G''G^(-1/2)`, respectively; in particular they are at most `p,q`.
They do not mean the Euclidean norm of a possibly nonsymmetric matrix
representation of `G^(-1)G'`.

Jacobi's determinant identity and differentiation of the inverse give

```
ell=j'/j = (1/2)tr A,
beta=ell' = (1/2)tr B-(1/2)tr(A^2),
t=j''/j = (1/2)tr B+F(A),
F(A)=(1/4)(tr A)^2-(1/2)tr(A^2).                    (2)
```

For a self-adjoint endomorphism in dimension three,

```
-(5/4)||A||op^2 <= F(A) <= (3/4)||A||op^2.          (3)
```

Here is a proof over the entire eigenvalue cube, rather than a bare list
of corner values. Scale the norm to one and write the eigenvalues
`x=(x1,x2,x3) in [-1,1]^3`, `S=sum xi^2`, and
`F(x)=(sum xi)^2/4-S/2`. The upper bound follows from the exact identity

```
3/4-F(x)=(3-S)/4 + sum_{i<j}(xi-xj)^2/4 >=0.
```

For the lower bound define nonnegative weights
`w_v(x)=product_i(1+v_i xi)/2`, `v in {-1,1}^3`.
They sum to one and satisfy the exact interpolation identity

```
F(x)=sum_v w_v(x)F(v)+(3-S)/4.
```

At a vertex `sum vi` is `+/-3` or `+/-1`, so `F(v)` is respectively
`3/4` or `-5/4`. Thus `F(x)>=-5/4`. Equality occurs at
`(1,1,-1)` and permutations/sign reversal. Consequently **5/4 is the
sharp absolute quadratic coefficient**. The one-sided coefficient `3/4`
is also sharp. The zero-norm case is immediate. These facts concern the
quadratic determinant term, not sharpness of the later quotient estimate.

Writing `U=Lambda^(3/2)`, `m=lambda^(3/2)`, (2)-(3) give the usable plain
area-density interface

```
j>=m,        J0=sup j <=U,
J1=sup|j'| <=(3/2)U p,
J2=sup|j''|<=U[(3/2)q+(5/4)p^2].                    (4)
```

It remains valid without a smallness condition on `a J1`.

## 3. Charge every pullback term before using P and Q

One exact local chart description is useful for obtaining the input budgets.
Let `X_e(tau,Omega)` be the **full** cut map into a fixed spacetime chart:
it includes the moving shell embedding and intrinsic observer flow, proper
clock, and transported angle map. Set `Z=partial_e X`, `W=partial_e^2 X`,
`T=d_Omega X`, `T1=d_Omega Z`, `T2=d_Omega W`, and
`H=g_e(X_e)`. Then `k=T^T H T` and, with all spacetime derivatives taken
in this fixed chart,

```
H1=g_e,e + Dg_e[Z],
H2=g_e,ee +2Dg_e,e[Z]+D^2g_e[Z,Z]+Dg_e[W],
k'=T^T H1 T +T1^T H T+T^T H T1,
k''=T^T H2 T +2T1^T H1 T+2T^T H1 T1
     +T2^T H T+T^T H T2+2T1^T H T1.                (5)
```

Use positive-definite auxiliary chart norms for the Lorentzian spacetime
bilinear form. If component derivative bounds are `g_ij>=||partial_e^i
D^j g||` for `(i,j)=(0,0),(1,0),(0,1),(2,0),(1,1),(0,2)`, and
`|Z|<=z`, `|W|<=w`, `||T_r||<=L_r` for `r=0,1,2` (domain norm `h0`), put

```
h1=g10+g01 z,
h2=g20+2g11 z+g02 z^2+g01 w,
P := h1 L0^2+2g00 L0 L1,
Q := h2 L0^2+4h1 L0 L1+2g00(L0 L2+L1^2).          (6)
```

In (6), the displayed right sides may be **chosen as** the certified `P,Q`;
they bound the actual derivative norms. Use a fixed
finite atlas covering the common tube and charge all chart/norm comparison
constants before taking uniform maxima. The identities are local tensor
component identities. Differentiating only the coordinate metric misses
the `Z,W,T1,T2` terms. This note does not bound the observer ODE or prove
that its common noncrossing tube exists.

For the scalar pullback `A_e^sc=phi_e(X_e)-phi_b`, the corresponding exact
identities are

```
(A^sc)'=phi_e,e+Dphi_e[Z],
(A^sc)''=phi_e,ee+2Dphi_e,e[Z]+D^2phi_e[Z,Z]+Dphi_e[W]. (7)
```

Thus compatible scalar derivative budgets give
`E1<=f10+f01 z`, `E2<=f20+2f11 z+f02 z^2+f01 w`.
Again their right sides may be chosen as upper budgets. The notation
`A^sc` distinguishes this scalar from the metric endomorphism `A` in (2).
There is no uncharged coordinate-time replacement of (5) or (7).

## 4. A direct normalized moment bound

Use the same fixed transported template `Y=n1^2-n2^2` and scalar
`A^sc` as the quotient follow-on. Define

```
D=int Y^2 j dOmega,    N=int Y A^sc j dOmega,    M=N/D,
E_r=sup|partial_e^r A^sc|, r=0,1,2,
R=(int |Y|j dOmega)/D.
```

The factor `c(tau)` cancels before differentiation. The exact sphere
integrals are `int|Y|=2pi`, `int Y^2=pi^2/3`; hence

```
R <= (6/pi)chi,   chi = sup_{e,tau} max_Omega j/min_Omega j
                     <=(Lambda/lambda)^(3/2).       (8)
```

Let `nu` be the probability measure `Y^2j dOmega/D`, `mu=E_nu ell`, and
`nu_t=E_nu t`. Direct quotient differentiation groups all four terms into

```
M''=int (Yj/D){(A^sc)''+2(A^sc)'(ell-mu)
            +A^sc[(t-nu_t)-2mu(ell-mu)]} dOmega.   (9)
```

This form is regular at `Y=0` and needs no division by `Y`.
Since `|ell|<=3p/2`, `|ell-mu|<=3p`. A useful correlated estimate is

```
osc_Omega(t-2mu ell) <= 3q+9p^2.                   (10)
```

To prove (10), the `tr B/2` contribution has range width at most `3q`.
For the other contribution scale the eigenvalues of `A` by `p` and set
`b=mu/p in [-3/2,3/2]`. The polynomial is
`p^2 [F(x)-b sum xi]`. At `b=3/2` its derivative in `xi` equals
`-[(1+xi)+(1-xj)+(1-xk)]/2<=0`; its extreme values are therefore
`21/4` and `-15/4`, a width of `9`. At `b=-3/2` use sign reversal.
For intermediate `b`, the polynomial is a convex combination of these
two endpoint polynomials, whose oscillation is at most the corresponding
convex combination of their widths. If `p=0` this term vanishes.
Thus (10) holds even when eigenbases vary with angle. Subtracting a
probability mean costs at most the range width. Equations (8)-(10) imply

```
K(a)=sup|M''| <= (6/pi)chi [E2+6p E1+(3q+9p^2)E0]. (11)
```

Because `A^sc_0=0`, integration in the parameter gives `E0<=a E1`.
Consequently

```
K(a) <= (6/pi)chi [E2+6p E1+a E1(3q+9p^2)].        (12)
```

This avoids repeated absolute-density denominator factors and is valid
under the spectral gap alone. It does not require `aJ1<=1/2`. One can
instead insert (4) into the older general quotient inequality; that is a
valid, usually looser interface. No claim compares (12) with every possible
specialized prior bound for every possible input budget.

## 5. Optional angular-oscillation refinement

Uniform angular rescaling of area should cancel. The exact normalized
density `rho=j/D` satisfies, with `x=ell-E_nu ell`,

```
rho'/rho=x,
rho''/rho=x^2-Var_nu(ell)+beta-E_nu beta.            (13)
```

Indeed differentiating `E_nu ell` gives `E_nu beta+Var_nu(ell)`.
Suppose stronger inputs certify
`osc_Omega ell<=delta` and `osc_Omega beta<=eta` uniformly. Since
`|x|<=delta`, `Var<=delta^2/4`, and
`|x^2-Var|<=delta^2`, (13) gives

```
K(a)<=(6/pi)chi [E2+2delta E1+(delta^2+eta)E0].     (14)
```

The variance bound follows from
`E[(ell-l)(u-ell)]>=0`, followed by
`(mean-l)(u-mean)<=(u-l)^2/4`; it is not essential to the looser
`|x^2-Var|<=delta^2` step. If `ell,beta` are angular constants their
contributions in (13) vanish exactly. Metric budgets alone imply
`delta<=3p`, `eta<=3q+3p^2/2`, but (11) is better than substituting
these generic values into (14). Use the smaller of (11) and (14) when
the respective hypotheses are certified.

Since `log j_0=0`, integration of angular differences also proves
`chi<=exp(a delta)`. Either this or (8) may bound the contrast. The
exponential estimate can be more useful when the metric includes a large
angularly uniform scale factor. All score budgets refer to the **same**
full pullback and admitted interval.

## 6. Exact sufficient guards, and what remains missing

For rational supplied budgets choose a positive rational `C` with

```
C^2 lambda^3 >= Lambda^3.                          (15)
```

This certifies `C>=(Lambda/lambda)^(3/2)` without rounding. With `p=P/lambda`,
`q=Q/lambda`, define the wholly rational bound

```
H=2C [E2+6p E1+a E1(3q+9p^2)].                    (16)
```

Here `6/pi<2`, with a weak `K<=H` also when all budgets vanish. After
independently certifying the common observer domain, spectral bounds,
scalar budgets and first coefficient `M'_0=S`, the guard

```
0<a<a_adm,       a H<106                            (17)
```

Use the retained actual proper-time interval
`I=[-rho_b*6145*ell_width/4096,-rho_b*6143*ell_width/4096]`, with its
background radius `rho_b>0` and retained width `ell_width=10^-6` (distinct
from the score `ell` in (2)). The guard preserves the retained
endpoint signs `S(t_L)>63`, `S(t_R)<-53` for
every `0<|e|<=a`: integral Taylor gives
`sup|M_e/e-S|<=aK/2<=aH/2<53`. The intermediate value theorem then
gives a scalar **moment** zero between those same actual proper-time
endpoints. It does not give a full scalar field zero or a unique root.

For (14), replace the bracket in (16) by
`E2+2delta E1+a E1(delta^2+eta)`. An optional rational certificate
for the contrast alternative is `a delta<n`, `n` a positive integer,
and `C>=(1-a delta/n)^(-n)`: the latter exceeds `exp(a delta)` by
integrating `1/(1-z)>=1` on `[0,a delta/n]`. These conditions are
sufficient, not necessary.

The accompanying standard-library Fraction checker validates the finite
identities, determinant jet formula including noncommuting rational
matrices, sharpness witnesses, correlated oscillation certificates,
normalized quotient identities, pullback coefficient bookkeeping, and
strict guard rejection. It runs identically under normal Python and
`python -O` and uses no assertions. Its supplied guard example is
**synthetic algebra**, never an evaluated HDBLAST amplitude. It does not
validate any missing analytic budget or observer domain.

The retained C2 theory gives no numerical `P,Q,E1,E2,a_adm`, spectral gap
or observer-flow constant here. Selected-family oddness `M''_0=0` remains
compatible with these absolute bounds but supplies no cubic remainder
under C2 regularity. The particular nonlinear completion, PDE/evolution
bounds, preparation and empirical interpretation remain separate inputs.

Sources read without modification: the private
`hdblast-observable-quotient-followon-04h/QUOTIENT_BOUND.md` and the public
`nonlinear/NONLINEAR_BRIDGE.md`, `nonlinear/GAP_LEDGER.md`, and
`nonlinear/inputs/INTRINSIC_OBSERVER_AND_SIGNS.md` in the retained odd-trace
tree. `SOURCE_AUTHENTICATION.json` pins the exact bytes observed. The
determinant, chain, quotient, variance, Taylor, and intermediate-value
arguments above are standard mathematics; the contribution is this explicit
conditional bookkeeping interface, not a new general theorem.
