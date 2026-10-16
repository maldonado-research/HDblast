# Actual first gauge wave: causal construction and an explicit certified comparator

16 October 2026. Conditional mathematical continuation of the retained actual
regular background and actual registered first master. The causal construction
below supplies their actual first gauge wave, an explicit evaluated comparator
and a continuous rational defect certificate. All work is separate from
previous frozen checkpoints.

On the full triangle, the actual raw amplitudes obey

| Amplitude | Uniform value | Each scaled first derivative |
|---|---:|---:|
| \(f\) | \(7.47\times10^{-8}\) | \(1.14\times10^{-6}\) |
| \(b\) | \(1.14\times10^{-12}\) | \(7.59\times10^{-13}\) |
| \(d\) | \(7.49\times10^{-20}\) | \(5.00\times10^{-20}\) |

Every displayed upper bound is strict. The radial comparator's uniform
value error is less than \(5.47\times10^{-9}\); its first-derivative
error in the weighted C1 norm is less than \(3.65\times10^{-9}\).
These are raw coordinate-amplitude bounds, with x/r derivatives in scaled
coordinates; t derivatives divide by \(\ell\), y derivatives by
\(\ell\rho\). They are not a covariant metric or scalar forcing norm.

Retain the actual smooth regular background and the actual registered first
master solution. Put \(x=t/\ell\), \(r=z/\ell\),
\(\ell=10^{-6}\), and
\(D_X=\{0\le x\le X,\ r\ge0,\ x+r\le8\}\), \(X\le3\).
The vector amplitudes are **raw contravariant** components

\[
\xi^y=fY,\qquad \xi^t=bY,\qquad
\xi^i=dD_\kappa^iY,\qquad Y=n_1^2-n_2^2.
\]

No radial amplitude is divided by \(\rho\). The saved operator gives

\[
U_{xx}-U_{rr}=A U_x+B U_r+D U+(S,0,0),\qquad U=(f,b,d),
\tag{1}
\]

with zero displacement and velocity at \(x=0\),
\(f(x,0)=0\), and \(b_r(x,0)=d_r(x,0)=0\). In particular the
tangential values \(b(x,0),d(x,0)\) are not set to zero. The authenticated
coefficient certificate supplies all entries of the scaled matrices. Its
source is the actual function

\[
S=P_s W_r+P_0W,\quad W(x,r)=w(\ell x,\ell r),
\]

\[
P_s=6\ell\rho C\cosh(\ell x)^{-3/2},\qquad
P_0=6\ell^2\rho^2C(a/2-h)\cosh(\ell x)^{-3/2},\quad
C=\alpha/\rho^2,\quad h=\mathcal V'(\phi_0)/p.
\tag{2}
\]

Here \(\mathcal V\) denotes the unchanged registered bulk potential,
called \(U\) in the inherited scalar model; \(U\) here is the vector.
Intervals are used to bound this actual master, not to replace it by unrelated
pointwise selections of source values.

## Causal Dirichlet/Neumann Volterra construction

For a continuous source \(F\) define the original half-line retarded kernels

\[
(K_{D/N}F)(x,r)=\frac12\int_0^x\int_0^\infty
\left[\mathbf1_{|r-s|\le x-v}\mp\mathbf1_{r+s\le x-v}\right]
F(v,s)\,ds\,dv.
\tag{3}
\]

The upper sign is Dirichlet. The second indicator's support is a subset of
the first. For either kernel,

\[
\|KF\|_\infty\le X^2\|F\|_\infty/2,\qquad
\|(KF)_x\|_\infty,\|(KF)_r\|_\infty\le X\|F\|_\infty.
\tag{4}
\]

The derivative estimates follow from the two direct characteristic integrals
and their boundary joins; coefficient functions are never reflected and
differentiated as if they had constant parity. \(K_DF\) has zero boundary
value and \(K_NF\) zero radial boundary derivative. Both have zero initial
displacement and velocity. Every effective predecessor obeys
\(v+s\le x+r\le8\), so (3) closes on the audited triangle, with no
artificial inner boundary or coefficient extension beyond its collar.

Use channel norms

\[
N_i(U)=\max\{2\|U_i\|_\infty/X,\|(U_i)_x\|_\infty,
                                      \|(U_i)_r\|_\infty\}.
\]

Let \(a_{ij},b_{ij},d_{ij}\) be nonnegative absolute coefficient caps and set

\[
M_{ij}(X)=X(a_{ij}+b_{ij}+Xd_{ij}/2),\qquad
q_X=\max_i\sum_jM_{ij}(X).
\tag{5}
\]

The operator \(\mathcal T(U)=K(AU_x+BU_r+DU+(S,0,0))\), with
\(K=\operatorname{diag}(K_D,K_N,K_N)\), is a contraction on the
closed C1 space with these initial/boundary rows when \(q_X<1\).
The global certificate gives \(q_3<0.011412\). Therefore

\[
U=\sum_{n\ge0}(K\mathcal L)^nK(S,0,0),\quad
\mathcal L U=AU_x+BU_r+DU,
\tag{6}
\]

converges in C1 and uniquely constructs the actual forced vector on the
whole observation triangle. The matrix resolvent
\((I-M)^{-1}=\sum_{n\ge0}M^n\) is nonnegative. It preserves the useful
channel distinction instead of bounding all three amplitudes by the radial
channel's much larger size.

The fixed point is initially a C1 mild solution. The authenticated separate
regularity argument proves that it is a classical C3 solution using the
actual compatible C2 coefficients and source, time-differentiated contractions,
and the original physical boundary trace terms. This is an analytic regularity
result, not a numerical bound for its C2/C3 jets. A C1 norm estimate alone is
not differentiated to assert higher derivatives, and no differentiated odd
extension of a nonzero Dirichlet source trace is used.

## An explicit comparator, including the characteristic crossing

Let \(Q(r)=q_\ell(\ell r)\), supported in \([1,2]\), and define two
different odd initial profiles

\[
Q_o(s)=\operatorname{sgn}(s)Q(|s|),\qquad
P_o(s)=sQ(|s|).
\]

The second is the odd extension of \(rQ(r)\). It is **not**
\(sQ_o(s)\), which would be even. Both profiles are smooth and flat near
zero. Put

\[
W_{\rm disp}=\tfrac12[Q_o(r+x)+Q_o(r-x)],\qquad
c_*={6\ell\over130},
\]

\[
\boxed{p_f={c_*\over4}\{P_o(r+x)+P_o(r-x)
                   -r[Q_o(r+x)+Q_o(r-x)]\},\quad p_b=p_d=0.}
\tag{7}
\]

This comparator has no quadrature or unknown boundary function. If \(V\) is
the free Dirichlet wave with initial profile \(rQ\) and zero velocity,
then \(p_f=c_*(V-rW_{\rm disp})/2\). Since
\((\partial_x^2-\partial_r^2)(rW_{\rm disp})=-2\partial_rW_{\rm disp}\),

\[
(p_f)_{xx}-(p_f)_{rr}=c_*\partial_rW_{\rm disp}.
\tag{8}
\]

Its initial value and velocity vanish on \(r\ge0\), and its Dirichlet
boundary value is zero. We use these auxiliary free profiles only to evaluate
(7) on the physical half-line; no parity is assigned to the actual variable
coefficients or coupled source.

For \(r\ge x\),
\(p_f=c_*x[Q(r+x)-Q(r-x)]/4\). For \(0\le r<x\),

\[
p_f={c_*\over4}[xQ(r+x)+(2r-x)Q(x-r)].
\]

These agree smoothly at \(r=x\): the potentially problematic endpoint is
zero, where all bump jets vanish. Since \(|Q|\le1\), \(|Q'|\le16\),

\[
\|p_f\|\le c_*X/2,\qquad
\|(p_f)_x\|,\|(p_f)_r\|\le c_*(1/2+8X).
\tag{9}
\]

At the boundary \((p_f)_{rr}=-c_*Q'(x)\) in general; one must not
pretend that this Dirichlet solution or its source has a globally smooth
odd extension at later boundary times.

## Actual-master defect bound and evaluated signs

The retained master has the exact decomposition into the zero-velocity free
wave above, its odd initial-velocity wave, the **actual** boundary lift,
and its potential Duhamel term. The scaled initial velocity has supremum
at most \(\ell\). Its radial derivative contribution is bounded directly
by the two endpoint velocities, hence by \(\ell\). If \(B_{y'}\)
bounds the actual boundary lift derivative and \(B_W\) the actual master,

\[
\|W_r-(W_{\rm disp})_r\|_{D_X}
\le\ell+B_{y'}+X\ell^2(102000)B_W.
\tag{10}
\]

For \(X\le1\), finite propagation gives the actual boundary function and
all its derivatives zero through the front, so \(B_{y'}\) is replaced by
zero. For \(X>1\) the authenticated actual-boundary derivative cap is
retained. Equation (10) is a derivative identity estimate, not a differentiated
value bound.

Let \(\Delta_s\ge|P_s-c_*|\), \(S_*\ge|P_s|\),
\(P_*\ge|P_0|\). The source defect is bounded by

\[
R_s=16\Delta_s+S_*[\ell+B_{y'}+X\ell^2(102000)B_W]+P_*B_W.
\]

Writing \(p_0=c_*X/2\), \(p_1=c_*(1/2+8X)\), define

\[
R_i=\delta_{i1}R_s+(a_{i1}+b_{i1})p_1+d_{i1}p_0.
\]

Then the exact defect certificate is

\[
\boxed{N_i(U-p)\le E_i=[(I-M(X))^{-1}XR]_i,
\quad\|U_i-p_i\|_\infty\le XE_i/2.}
\tag{11}
\]

This also gives point enclosures by adding the signed explicit value (7).
For example, (7) evaluates exactly to \(+c_*/8\) at
\((x,r)=(1/2,1)\) and \(-c_*/8\) at \((1/2,2)\); these evaluations
use only \(Q(3/2)=1\) and zero values outside the support. The small
local-time defect proves opposite signs for the **actual radial gauge-vector
amplitude** at these two points. The final exact receipt records the rational
error margins and additional evaluated points, including one on \(r=x\).
For three convenient strict actual enclosures, all endpoints below are exact
terminating rational decimals:

| Scaled \((x,r)\) | Actual \(f(x,r)\) interval |
|---|---:|
| \((1/2,1)\) | \((5.6,6.0)\times10^{-9}\) |
| \((1/2,2)\) | \((-6.0,-5.6)\times10^{-9}\) |
| \((3/4,3/4)\) | \((8.3,9.0)\times10^{-9}\) |

The third lies on the characteristic crossing \(r=x\). The receipt checks
each displayed endpoint against the full rational actual interval, rather
than deciding a sign with a float. These are coordinate gauge-vector
amplitudes; they do not describe an invariant force or a second-order response.

Componentwise Picard refinement from this comparator has error cone
\(M(X)^nE\) after n further applications of \(\mathcal T\), and the
scalar bound \(q_X^n\max_iE_i\) follows. Thus the representation has a
quantitative stopping criterion even before any finer numerical integral
evaluation is added.

## Authentication, exact checks and scope

`INPUT_PINS.json` authenticates the whole-collar coefficient/source report,
the actual master/boundary-jet report, and the operator and regularity
derivations. The coefficient report's explicit master hash is checked against
the separately authenticated master input. The old smaller shell-box collar
is not used as if it covered this triangle: the coefficient report includes
the fresh whole-\(8\ell\) first-exit proof and its provenance.

`check_actual_gauge_solution.py` uses only the Python standard library. It
evaluates every matrix inverse, error cone, point sign and stopping constant
with exact Fractions; floats are display-only. Forty-two finite exact checks
exercise the physical D/N kernel masses, comparator Cauchy and Dirichlet
rows, and its commutator source, including \(r=x\). Four controls flip a
kernel image sign, replace the odd product profile by the wrong even one,
or halve the commutator factor. Both ordinary and optimized Python must
reject each for its intended reason and leave no output report.

```sh
python check_actual_gauge_solution.py --run-controls --output ACTUAL_GAUGE_SOLUTION_NORMAL.json
python -O check_actual_gauge_solution.py --run-controls --output ACTUAL_GAUGE_SOLUTION_OPTIMIZED.json
```

Paths are relative to the checker directory by default, or to `--input-root`;
`--output` may select an external path. No prior scientific producer module
is imported or executed. The operator, Green and contraction arguments give
the all-point proof; finite witnesses check specified identities only.

The construction solves the **first** gauge-vector wave for the actual
registered source on the controlled finite triangle. It does not evaluate
the transformed first metric/scalar fields or the transformed quadratic bulk
and physical boundary forcing, solve second evolution, certify a nonlinear
finite-amplitude family or remainder, or compute an invariant force. Those
steps must retain tangential shell/clock transport and every quadratic
covariance term in the same representative. Neumann tangential rows do not
mean zero tangential gauge-vector values.
