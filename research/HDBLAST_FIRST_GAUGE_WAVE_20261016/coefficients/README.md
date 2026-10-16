# Whole-collar coefficients and actual gauge-wave source

16 October 2026. This directory certifies the matrices and source of the
raw three-amplitude first relative-harmonic gauge system. It uses Python's
standard library and exact rational/outward interval arithmetic. It does
not import any previous scientific program or the annulus-only jet report.

The verified background rectangle is

\[
0\le z\le8\ell,\qquad0\le t\le3\ell,\qquad\ell=10^{-6}.
\]

The associated public vector-wave solution theorem can use only the
backward-closed scaled triangle
\(D=\{0\le x\le3,\ r\ge0,\ x+r\le8\}\), where
\(x=t/\ell\), \(r=z/\ell\). No rectangle extension of the master
theorem is needed for that application. The additional source coverage on
the full rectangle uses its inherited smooth finite propagation, as
explained below.

The exact bound report is `GLOBAL_GAUGE_COEFFICIENTS_NORMAL.json`. The
optimized report has identical bytes. Its strict readable bounds include:

| Quantity | Strict upper bound |
| --- | ---: |
| Maximum absolute row sum of A | 0.003761 |
| Maximum absolute row sum of B | 0.00007236 |
| Maximum absolute row sum of D | 0.00000003177 |
| Absolute scaled source S | 0.0000007392 |
| Every first scaled source derivative | 0.00001183 |
| Every second scaled source derivative | 0.000308 |
| Weighted C1 contraction row sum | 0.011412 |
| Absolute deviation of P_s from 6ell/130 | 0.00000000007062 |
| Absolute P_0 | 0.00000000001280 |

All decimals in this table are exact rational numbers, checked strictly
against the computed Fraction bounds. Elementwise matrix bounds and every
mixed coefficient derivative through order two are retained in the report.

## 1. Actual global background and the collar-length issue

Use the authenticated actual potential, written in \(\eta=\phi-1\),

\[
U(\eta)=-\tfrac2{27}+\tfrac{14}9\eta^2+\tfrac{50}{27}\eta^3
        -\tfrac16\eta^4-\tfrac49\eta^5-\tfrac2{27}\eta^6.
\]

The shell state and retained outer box are pinned copies of the registered
branch at \(\delta=1/1000\), with its exact registered value of the
tension parameter. The actual static ODE in inward tortoise depth is

\[
\rho_z=-\rho^2a,\quad \eta_z=-\rho p,\quad
p_z=-\rho(U_{,\phi}-4ap),\quad
a_z=\rho^{-1}+\rho p^2/3.
\tag{1}
\]

The retained global outer rectangle is

\[
129\le\rho\le131,\quad-9/10^5\le\eta\le-1/12500,\quad
-7/50000\le p\le-3/25000,\quad111/1000\le a\le14/125.
\]

Let \(I\) be this box, \(I_b\) the authenticated shell box, and
\(Z=8\ell\). The checker independently evaluates the vector-field
interval \(F(I)\) and proves

\[
J=I_b+[0,Z]F(I)\subset\operatorname{interior}(I).
\tag{2}
\]

Its exact four pairs of face margins equal the separate fresh-referee
certificate. They exceed the following positive rational numbers:

| Component | Lower margin | Upper margin |
| --- | ---: | ---: |
| rho | 9/10 | 107/100 |
| eta | 27/10^7 | 69/10^8 |
| p | 28/10^7 | 41/10^7 |
| a | 377/10^6 | 622/10^6 |

For a smooth static ODE, (2) is a first-exit/continuation proof over the
entire interval \([0,Z]\). Before a first exit, the integral ODE puts the
actual state in J, strictly inside I; hence such an exit is impossible.
The compact image and smooth vector field permit continuation to Z.
The original actual branch and shell-box authentication remain retained
hypotheses. This does not independently reprove their existence.

In particular, the older prose's short collar of length \(8\times10^{-10}\)
is not treated as if it already covered \(8\times10^{-6}\). The larger
interval is justified by (2), which is replayed here and compared with the
fresh independent referee margins. The narrower tube J, which covers the
whole collar, is then used to sharpen the new derivative bounds. No
initial-annulus estimate is extrapolated.

## 2. Global alpha and normalized derivative jets

Define

\[
h=U_{,\phi}/p,\quad k=\rho(5a/2-h),\quad
\alpha_z=k\alpha,\quad\alpha(0)=1,\quad C_{\rm rec}=\alpha/\rho^2.
\tag{3}
\]

The reconstruction coefficient \(C_{\rm rec}\), called `C` in this
report, is not the older shell report's scalar master drift named C.
The nonzero p gap is certified on the whole collar.

The interval image of k on J is negative. With \(K=\sup_J|k|\),

\[
1-KZ\le e^{-KZ}\le\alpha(z)\le1,\qquad0\le z\le Z.
\tag{4}
\]

The checker proves \(1-KZ>0\). Equation (4) integrates from the actual
normalization point zero over the full collar, including the region before
the initial packet annulus. It does not reuse an annulus alpha upper bound
that would exclude \(\alpha(0)=1\).

At any point \(z_0\), let \(j_n=D_z^nf(z_0)/n!\). Products use finite
convolution and inverses use the exact reciprocal recursion. For each
static state variable, the normalized ODE recurrence is

\[
j_{n+1}=\frac{[s^n]F(j_0+j_1s+\cdots+j_ns^n)}{n+1}.
\tag{5}
\]

This is automatic differentiation at every allowed point, through order
three. It is not a sampled trajectory, a Taylor approximation to the
solution, or a claim about an unbounded Taylor remainder. The zeroth jets
range over J and (4). Interval operations then enclose all normalized
derivatives. The proper-radial derivative needed by the operator is

\[
a_y=-\rho^{-2}-p^2/3,
\]

and its z derivatives are obtained by differentiating this identity.
The report includes \(\rho,a,a_y,\alpha,C_{\rm rec},C_{{\rm rec},z}\),
the other static auxiliaries, and their available derivative jets.
Reported `static_z_jets` are normalized by n!; the accompanying
`static_actual_derivative_caps` restore n!.

An independent sparse Laurent-polynomial total derivative verifies every
static jet identity through order three. It uses the exact ODE vector
field and does not rely on interval overlap as evidence of equality.
After the exact first-exit check, each interval operation is rounded
outward on a rational grid of denominator \(10^{36}\). Both endpoints
are rounded by integer floor/ceiling. No floating arithmetic is used.

## 3. Time factors on the full time interval

Put \(H=\tanh t\), \(q=\cosh(t)^{-3/2}\), and \(T=3\ell\).
Elementary differential inequalities give

\[
0\le H\le T,\quad1-T^2\le\operatorname{sech}^2t\le1,\qquad
1-\tfrac34T^2\le q\le1.
\tag{6}
\]

For the last inequality,
\(\log q=-\frac32\int_0^t\tanh s\,ds\ge-3T^2/4\), followed by
\(e^{-v}\ge1-v\). The normalized jet recurrences for
\(H_t=1-H^2\) and \(q_t=-3Hq/2\) enclose every time derivative used.
Thus no transcendental floating-point evaluation is needed.

## 4. Matrix and actual source interface

The operator agent's raw amplitudes are
\(U=(f,b,d)=(\xi^y/Y,\xi^t/Y,\hbox{angular gradient amplitude})\).
There is no extra rho rescaling. The convention is

\[
U_{xx}-U_{rr}=A U_x+B U_r+D U+(S,0,0).
\tag{7}
\]

The checker evaluates the matrices supplied by the explicit gauge-vector
operator with \(\lambda=8\). It reports actual mixed scaled derivatives
\(D_x^iD_r^j\) through \(i+j=2\), including all ell factors, and
their elementwise and row-sum bounds.

For the known first master, put \(W(x,r)=w(\ell x,\ell r)\). The
actual longitudinal amplitude is
\(u=C_{\rm rec}(z)q(t)w(t,z)\). Since
\(D_y=-\rho^{-1}D_z\), the radial gauge source is exactly

\[
-6\ell^2\rho^2(u_y+4au)
=P_s W_r+P_0W,
\]
\[
P_s=6\ell\rho C_{\rm rec}q,\qquad
P_0=6\ell^2\rho(C_{{\rm rec},z}-4\rho a C_{\rm rec})q.
\tag{8}
\]

The source sign, derivative and ell factors are independently checked as
exact Laurent-polynomial identities. Equivalently,
\(C_{{\rm rec},z}=\rho(9a/2-h)C_{\rm rec}\), so
\(P_0=6\ell^2\rho^2C_{\rm rec}(a/2-h)q\).

The pinned master certificate supplies strict uniform mixed derivative
caps by total scaled order zero through three:

\[
(M_0,M_1,M_2,M_3)
=(1.000193,16.000488,256.002081,6665.267480).
\]

The checker verifies each decimal exceeds the report's exact Fraction
bound. The derivative source bounds follow from the full two-variable
Leibniz rule; a source derivative of order j needs master order j+1.
This supplies every mixed source derivative through order two without
assuming any new master derivative theorem.

The master bounds were proved on D. On the larger rectangle the actual
packet's finite propagation gives support \(z\le2\ell+t\). Wherever
the packet or any of its smooth derivatives is nonzero,
\(z+t\le2\ell+2t\le8\ell\). All remaining points have zero jets,
including the smooth flat support front. This is the precise additional
support argument for the rectangle source bounds; the vector solution
construction can stay on D.

## 5. Conditional Volterra row bounds and scope

For the separate mixed Dirichlet/Neumann solution theorem, use

\[
N_i=\max\{(2/3)\|U_i\|_\infty,
                  \|U_{i,x}\|_\infty,\|U_{i,r}\|_\infty\}.
\]

If the retarded kernels have weighted norm at most 3, their coefficient
matrix majorant is

\[
M_{ij}=3(\sup|A_{ij}|+\sup|B_{ij}|+\tfrac32\sup|D_{ij}|).
\]

Its row sum is less than 0.011412. Exact rational inversion of I-M is
nonnegative and bounds the three channels by
\(3(I-M)^{-1}_{i1}\sup|S|\). These are conditional algebraic outputs;
the kernel bound, solution construction, boundary compatibility and
regularity are proved separately by the solution agent. This coefficient
certificate alone is not an evolved gauge vector, second-order gauge
conversion, or full coupled shell-contact theorem.

## 6. Portable replay and controls

All inputs are copied under `inputs/`, hashed in the checker, and traced
with relative `sourcePath` entries in `INPUT_PROVENANCE.json`. Programs
write only to a new explicit `--output` outside this coefficient input
tree. They reject overwriting. For example, from this directory, after
creating a new external output directory:

```sh
python -I -B check_global_gauge_coefficients.py --output /tmp/new-coefficient-replay/normal.json
python -I -B -O check_global_gauge_coefficients.py --output /tmp/new-coefficient-replay/optimized.json
python -I -B run_coefficient_controls.py --output /tmp/new-coefficient-replay/controls.json
```

Normal and optimized bound reports are byte-identical. The control runner
executes five semantic corruptions in each Python mode: wrong Taylor
factorial, wrong radial-source sign, annulus-length alpha integral, alpha
interval excluding its origin normalization, and omitted ell in A. Each
must fail for its named reason with an explicit exception, including under
optimization. The complete roster is `COEFFICIENT_REPLAY_ROSTER.json`.
The ordinary/optimized reports, control receipt and file hashes are
recorded in `COEFFICIENT_VERIFICATION_RECEIPT.json`.
