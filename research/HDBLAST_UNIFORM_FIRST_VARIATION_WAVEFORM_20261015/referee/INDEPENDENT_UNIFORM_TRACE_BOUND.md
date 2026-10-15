# Independent uniform waveform-error bound for the same registered state

15 October 2026 UTC. Internal independent analytic review. The retained actual
background and exact-family assumptions remain conditional inputs. This note
does not deduce a pointwise bound from an `L2` error estimate.

## 1. Fixed state and exact output

Retain the preceding response checkpoint's exact-decimal near-`phi=+1`
background at `delta=1/1000`, and the same normalized bump width
`ell=1/10^10`. The conjugated master has initial data
`w(0,z)=q_ell(z)`, `w_t(0,z)=exp(-z)q_ell(z)`, and obeys

\[
w_{tt}-w_{zz}+Ww=0,\qquad
-w_z+\mu w_{tt}+\beta(t)w=0\quad(z=0),\qquad\mu>0.
\tag{1}
\]

Write `gamma=C_b/2` and `D(t)=w_z(t,0)-gamma w(t,0)`. The physical
linear scalar harmonic amplitude is

\[
S(t)=\frac{3\cosh^{-3/2}t}{\rho_b^3p_b}D(t).
\tag{2}
\]

Use `T=3ell` and the shrinking domain `0<=z<=L-t`, `L=8ell>T`.
The retained actual-background enclosure supplies, throughout that domain,

\[
|W|\le M=10^6,\quad |\beta|\le B=1000,\quad
\mu^{-1}\le a=1000,\quad |\gamma|\le250.
\tag{3}
\]

The prior analytic collar JSON is pinned at
`8ccba65a1c7fd66e175703cd4dfb8e1adb029b52d455a2b9a440f346bf2e4bc8`.
These numerical bounds concern the actual analytic branch, not sampled
historical floating profiles. The full nonlinear constraint correction
and observation-time margins remain separately required.

## 2. Comparator, sources and characteristic formulas

Let `u` be the same smooth free Dirichlet comparator from the previous
quantitative checkpoint: extend `q_ell` and `exp(-z)q_ell` oddly across
zero and use the whole-line d'Alembert formula. Then

\[
u(t,0)=u_{tt}(t,0)=0,\quad
\widehat D(t)=u_z(t,0)=q_\ell'(t)+e^{-t}q_\ell(t),
\quad |u(t,z)|\le1+T.
\tag{4}
\]

For `e=w-u`, all initial error data are zero, and

\[
e_{tt}-e_{zz}+We=F=-Wu,\qquad
-e_z+\mu e_{tt}+\beta e=G=\widehat D.
\tag{5}
\]

In particular the dynamic boundary defect is present with the **positive**
source sign shown. Safe uniform source bounds are

\[
|F|\le F_0=M(1+T),\qquad |G|\le G_0=16/\ell+1.
\tag{6}
\]

Put `j_+=e_t+e_z`, `j_-=e_t-e_z`, `q_b=e(t,0)`, and `v=q_b'`.
Their exact equations are

\[
(\partial_t\mp\partial_z)j_\pm=F-We,
\qquad \mu v'+v+\beta q_b=j_+(t,0)+G,
\quad j_-(t,0)=2v-j_+(t,0).
\tag{7}
\]

Let `E` and `V` denote the actual global suprema of `|e|` on the domain
and `|v|` on `[0,T]`. They are finite for the retained smooth solution.
Backward incoming characteristics reach the initial edge and give

\[
|j_+(t,z)|\le t(F_0+ME).
\tag{8}
\]

An outgoing characteristic either reaches the initial edge or the shell.
In the shell case, combine the incoming integral up to the reflection time
with the outgoing source integral after reflection. Those time intervals
are disjoint, so

\[
|j_-(t,z)|\le 2V+t(F_0+ME).
\tag{9}
\]

Thus `|e_t(t,z)|<=(V+t(F0+ME))`. Integrating vertically at fixed `z`
stays inside the shrinking domain and yields

\[
E\le TV+\frac{T^2}{2}(F_0+ME).
\tag{10}
\]

The damped boundary variation-of-constants formula, with zero initial
velocity, is

\[
v(t)=\mu^{-1}\int_0^t e^{-(t-s)/\mu}
[j_+(s,0)+G(s)-\beta(s)q_b(s)]ds.
\]

Use its positive kernel, `|q_b|<=E`, `mu^-1<=a`, and (8). Discarding
damping gives the conservative inequality

\[
V\le aT[G_0+TF_0+(TM+B)E].
\tag{11}
\]

No trace theorem, residual differentiation or boundary acceleration
estimate is hidden in (8)-(11). No artificial outer-boundary condition
is imposed: every characteristic used above stays in the stated domain.

## 3. Explicit closed bound and exact-rational consequence

Define

\[
\theta=aT^2(TM+B)+MT^2/2.
\tag{12}
\]

If `theta<1`, substituting (11) into (10) gives the explicit bounds

\[
E\le E_*:=\frac{T^2[aG_0+(aT+1/2)F_0]}{1-\theta},
\quad
V\le V_*:=aT[G_0+TF_0+(TM+B)E_*].
\tag{13}
\]

Since `e_z=j_+-v` and `u_b=0`, the actual observed derivative trace obeys

\[
\boxed{\|\ell(D-\widehat D)\|_{L^\infty(0,T)}
\le\ell[T(F_0+ME_*)+V_*+250E_*]
<\frac1{200000}=5\,10^{-6}.}
\tag{14}
\]

The independent program `independent_uniform_trace_audit.py` evaluates only
these rational inequalities, authenticates the prior collar JSON, and checks
the required coefficient gates. For reference, the rational `theta` is
`1.35000027e-13` when displayed as a decimal; the computed uniform error
is approximately `4.800000390030874e-6`. The rational values and strict
comparisons, not these decimal displays, decide pass/fail.

### A sharper bound using the known signed boundary forcing

The known source `G=q_ell'+exp(-t)q_ell` permits an exact integration by
parts inside the damped boundary formula. Since `q_ell(0)=0` and
`0<=q_ell<=1`,

\[
\mu^{-1}\int_0^t e^{-(t-s)/\mu}q_\ell'(s)ds
=q_\ell(t)/\mu-\mu^{-2}\int_0^t e^{-(t-s)/\mu}q_\ell(s)ds.
\tag{14a}
\]

Both nonnegative terms on the right lie between zero and `1/mu`.
Their difference therefore has absolute value at most `1/mu`, without
any estimate involving `sup|q_ell'|`. The remaining initial velocity term
has convolution magnitude at most `T/mu`. Using `|q_b(s)|<=sV` instead
of `|q_b|<=E` gives

\[
(1-s_*)V\le aU+aT^2(F_0+ME)/2,\quad
U=1+T,\ r=MT^2/2,\ s_*=aBT^2/2.
\tag{14b}
\]

Together with (10) and `F0=MU`, this is the positive two-variable system

\[
(1-r)E-TV\le rU,\qquad
-arE+(1-s_*)V\le aU(1+r).
\tag{14c}
\]

For `r<1`, `s_*<1`, and
`d=(1-r)(1-s_*)-arT>0`, its inverse has nonnegative entries. Hence

\[
E\le E_\sharp=
\frac{U[r(1-s_*)+aT(1+r)]}{d},\qquad
V\le V_\sharp=\frac{aU}{d}.
\tag{14d}
\]

In particular `U+E_sharp=U(1-s_*+aT)/d`. The observed error now has the
sharper bound

\[
\boxed{\|\ell(D-\widehat D)\|_\infty\le
\eta_\sharp:=\ell\{MT(U+E_\sharp)
                    +(1+250T)V_\sharp\}
<1.0000004\,10^{-7}.}
\tag{14e}
\]

The independent rational audit checks all the positive-inverse gates and
this strict outward bound as well. It obtains the same sharp determinant
formula as the separately derived direct-master route through an error-field
characteristic argument. The cancellation is permitted because `mu` is a
positive constant and the source is the explicitly specified smooth bump;
it is not valid for an arbitrary source merely bounded by `G0`.

## 4. Continuous waveform interpretation and cautions

Set `x=t/ell` and `Q(x)=q_ell(ell x)`. On `0<=x<=3`, (14) gives the
continuous uniform envelope

\[
H(x)=\ell D(\ell x),\qquad
\widehat H(x)=Q'(x)+\ell e^{-\ell x}Q(x),\qquad
|H-\widehat H|<5\,10^{-6}.
\tag{15}
\]

The exact linear trace is zero on `0<=x<=1` by finite propagation and
zero initial shell data. All support-edge jets are zero; there is no jump
at `x=1`. The comparator is zero for `2<=x<=3`, but the actual dynamic
boundary trace need not be: only the uniform small tail bound follows.

To display a **certified continuous curve**, enclose `Q'` and
`ell exp(-ell x)Q` over every whole cell, then widen by (14). Values at
cell centers or extrema inferred from a floating grid are insufficient.
The physical scaled amplitude is

\[
\ell S(\ell x)=
\frac{3\cosh^{-3/2}(\ell x)}{\rho_b^3p_b}H(x).
\tag{16}
\]

The gain is negative on this branch because `p_b<0`; interval multiplication
must preserve this sign and its uncertainty. Its `rho_b,p_b` intervals must
come from the same actual-background box. A plot of (16) is a bound for the
first parameter derivative's harmonic amplitude, not the scalar at finite
nonlinear epsilon. The full scalar variation has angular factor `Y`.
An exact nonlinear field bound requires an effective Taylor remainder,
which the present argument does not supply.

For the same `ell`, changing to `x` gives the exact norm relation
`int_0^(3ell) D(t)^2 dt = ell^-1 int_0^3 H(x)^2 dx`.
Any sharper response coefficient interval must preserve this Jacobian and
the proper-time prefactor. The waveform envelope and its model-coordinate
axis provide no measured physical time or cosmological-origin inference.
