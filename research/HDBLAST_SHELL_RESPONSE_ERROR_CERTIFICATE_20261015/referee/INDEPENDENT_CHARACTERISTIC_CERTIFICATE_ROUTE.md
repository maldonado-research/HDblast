# Independent referee route: residuals to the physical shell coefficient

15 October 2026 UTC. Independent analytic derivation, not a numerical
certificate. This note concerns a conditional bound for the already retained
smooth linear shell problem. It does not certify the retained nonlinear
existence theorem, provide an actual-background coefficient enclosure, or
produce a numerical amplitude.

## 1. Retained exact problem and observable

Read against the published `PHYSICAL_SHELL_RESPONSE.md`, equations (9),
(11a), (13), and (20). In the inward coordinate the conjugated master obeys

\[
w_{tt}-w_{zz}+W(t,z)w=0,
\qquad -w_z(t,0)+\mu w_{tt}(t,0)+\beta(t)w(t,0)=0,
\qquad \mu=(\rho_b B_b)^{-1}>0.
\tag{1}
\]

The specified initial data are `w(0,z)=q_ell(z)` and
`w_t(0,z)=exp(-z) q_ell(z)`. Let `c=C(0)` and
`A(t)=cosh(t)^(-3/2)`. Since `alpha(0)=1` and
`alpha_z(0)=-c/2`, the physical trace is exactly

\[
S(t)=\frac{3}{\rho_b^3p_b}A(t)H(t),\qquad
H(t)=w_z(t,0)-\frac c2 w(t,0).
\tag{2}
\]

Consequently the proper-time response coefficient simplifies to

\[
\boxed{\mathcal C_T=\frac{3\pi^2}{\rho_b^2p_b^2}
                  \int_0^T H(t)^2\,dt.}
\tag{3}
\]

The `cosh^3(t)` clock-volume weight cancels the squared conjugating
factor; `rho_b^4` still enters before combining it with the physical trace.
Omitting that proper-time factor would incorrectly replace `rho_b^2` in
the denominator by a different power. The normal embedding contribution
has already been accounted for by the published straight-shell
representative and its gauge bridge.

## 2. Domain, approximation, and true residuals

Fix `Z>T>0` and use only the characteristic triangle

\[
D=\{(t,z):0\le t\le T,\ 0\le z\le Z-t\}.
\tag{4}
\]

All backward characteristics used below stay in `D` and reach either its
initial edge or the physical shell. No condition on an artificial boundary
at `z=Z` is needed. The actual positive-cone collar and retained existence
domain must contain `D`; this is an input, not a consequence of a residual
calculation. A smaller causal triangle can be used for a shorter shell
observation interval.

Let `w_hat` be a specified `C^2` approximation on `D`, with the regularity
needed for the displayed traces. Define residuals with the **true**
coefficients:

\[
R=\widehat w_{tt}-\widehat w_{zz}+W\widehat w,
\qquad r=-\widehat w_z+\mu\widehat w_{tt}+\beta\widehat w
\quad(z=0).
\tag{5}
\]

For `e=w-w_hat`, put `e0=e(0,.)`, `e1=e_t(0,.)`, and choose finite
nonnegative enclosures

\[
\begin{split}
a_0&\ge\|e_0\|_\infty,\\
a_1&\ge\max\{\|e_1+e_0'\|_\infty,
                    \|e_1-e_0'\|_\infty\},\\
M(t)&\ge\sup_{0\le z\le Z-t}|W(t,z)|,\quad
f(t)\ge\sup_{0\le z\le Z-t}|R(t,z)|,\\
B(t)&\ge|\beta(t)|,\quad b(t)\ge|r(t)|.
\end{split}\tag{6}
\]

Supremum, residual, coefficient, and initial derivative bounds must be
proved over their continuous domains. A finite mesh sample is insufficient.
An interval polynomial or spline representation is one possible route; a
piecewise approximation needs matching regularity or explicit interface
residual terms.

For uncertain but rigorously enclosed background parameters suppose
`0<mu_minus <= mu <= mu_plus`. The fixed-coefficient case sets both
bounds equal to `mu`.

## 3. Characteristic bound that avoids acceleration loss

Set `J_+=e_t+e_z`, `J_-=e_t-e_z`, `q=e(t,0)`, and `v=q'`. The error
problem gives

\[
(\partial_t-\partial_z)J_+=-We-R,
\quad(\partial_t+\partial_z)J_-=-We-R,
\quad \mu v'+v+\beta q=J_+(t,0)-r,
\quad J_-(t,0)=2v-J_+(t,0).
\tag{7}
\]

In particular the trace identity is `e_z(t,0)=J_+(t,0)-v(t)`.
Boundary acceleration need not be estimated. Let

\[
E_*(t)=\sup_{0\le z\le Z-t}|e(t,z)|,
\quad I_*(t)=\int_0^t[M(s)E_*(s)+f(s)]\,ds,
\quad Q_*(t)=|q(t)|+\mu|v(t)|.
\tag{8}
\]

The incoming characteristic formula yields `|J_+(t,z)| <= a1+I_*(t)`.
The upper right Dini derivative obeys

\[
D^+Q_*\le BQ_*+a_1+I_*+b.
\tag{9}
\]

Indeed, the `+|v|` from differentiating `|q|` cancels the damping
`-|v|` in `mu D^+|v|`. This cancellation also holds at zeros by the
usual upper Dini derivative inequality.

Define a nonnegative cooperative majorant by

\[
\boxed{\begin{aligned}
E'&=a_1+I+Q/\mu_-,\\
I'&=M E+f,\\
Q'&=B Q+a_1+I+b,
\end{aligned}\qquad
E(0)=a_0,\ I(0)=0,\
Q(0)=|e_0(0)|+\mu_+|e_1(0)|.}
\tag{10}
\]

All three majorants are nondecreasing. To verify the first row, an
outgoing characteristic either starts at the initial edge or at a shell
time `s<=t`. Its magnitude is at most
`a1+I_*(t)+2 sup_{s<=t}|v(s)|`. Integrate
`e_t=(J_++J_-)/2` at fixed depth; fixed-depth vertical segments also
stay in `D`. Since `|v(s)|<=Q(s)/mu_minus<=Q(t)/mu_minus`, this gives
the first row as an integral supersolution. The other rows follow from
(8)-(9). Standard comparison for positive linear Volterra systems gives

\[
E_*\le E,\quad I_*\le I,\quad Q_*\le Q,
\quad |e_z(t,0)|\le a_1+I(t)+Q(t)/\mu_-.
\tag{11}
\]

This can also be proved by Picard iteration of the nonnegative integral
inequalities, avoiding any assumption that a supremum is differentiable.
The positive majorant need not be optimal; it must only be a certified
supersolution of (10).

## 4. An explicit fallback with only scalar exponential bounds

For uniform bounds `M<=M0`, `B<=B0`, `f<=f0`, `b<=b0`, set

\[
L=\max\{M_0,2,B_0+1/\mu_-\},\quad
g=2a_1+f_0+b_0,\quad
Y_0=a_0+|e_0(0)|+\mu_+|e_1(0)|.
\tag{12}
\]

As `Y=E+I+Q` obeys `Y'<=L Y+g` and `L>=2`,

\[
E(t),I(t),Q(t)\le
\overline Y(t)=e^{Lt}Y_0+\frac gL(e^{Lt}-1).
\tag{13}
\]

This is deliberately coarse but supplies a closed, computable alternative
to a validated ODE solve. Exponential upper bounds still require directed
rounding or a proved rational Taylor remainder. Ordinary floating-point
evaluation does not become certified by inserting it into this formula.

## 5. Coefficient error and a certified coefficient interval

Choose a real surrogate `c_hat`, define
`H_hat=w_hat_z(t,0)-c_hat*w_hat(t,0)/2`, and enclose
`|c|<=c_abs`, `|c-c_hat|<=delta_c`. Equations (2) and (11) yield

\[
|H-\widehat H|\le d(t):=
a_1+I(t)+\left(\mu_-^{-1}+c_{\rm abs}/2\right)Q(t)
                +\frac{\delta_c}{2}|\widehat w(t,0)|.
\tag{14}
\]

An independently certified bound `D >= (int_0^T d(t)^2 dt)^(1/2)`
and certified `N_minus <= ||H_hat||_L2 <= N_plus` give

\[
\max\{0,N_- -D\}\le\|H\|_{L^2}\le N_++D.
\tag{15}
\]

If `K_minus <= 3*pi^2/(rho_b^2*p_b^2) <= K_plus`, with
`0<K_minus<=K_plus`, then

\[
\boxed{K_-\max\{0,N_--D\}^2\le\mathcal C_T
                 \le K_+(N_++D)^2.}
\tag{16}
\]

A **positive numerical lower bound** requires `N_minus>D`. The prior
qualitative nonvanishing theorem does not supply this numerical inequality.
Intervals for `rho_b` must be strictly positive and those for `p_b` must
exclude zero. Use the minimum/maximum absolute values of the latter interval
when bounding `K`; squaring an interval that crosses zero is invalid here.
Validated bounds for `pi` are also needed in a rational certificate.

When the approximation uses coefficients `W_hat,mu_hat,beta_hat`, compute
residuals with them only if the following additional terms are enclosed:

\[
|R|\le|R_{\rm hat}|+|W-W_{\rm hat}|\,|\widehat w|,
\tag{17}
\]
\[
|r|\le|r_{\rm hat}|+|\mu-\mu_{\rm hat}|\,|\widehat w_{tt}|
                       +|\beta-\beta_{\rm hat}|\,|\widehat w|.
\tag{18}
\]

These are true coefficient errors, not optional discretization corrections.
The approximation's boundary second derivative occurs here as an explicitly
known function; it does not demand a bound on the exact boundary acceleration.

## 6. Referee blockers and scope

1. A bulk `H1 x L2` energy estimate without an additional normal trace or
   boundary acceleration argument cannot directly certify this observable.
2. A small residual measured only at collocation nodes is insufficient.
3. Boundary residuals, initial derivative errors, and coefficient errors
   cannot be silently set to zero.
4. Finite-domain reflections require either explicit residual accounting or
   a proven domain of dependence that excludes them. Triangle (4) supplies
   the latter without assuming an exact artificial boundary condition.
5. A smooth background solution with rigorously enclosed `W,beta,mu,c,rho_b,p_b`
   on the actual collar is required for a numerical HDBLAST interval.
6. Toy coefficients, sampled real coefficients, high precision without
   rounding control, and ordinary convergence plots are not certificates.
7. The strict interior observation margin and physical proper-time clock
   remain those of the retained exact-family construction.
8. A bound for `C_T` alone supplies no explicit nonzero amplitude radius for
   the nonlinear family; that requires a bound on its Taylor remainder.

The theorem above is a residual-to-observable framework. It can be stated
conditionally before any numerical data are available, with an explicit
list of missing numerical inputs. It should not be advertised as a
computed response magnitude or as an unconditional exact-solution theorem.
