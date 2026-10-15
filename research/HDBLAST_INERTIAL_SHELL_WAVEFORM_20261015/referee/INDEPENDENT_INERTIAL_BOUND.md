# Independent inertial comparator and residual estimate

This note concerns a conditional linear first variation on the retained
background. It supplies no finite-amplitude nonlinear solution or remainder.
All source files from earlier checkpoints remain unchanged.

Let `l=ell`, `T=3l`, and assume the actual causal triangle lies in the certified
collar. Set `a0=35/4` and write

\[
W(t,z)=V(z)-a_0\tanh^2t,\qquad
\beta(t)=\beta_0-a_0\mu\tanh^2t,\qquad
\beta_0=\gamma+5\mu/2.
\]

The actual field has initial data `f=q_l`, `g=e^{-z}q_l`, where `q_l` is smooth,
supported in `[l,2l]`, and `0<=q_l<=1`. Assume throughout the required causal
triangle

\[
\mu\ge m>0,\quad |W|\le M,\quad |V'|\le K,
\quad |\gamma|\le\Gamma,\quad |\beta_0|\le B_0,
\quad |w|\le A,\quad |\partial_tw(t,0)|\le B_v.
\tag{1}
\]

Here the previously established positive-mass majorants can provide `A,B_v`.
No source derivative or estimate of a second time derivative is needed below.

## 1. Free bulk comparator with dynamic boundary

Let `u` solve the free bulk equation, retain the same initial data and satisfy
the constant dynamic boundary `-u_z+mu u_tt+beta_0 u=0`. Its boundary value
`b(t)=u(t,0)` satisfies exactly

\[
\mu b''+b'+\beta_0b=j_0(t),\qquad
j_0=q_l'+e^{-t}q_l,\quad b(0)=b'(0)=0.
\tag{2}
\]

The free physical trace is therefore

\[
D_u=j_0-b'-\gamma b.
\tag{3}
\]

Positive inertia must be retained. In particular, replacing `(2)` by a
Dirichlet boundary discards the response being modeled.

## 2. Incoming Dirichlet field and exact cancellation

Let `u_D` be the free Dirichlet solution with the same data. For odd extensions
`f_o,g_o`, along the incoming characteristic ending at `(t,0)` it is

\[
u_D(s,t-s)=\frac{f(t)+f_o(t-2s)}2
 +\frac12\int_{t-2s}^{t}g_o(r)\,dr.
\tag{4}
\]

Consequently the displacement-only and velocity-only integrals obey

\[
\int_0^t u_D^{(f)}(s,t-s)\,ds=\frac{t f(t)}2,
\qquad
\int_0^t u_D^{(g)}(s,t-s)\,ds=
\frac12\int_0^t r e^{-r}q_l(r)\,dr.
\tag{5}
\]

The first is exactly zero for `t>=2l`. The second is positive and generally
nonzero after the pulse. This is not cancellation of the complete field.
Moreover `u_D^(g)(s,t-s)>=0`, so absolute values inside its integral do not
lose the moment estimate.

For the potential-weighted displacement, the odd part admits the exact
paired formula

\[
\begin{split}
\int_0^t W(s,t-s)u_D^{(f)}(s,t-s)\,ds
={}&\frac{q_l(t)}2\int_0^t W(s,t-s)\,ds\\
&+\frac14\int_0^t\left[
W\left(\frac{t-r}2,\frac{t+r}2\right)
-W\left(\frac{t+r}2,\frac{t-r}2\right)\right]q_l(r)\,dr.
\end{split}
\tag{6}
\]

The bracket is bounded by `(K+a0 t)r`: use the radial Lipschitz bound and
`|tanh^2(b)-tanh^2(a)|<=b^2-a^2` for `0<=a<=b`. Since
`integral r q_l(r) dr <=3l^2/2`, equations (5)-(6) imply

\[
\left|\int_0^t W u_D\,ds\right|
\le \frac{Mtq_l(t)}2+R_0,
\quad R_0=\frac38(K+a_0T)l^2+\frac34Ml^2.
\tag{7}
\]

Only the first term can be order `Ml`. It is absent after `2l`.

## 3. Boundary lift and potential scattering cannot be omitted

The actual field has the exact free-Dirichlet representation

\[
w(t,z)=u_D(t,z)+w_b(t-z)\mathbf1_{t\ge z}
 +\mathcal R[-Ww](t,z),\qquad w_b(t)=w(t,0).
\tag{8}
\]

The half-line Dirichlet retarded operator, obtained by odd reflection, obeys
`|R[-Ww](s,z)|<=MA s^2/2`. This holds with the one-sided source regularity
inherited from the smooth actual solution; a globally smooth odd source is
not required. All past integration domains remain inside the causal triangle.

Set

\[
I(t)=\int_0^t W(s,t-s)w(s,t-s)\,ds.
\]

The boundary lift contributes at most `M B_v t^2/4`, since on this
characteristic its integral is `1/2 integral_0^t W((t+r)/2,(t-r)/2)w_b(r)dr`
and `|w_b(r)|<=B_v r`. The scattered field contributes at most
`M^2 A t^3/6`. Hence

\[
|I(t)|\le \frac{Mtq_l(t)}2+R_0+
 \frac{M B_vt^2}{4}+\frac{M^2At^3}{6}.
\tag{9}
\]

On `[2l,T]` a uniform bound is

\[
I_{\rm tail}=R_0+\frac{M B_vT^2}{4}+\frac{M^2AT^3}{6}.
\tag{10}
\]

The whole history is bounded more sharply than `T sup|I|`:

\[
\int_0^T|I(t)|\,dt\le J_1
=\frac34Ml^2+TR_0+\frac{M B_vT^3}{12}
 +\frac{M^2AT^4}{24}.
\tag{11}
\]

Here the first term follows by integrating `Mtq_l(t)/2` over its actual
support. Its memory remains relevant after the forcing stops.

## 4. Actual-minus-comparator boundary error

Set `e=w-u`, `q=e(t,0)`, `v=q'`. Incoming transport gives
`e_t+e_z=-I` at the shell. Subtracting the two exact boundary equations gives

\[
\mu v'+v+\beta_0q=-I(t)+a_0\mu\tanh^2t\,w_b(t).
\tag{12}
\]

This uses the actual boundary value on the right; no coefficient discrepancy
has been lost. With zero error initial traces and the positive exponential
kernel, if `s0=B_0T^2/(2m)<1` then

\[
\|v\|_{\infty,[0,T]}\le E_v
=\frac{J_1/m+a_0B_vT^4/4}{1-s_0},
\qquad \|q\|_\infty\le TE_v.
\tag{13}
\]

The time-coefficient forcing uses
`integral tanh^2(s)|w_b(s)| ds <= B_v T^4/4`.
The physical normalized trace error satisfies

\[
|H(x)-H_u(x)|\le l\{|I(lx)|+(1+\Gamma T)E_v\}.
\tag{14}
\]

Thus a global bound uses (9), and the postpulse bound uses (10). An estimate
using only postpulse `I` inside the velocity convolution is invalid: it
would omit the accumulated on-pulse error.

## 5. Leading free postpulse coefficient

For fixed dimensionless `x>=2`, put `I_Q=integral_1^2 Q(r)dr>0`. Expansion of
the exact constant-coefficient comparator at small `l` gives

\[
b(lx)=\frac{l}{\mu}I_Q+O(l^2),\qquad
b'(lx)=\frac{l}{\mu}\left(1-\frac1\mu\right)I_Q+O(l^2),
\]

and therefore

\[
H_u(x)=\frac{l^2 I_Q}{\mu}
\left(\frac1\mu-\gamma-1\right)+O(l^3).
\tag{15}
\]

The `-1` comes from the prescribed nonzero initial velocity. With zero
initial velocity it would be absent. Positivity of `1/mu-gamma` alone is
not the exact leading-coefficient condition for the present data. A finite
width positivity claim requires an explicit remainder or full interval
enclosure; a formal leading term by itself is insufficient.

## 6. Independent derivative identity

For `h=U_eta/p`, the retained background ODE implies
`h_z=rho(h^2-4ah-U_etaeta)`. Total differentiation of the derivative-free
potential gives

\[
\begin{split}
V_z={}&\rho^3\left[4h^3-30ah^2+60a^2h-\frac{35}{2}a^3
+(12a-4h)U_{\eta\eta}
+\left(\frac{29}{6}a-3h\right)p^2+pU_{\eta\eta\eta}\right]\\
&+\rho\left(\frac{35}{2}a-10h\right).
\end{split}
\tag{16}
\]

This identity is a cross-check for an outward interval `K` calculation.
Neither an unsampled floating derivative nor a finite difference is a
replacement for that calculation.

## 7. Concrete outward conclusions for the new wider state

The authenticated derivative calculation and independent coefficient audit supply
`K<30000000`, while the retained wider rectangle gives
`m=.002, M=102000, B_0=250, Gamma=250`. For `l=10^-6`, the exact rational
checker accompanying this note proves, with all on-pulse error memory kept,

\[
\sup_{0\le x\le3}|H-H_u|<1.022\times10^{-7},\qquad
\sup_{2\le x\le3}|H-H_u|<1.645\times10^{-10}.
\tag{17}
\]

For an additional independent coarse tail-sign proof, positive energy of the
constant-coefficient Green kernel `k` yields `|k'|<=1/mu`, `|k|<=t/mu`.
The integration-by-parts kernel is

\[
P(r,s)=(1/\mu-\gamma-e^{-s})k'(r)
 +(\beta_0/\mu-\gamma e^{-s})k(r).
\]

Its value at `r=0` is at least `59600`, using `mu<1/400` and `gamma<250`.
The bound

\[
|P_r|\le K_3+(\Gamma+1)K_2+\Gamma/m,
\quad K_2=(1+B_0T)/m^2,\quad K_3=(K_2+B_0/m)/m
\]

then bounds `P` below on the entire needed rectangle. Also
`integral Q>1/6`: on `[4/3,5/3]`, `Q>=exp(-1/2)>1/2`.
Combining these inequalities with (17) proves continuously

\[
\boxed{H(x)>9.64\times10^{-9}\quad(2\le x\le3).}
\tag{18}
\]

This is an independently derived conservative sign bound. The physical
scalar first variation has negative gain on the retained shell, so its
harmonic amplitude has the opposite sign. The coefficient of the angular
harmonic is being discussed; the full angular field also contains the
sign-changing quadrupole harmonic.

The companion checker authenticates its two data inputs before reading them,
uses exact `Fraction` decisions, and passes 21 groups in normal and optimized
execution modes. Its floating display fields are not decision endpoints.
