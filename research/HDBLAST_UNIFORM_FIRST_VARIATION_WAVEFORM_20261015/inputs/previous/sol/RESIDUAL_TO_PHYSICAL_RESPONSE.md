# A characteristic residual certificate for the physical shell response

15 October 2026 UTC. Independent Sol derivation. Conditional analytic result.

This note takes the coupled master reconstruction, straight-shell/gauge
identification, exact-family hypotheses and physical observable of
`shell-response-20261015/source/HDBLAST_PHYSICAL_SHELL_RESPONSE_20261015/PHYSICAL_SHELL_RESPONSE.md`
as inputs. It supplies an explicit residual estimate for that master's physical
normal-derivative observation. It neither replaces the master by a test scalar
nor supplies a numerical background or an effective nonlinear amplitude radius.
No other independent response branch was used to obtain the estimate below.

## 1. The observation after conjugation

Write

\[
 V=\alpha(z)h(t)w,\qquad h=\cosh^{-3/2}t,\qquad
 \alpha(0)=1,\quad \alpha_z(0)=-\gamma,\quad \gamma=C_b/2,
\]
\[
 k(t)=\frac{35}{4}\operatorname{sech}^2t-\frac{25}{4},\qquad
 \mu=(\rho_b B_b)^{-1}>0,\qquad \beta(t)=\mu k(t)+\gamma.
\]

The inherited coupled master satisfies

\[
 w_{tt}-w_{zz}+W(t,z)w=0,\qquad
 -w_z+\mu w_{tt}+\beta(t)w=0\quad(z=0),                 \tag{1}
\]

where the source notation is \(W=\mathcal W\). Define

\[
 q=w(t,0),\quad v=q',\quad j_\pm=w_t\pm w_z,\qquad
 R=w_z(t,0)-\gamma q=j_+(t,0)-v-\gamma q.               \tag{2}
\]

The physical scalar is \(S=3hR/(\rho_b^3p_b)\), including the physical
embedding/gauge identification retained from the source. Consequently

\[
 \boxed{\mathcal C_T=K\|R\|_{L^2(0,T)}^2,\qquad
 K=\frac{3\pi^2}{\rho_b^2p_b^2}.}                     \tag{3}
\]

There is no remaining hyperbolic-cosine weight: \(\cosh^3t\,h(t)^2=1\)
exactly. The boundary equation also gives
\(R=\mu(q''+kq)\), but using that identity for error control would require
an acceleration error. Formula (2) instead uses one incoming characteristic,
the boundary velocity and the boundary value.

The specified initial state becomes exactly
\(w_0=q_\ell\), \(w_1=e^{-z}q_\ell\). In particular the background factor
\(\alpha\) cancels from these initial data. An approximation formed from an
approximate \(\alpha\) must still be certified against the true coefficient
enclosures in (1); approximate conjugation alone is not a proof.

## 2. Exact error equations, including coefficient uncertainty

Fix \(T>0\) and work on the closed causal triangle

\[
 D_T=\{(t,z):0\le t\le T,\ 0\le z\le T-t\}.
\]

The physical collar must contain this triangle with the required strict margin.
One may use a larger triangle if desired. There is no artificial inner boundary
condition: both backward characteristics remain in the triangle and meet the
initial segment or the physical shell.

Let \(u\) be a smooth, explicitly reconstructible approximate master; all of its
residuals and initial derivatives below must have rigorous uniform enclosures.
Let \(\widehat W,\widehat\mu,\widehat\beta,\widehat\gamma\) be chosen approximate
coefficients. Define signed defects

\[
 r=u_{tt}-u_{zz}+\widehat W u,\qquad
 b=-u_z+\widehat\mu u_{tt}+\widehat\beta u\quad(z=0).
\]

For the convention \(e=w-u\), put

\[
 f=-r-(W-\widehat W)u,\qquad
 g=-b-(\mu-\widehat\mu)u_{tt}-(\beta-\widehat\beta)u_b.
                                                               \tag{4}
\]

Then

\[
 e_{tt}-e_{zz}+We=f,\qquad
 -e_z+\mu e_{tt}+\beta e=g.                              \tag{5}
\]

Writing \(q_e=e_b\), \(v_e=q_e'\), \(j_\pm^e=e_t\pm e_z\), the boundary
and characteristic identities have the fixed signs

\[
 (\partial_t\mp\partial_z)j_\pm^e=f-We,
 \qquad \mu v_e'+v_e+\beta q_e=j_+^e(t,0)+g,
\]
\[
 q_e'=v_e,\qquad j_-^e(t,0)=2v_e-j_+^e(t,0).             \tag{6}
\]

In particular a boundary defect enters the acceleration equation with the
**plus** sign \(+g\) under definition (4). No energy bound on \(e\) has been
substituted for a bound on its normal derivative.

Useful admissible source bounds are

\[
 F(t)\ge\sup_{0\le z\le T-t}|r(t,z)|
       +\sup_{0\le z\le T-t}|(W-\widehat W)u(t,z)|,
\]
\[
 G(t)\ge |b(t)|+|\mu-\widehat\mu|\,|u_{tt}(t,0)|
                    +|\beta-\widehat\beta|\,|u(t,0)|.   \tag{7}
\]

Every coefficient discrepancy in these inequalities is an interval bound for
the retained true background, not a discrepancy between two unvalidated
background samples. In particular uncertainty in \(\mu\) cannot simply be
dropped from the boundary residual.

### A sharper coherent boundary defect

The exact relation \(\beta=\mu k+\gamma\) should be preserved where possible.
Choose \(\widehat\mu>0\), \(\widehat\beta=\widehat\mu k+\widehat\gamma\), and define

\[
 \widehat R=u_z-\widehat\gamma u_b,\qquad
 d=u_{tt}(t,0)+k(t)u_b-\widehat R/\widehat\mu.
\]

Here \(b=\widehat\mu d\). Exact rearrangement of (4) gives

\[
 \boxed{g=-\mu d-
       \frac{\mu-\widehat\mu}{\widehat\mu}\widehat R
       -(\gamma-\widehat\gamma)u_b.}                    \tag{8}
\]

Thus, if \(0<\mu_-\le\mu\le\mu_+\),

\[
 G(t)=\mu_+|d(t)|+
     \frac{\delta_\mu}{\widehat\mu}|\widehat R(t)|
     +\delta_\gamma|u(t,0)|                            \tag{9}
\]

is valid whenever \(\delta_\mu\ge|\mu-\widehat\mu|\) and
\(\delta_\gamma\ge|\gamma-\widehat\gamma|\).
This uses the correlation between boundary coefficients, avoiding an
uncertainty factor times a separately bounded acceleration. The measured defect
\(d\) still contains \(u_{tt}\): claiming it is small requires an acceleration
residual enclosure or a rigorously consistent first-order reconstruction.

## 3. An explicit global finite-time bound

Assume continuous bounds on the compact triangle,

\[
 |W|\le M,\quad |\beta|\le B,\quad
 0<\mu_-\le\mu,\quad a=\mu_-^{-1},\quad F,G\ge0.
\]

For the source coefficients one can always use
\(B=(25/4)\mu_++\bar\gamma\), since \(|k|\le25/4\), and
\(M=25/4+\|C_z/2+C^2/4+2\rho^2p^2/3\|_\infty\).
Sharper interval bounds on the actual compact triangle are preferable.
The proof below needs no time derivatives of \(W,\beta,F,G\); its stronger
initial derivative norms and continuous residual bounds are stated explicitly.

Let the initial error norms on \([0,T]\) be

\[
 A_0=\|e(0,\cdot)\|_\infty,\quad
 b_+=\|e_t(0,\cdot)+e_z(0,\cdot)\|_\infty,\quad
 b_-=\|e_t(0,\cdot)-e_z(0,\cdot)\|_\infty,
\]
\[
 J_0=\max(b_+,b_-),\qquad V_0=|e_t(0,0)|.              \tag{10}
\]

For the registered data, these are errors relative to \(q_\ell\) and
\(e^{-z}q_\ell\), including the derivative of the displacement. Initial
compatibility for the true smooth solution remains a retained hypothesis.

Solve the following **three-dimensional nonnegative majorant system**:

\[
 \begin{aligned}
 E'&=J_0+I+V,&E(0)&=A_0,\\
 I'&=ME+F,&I(0)&=0,\\
 V'&=aBE+aI+a(b_++G),&V(0)&=V_0.
 \end{aligned}                                             \tag{11}
\]

All variables on the right of (11) are bounds, not fields of the physical
model. Its explicit representation is

\[
 y(t)=e^{At}y_0+\int_0^t e^{A(t-s)}c(s)\,ds,
\quad A=\begin{pmatrix}0&1&1\\M&0&0\\aB&a&0\end{pmatrix},
\quad c=(J_0,F,a(b_++G))^T.                              \tag{12}
\]

The matrix is entrywise nonnegative, so validated upper bounds from a positive
series, a rational Taylor remainder bound, or interval ODE integration preserve
the inequality. No spectral-stability assumption or positive potential is
needed. The bound can grow exponentially on a finite interval.

**Claim.** At every point of \(D_T\),

\[
 |e(t,z)|\le E(t),\quad |v_e(t)|\le V(t),\quad
 |j_+^e(t,z)|\le b_++I(t).                              \tag{13}
\]

**Proof.** Let \(E_*(t)\) be the supremum of \(|e|\) up to time \(t\) inside
the triangle, \(V_*(t)=\sup_{s\le t}|v_e(s)|\), and
\(I_*(t)=\int_0^t(ME_*(s)+F(s))ds\). Backward integration of \(j_+^e\)
always reaches the initial line, giving \(|j_+^e|\le b_++I_*\). An outgoing
\(j_-^e\) either starts on the initial line or starts at the shell. In the
second case (6) and the incoming characteristic give

\[
 |j_-^e(t,z)|\le b_++2V_*(t)+I_*(t).
\]

The source integrals before and after reflection are disjoint in time, so there
is only one \(I_*\) in this inequality. Both cases imply

\[
 |e_t(t,z)|\le J_0+V_*(t)+I_*(t),\quad
 E_*(t)\le A_0+\int_0^t(J_0+V_*+I_*)ds.                \tag{14}
\]

The shell ODE has the exact damped variation-of-constants formula

\[
 v_e(t)=e^{-t/\mu}v_e(0)
   +\mu^{-1}\int_0^t e^{-(t-s)/\mu}
       [j_+^e(s,0)+g(s)-\beta(s)q_e(s)]ds.
\]

Using \(|q_e|\le E_*\), the positive damping kernel and \(\mu^{-1}\le a\)
yields

\[
 V_*(t)\le V_0+a\int_0^t[b_++I_*(s)+G(s)+BE_*(s)]ds.   \tag{15}
\]

Equations (14), (15), and the definition of \(I_*\) form a positive linear
Volterra inequality. Iteration or the nonnegative resolvent series bounds them
by (11). This proves (13). \(\square\)

A coarser scalar version is useful when no matrix computation is desired.
For constant upper bounds \(F\le F_*\), \(G\le G_*\), put

\[
 L=\max\{2,M,a(B+1)\},\quad
 Y_0=\max\{A_0,V_0\},\quad
 H=\max\{J_0,F_*,a(b_++G_*)\}.
\]

Then all three components of (11) are bounded by

\[
 Y(t)=Y_0e^{Lt}+H\frac{e^{Lt}-1}{L}.                   \tag{16}
\]

This explicit bound is conservative but finite for every fixed \(T\).
Sharper boundary damping or sign information can improve it; neither is needed
for validity.

## 4. A short-slab characteristic contraction with explicit constants

For completeness, take a slab of length \(h\) and two characteristic iterates
with the same initial data, sources and coefficients, whose input masters differ
by at most \(D\). Integrate the two \(j_+\) equations, solve the boundary ODE
with the same starting \((q,v)\), reflect using (6), and reconstruct the new
master by integrating \((j_++j_-)/2\). The changes satisfy

\[
 |\delta j_+|\le MhD,\qquad
 \sup|\delta v|\le \frac{aMh^2}{2}D
                  +\frac{aBh^2}{2}\sup|\delta v|.
\]

If \(\theta=aBh^2/2<1\), then

\[
 \|\delta w_{\rm new}\|_\infty\le\kappa(h)D,
 \qquad
 \kappa(h)=\frac{Mh^2}{2}
       +\frac{aMh^3}{2(1-aBh^2/2)}.                    \tag{17}
\]

For example the fully explicit positive choice

\[
 h\le\min\left\{1,\ [a(B+1)]^{-1/2},\
                  [2(M+1)]^{-1/2},\ [4a(M+1)]^{-1/3}\right\}
                                                               \tag{18}
\]

gives \(\theta<1/2\) and \(\kappa<1/2\).
If \(M=0\) the input master does not enter the characteristic source and the
map has contraction constant zero. At subsequent slabs the same argument uses
the inherited compatible initial characteristic and shell values; finitely
many steps cover every fixed interval. Derivative identities propagate from
consistent initial derivatives as in the source construction. This argument
does not create new background or metric-evolution hypotheses.

## 5. Observation error and a certified coefficient interval

Let \(\bar\gamma\ge|\gamma|\),
\(\delta_\gamma\ge|\gamma-\widehat\gamma|\). The exact observation identity is

\[
 R-\widehat R=j_+^e(t,0)-v_e(t)-\gamma q_e(t)
                   -(\gamma-\widehat\gamma)u_b(t).
\]

Consequently the computable envelope

\[
 D_R(t)=b_++I(t)+V(t)+\bar\gamma E(t)
                   +\delta_\gamma|u_b(t)|             \tag{19}
\]

satisfies \(|R-\widehat R|\le D_R\). Any outward bound

\[
 \varepsilon_R\ge\left(\int_0^T D_R(t)^2dt\right)^{1/2}           \tag{20}
\]

is an observation-norm error certificate. A simple fallback is
\(\varepsilon_R=\sqrt T\,[b_++I(T)+V(T)+\bar\gamma E(T)+
\delta_\gamma\|u_b\|_\infty]\), since the majorants are nondecreasing.
For exact \(\gamma\), the last term is zero. A weighted scalar-trace bound is
also immediate:

\[
 \|\cosh^{3/2}t(S-\widehat S)\|_2
       \le \frac{3}{|\rho_b^3p_b|}\varepsilon_R,
 \quad \widehat S=\frac{3h}{\rho_b^3p_b}\widehat R,       \tag{21}
\]

when the physical prefactor is exact. If it is interval-valued, certify the
coefficient using (3) directly, rather than treating an uncertain prefactor as
exact in (21).

Suppose rigorous quadrature gives
\(\|\widehat R\|_2\in[N_-,N_+]\), and validated background intervals give
\(K\in[K_-,K_+]\) with \(K_->0\). Then

\[
 \boxed{\mathcal C_T\in
 \left[K_-\max\{0,N_--\varepsilon_R\}^{\,2},
       K_+(N_++\varepsilon_R)^2\right].}                \tag{22}
\]

For \(\rho_b\in[\rho_-,\rho_+]\) and
\(|p_b|\in[p_-,p_+]\), all lower endpoints positive,
\(K_-=3\pi_-^2/(\rho_+^2p_+^2)\) and
\(K_+=3\pi_+^2/(\rho_-^2p_-^2)\) are valid outward bounds, with a rigorous
interval for \(\pi\). If exact rational \(L\le\|\widehat R\|_2^2\le U\) is
available, outward square-root bounds for \(L,U\) suffice; exact-square rational
bounds avoid transcendental square-root evaluation entirely.

The certificate **decides positive** if \(N_->\varepsilon_R\). A zero lower
endpoint means the numerical certificate is inconclusive; it is not evidence
that the physical response vanishes. The source's qualitative activation
theorem remains conditional evidence of positivity. In particular its
contradiction-only auxiliary whole-line lower bound must not be inserted as a
lower bound for the actual dynamic-boundary trace. Near the flat arrival front
the true response can be extremely small, so refinement or a later observation
window may be necessary.

One may replace \([0,T]\) in the norm lower bound by any certified collection of
observation windows, using the same local error envelope, because their response
norm is a lower bound for the total norm. Such a lower bound remains tied to the
specified coupled initial data and actual background coefficients.

## 6. A derivative-loss control and what has not been certified

It is invalid to turn a small \(\|w-u\|_\infty\) alone into a small observation.
As an exact diagnostic for the abstract dynamic-boundary operator with
\(W=\beta=\gamma=0\), let

\[
 u_\omega(t,z)=\left[\frac{\cos(\omega z)}{\mu\omega^2}
                    -\frac{\sin(\omega z)}{\omega}\right]
                    \cos(\omega t),\qquad \mu>0.
\]

It satisfies both zero bulk residual and \(-u_z+\mu u_{tt}=0\). Its supremum
on any fixed finite collar is at most \(1/(\mu\omega^2)+1/\omega\to0\), while
\(u_z(t,0)=-\cos(\omega t)\). The observation does not tend to zero in
\(L^2(0,T)\): its squared norm is
\(T/2+\sin(2\omega T)/(4\omega)\to T/2\). The initial characteristic
derivative errors in (10) detect this example. It is a negative control for an
estimate, not a claimed physical HDBLAST background or a substitute state.

Equations (11), (19), and (22) are a rigorous conditional analytic route to a
validated number. They do not themselves provide \(M,B,\mu_\pm,\gamma\), a
physical waveform, certified residual enclosures, a physical quadrature result,
or an effective nonlinear radius. Smooth sampled residuals are insufficient:
the supremum bounds, coefficient enclosures, initial errors and quadrature must
be certified on their entire cells. The second exact-family correction does
not enter \(\mathcal C_T\), but this linear certificate does not bound the
nonlinear Taylor remainder.

The accompanying `algebra/` audit checks exact identities and rejects specified
wrong controls. Finite symbolic checks support bookkeeping; the characteristic
comparison proof above supplies the analytic estimate.

`rational_majorant.py` evaluates (12) and the constant-envelope fallback in
(20) using `fractions.Fraction` only. It augments the inhomogeneous ODE by a
constant component and encloses its positive exponential series with a geometric
tail once the successive norm ratio is below one. Its default self-checks and
negative controls are saved in `rational_majorant_checks.json`. The examples
exercise the operator estimate and are explicitly not physical background data.
The script's lower and upper endpoints enclose the majorant ODE solution; its
lower endpoints are not lower bounds on the unknown actual error.
