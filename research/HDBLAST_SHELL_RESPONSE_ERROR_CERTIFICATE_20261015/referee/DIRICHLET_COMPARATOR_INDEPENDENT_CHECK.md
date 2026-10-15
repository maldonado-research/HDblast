# Independent check of the short-bump Dirichlet comparator

15 October 2026 UTC. This is an analytic conditional estimate, not a
numerical HDBLAST background certificate.

## Comparator and exact trace

For the already specified normalized bump `q=q_ell`, extend both `q` and
`v(z)=exp(-z)q(z)` oddly across zero. Their odd extensions are smooth
because the data vanish on a whole neighborhood of zero. The free
whole-line d'Alembert solution, restricted to positive depth, is a smooth
Dirichlet comparator `w_hat`. It has exactly the registered master initial
data, `w_hat(t,0)=0`, and

\[
D(t):=\widehat w_z(t,0)=q'(t)+e^{-t}q(t),
\qquad |\widehat w(t,z)|\le 1+t.
\tag{1}
\]

Put `x=(z-ell)/ell` and `u=1/[x(1-x)]>=4`. On the support,

\[
|q'(z)|=\ell^{-1}u^2|1-2x|e^{4-u}
\le\ell^{-1}u^2e^{4-u}\le16/\ell,
\tag{2}
\]

since `u^2 exp(4-u)` is decreasing for `u>=4`. Thus
`|D|<=16/ell+1`, and `D` has support in `[ell,2ell]`.

On `[0,3ell]`, integration by parts with zero endpoint data gives

\[
\int D^2=\int(q')^2+\int e^{-t}q^2+\int e^{-2t}q^2
\ge\int(q')^2\ge4/\ell.
\tag{3}
\]

The last inequality applies Cauchy-Schwarz separately to the rise from
zero to one on a segment of length `ell/2` and the fall back to zero on
another such segment. In particular

\[
2/\sqrt\ell\le\|D\|_{L^2(0,3\ell)}
\le\sqrt\ell(16/\ell+1).
\tag{4}
\]

For the error `e=w-w_hat`, the bulk source in
`e_tt-e_zz+We=F` is `F=-W w_hat`. In the convention
`-e_z+mu e_tt+beta e=G`, its boundary source is `G=D`.
The error has identically zero initial displacement and velocity.

## Independent shifted-energy and trace check

On `0<=z<=L-t`, with `L>T`, let

\[
X^2=\int(e_t^2+e_z^2+\lambda e^2)dz+
                   \mu v^2+\nu q^2,
\quad q=e(t,0),\quad v=q',\quad\lambda,\nu>0.
\tag{5}
\]

The moving right edge contributes
`-[(e_t-e_z)^2+lambda e^2]/2` to the derivative of `X^2/2`.
At the shell the bulk flux and the boundary kinetic energy cancel; their
remaining term is `(nu-beta)qv+Gv`. Consequently

\[
X'\le(a/2)X+\sqrt{f(t)^2+g(t)^2/\mu},\quad
a=\max\left\{\sup_D\frac{|\lambda-W|}{\sqrt\lambda},
                   \sup_t\frac{|\nu-\beta|}{\sqrt{\mu\nu}}\right\},
\tag{6}
\]

where `f(t)>=||F(t,.)||_L2(0,L-t)` and `g(t)>=|G(t)|`.
An integral formulation, or regularization by `sqrt(X^2+eta^2)`, justifies
this inequality when `X=0`.

For `J_+=e_t+e_z`, the incoming characteristic formula gives

\[
J_+(t,0)=\int_0^t(F-We)(s,t-s)ds.
\tag{7}
\]

Minkowski's integral inequality, with the change of variable `z=t-s`,
therefore implies, for `M>=sup|W|`,

\[
\|J_+(.,0)\|_{L^2(0,T)}
\le\int_0^T[f(s)+M X(s)/\sqrt\lambda]ds.
\tag{8}
\]

Using `e_z=J_+-v`, one obtains the physical master trace error

\[
\|e_z-(c/2)q\|_{L^2(0,T)}\le
\int_0^T[f+MX/\sqrt\lambda]dt+
\sqrt{\mu^{-1}+c^2/(4\nu)}\,\|X\|_{L^2(0,T)}.
\tag{9}
\]

This controls a derivative trace without assuming bulk `H1` energy itself
has a classical normal trace. The separate characteristic argument is
essential. All characteristics lie in the causal triangle, so there is no
artificial outer-boundary condition in this argument.

## A deliberately coarse explicit small-depth criterion

Assume fixed verified background bounds on a rectangle containing
`0<=t<=3ell`, `0<=z<=8ell`: `|W|<=M`, `|beta|<=B`,
`mu>=mu_minus>0`, and `|c|<=c_abs`. Set

\[
a_*=\max\{1+M,(1+B)/\sqrt{\mu_-}\},\qquad
k=\sqrt{\mu_-^{-1}+c_{\rm abs}^2/4}.
\tag{10}
\]

Use `lambda=nu=1`, `T=3ell`, and `L=8ell`. For `ell<=1/3`, valid
uniform forcings are `f<=6M sqrt(ell)` and `g<=17/ell`. Equations
(6)-(9) give

\[
X(t)\le e^{a_*T/2}t
       [6M\sqrt\ell+17/(\ell\sqrt{\mu_-})],
\tag{11}
\]

\[
\Delta\le18M\ell^{3/2}+e^{3a_*\ell/2}
\left[(9/2)M\ell^2+3k\ell^{3/2}\right]
\left[6M\sqrt\ell+17/(\ell\sqrt{\mu_-})\right].
\tag{12}
\]

Here `Delta` bounds the `L2` difference between the actual physical master
trace `H=w_z-cw/2` and `D`. Define

\[
R_*=9M+27M^2+\frac{153M}{2\sqrt{\mu_-}}
                     +18kM+\frac{51k}{\sqrt{\mu_-}}.
\tag{13}
\]

If

\[
0<\ell\le\min\left\{1/3,\ 1/(3a_*),\ 1/[4(1+R_*)]\right\},
\tag{14}
\]

then `exp(3a_*ell/2)<=exp(1/2)<2`, and direct expansion of (12) yields

\[
\frac{\sqrt\ell}{2}\Delta\le\ell R_*\le1/4,
\quad\hbox{so}\quad \Delta\le1/(2\sqrt\ell).
\tag{15}
\]

The power reductions used here all have `ell<=1`: for example
`ell^(3/2)<=ell` and `ell^3<=ell`. Combining with (4),

\[
\frac{3}{2\sqrt\ell}\le\|H\|_{L^2(0,3\ell)}
                 \le\frac{35}{2\sqrt\ell}.
\tag{16}
\]

Thus, if `K_minus<=3*pi^2/(rho_b^2*p_b^2)<=K_plus`,

\[
\boxed{\frac{9K_-}{4\ell}\le\mathcal C_{3\ell}
                    \le\frac{1225K_+}{4\ell}.}
\tag{17}
\]

This criterion is intentionally loose. It supplies a conditional analytic
bound for a newly fixed small-depth state. It becomes a numerical bound
only after all background/parameter/collar enclosures and the chosen
`ell` are explicitly validated. The registered exact-data and nonlinear
evolution assumptions must still apply to that state's full correction
annulus, with a strict observation-time margin. They are not supplied by
the comparator estimate. Since the initial data develop large derivatives
as `ell` shrinks, a uniform nonlinear amplitude radius does not follow.

The constants in (13)-(17) were independently expanded and checked by hand;
an exact rational audit can further guard against transcription errors.
