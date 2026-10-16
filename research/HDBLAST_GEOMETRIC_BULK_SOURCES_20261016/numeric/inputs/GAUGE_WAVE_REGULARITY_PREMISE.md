# Constructive mixed-boundary gauge wave: continuous solution and C3 bootstrap

This is an independent proof for the actual three-amplitude first gauge wave
on the controlled causal triangle. It establishes a continuous C1 solution
by contraction and, with the retained compatible C2 coefficients and source,
a classical C3 solution. It is not a high-Sobolev theorem, a nonlinear
Einstein initial-boundary theorem, a second evolution, or a computation of
the transformed quadratic forcing.

The raw amplitudes are `U=(f,b,d)=(xi^y,xi^t,d)`, with angular factors
`Y,Y,grad Y` respectively. No normalization of these coordinate amplitudes
is implicit. Put `x=t/ell`, `r=z/ell`, `ell=10^-6`, and

\[
 D=\{0\le x\le X=3,\quad r\ge0,\quad x+r\le8\}.
\]

## 1. The physical half-line system and an actual C1 bound

Write `H=tanh(ell x)`, `J=sech^2(ell x)`, `lambda=8` and use background
coefficients at `z=ell r`. The published provisional operator becomes

\[
 U_{xx}-U_{rr}=P U_x+Q U_r+R U+s,\qquad
 U(0,r)=U_x(0,r)=0,
\]

with `f(x,0)=0`, `b_r(x,0)=d_r(x,0)=0`, and

\[
 P=\ell\begin{pmatrix}-3H&-2\rho^2a&0\\-2a&-3H&0\\0&0&-5H\end{pmatrix},
 \quad Q=-\ell\rho a\operatorname{diag}(3,5,5),
\]

\[
 R=\ell^2\begin{pmatrix}
 -\lambda J-4\rho^2a_y-8\rho^2a^2&-6\rho^2aH&2\lambda\rho^2a\\
 0&3+3H^2-\lambda J&-2\lambda H\\
 2aJ&2HJ&(4-\lambda)J
 \end{pmatrix},
 \qquad s=(-6\ell^2\rho^2(u_y+4au),0,0).
\]

Here `u=C cosh^(-3/2)(t) w`, where C is the retained reconstruction
coefficient, not a gauge-vector component. The signs follow from
`D_y=-rho^-1 D_z` and `D_y^2=rho^-2 D_z^2+(a/rho)D_z`.
This yields the `-3rho a` radial drift for f and `-5rho a` for b,d.

For a matrix row, sum absolute values of all entries. Set

\[
 L_1=\sup_D\max_i\sum_j(|P_{ij}|+|Q_{ij}|),\qquad
 L_0=\sup_D\max_i\sum_j|R_{ij}|.
\]

The exact companion checker authenticates the retained whole-collar actual
background and global first-master certificate. Their conservative rational
boxes prove

\[
 L_1<0.003889,\quad L_0<3.3\,10^{-8},\quad
 q:=X(L_1+XL_0/2)<0.012,
 \qquad \|s\|_\infty<8\,10^{-7}.
\]

Consequently the construction below gives the strict continuous bounds

\[
 \boxed{\max_i\|U_i\|_\infty<3.644\,10^{-6},\qquad
 \max_i\{\|U_{i,x}\|_\infty,\|U_{i,r}\|_\infty\}
       <2.430\,10^{-6}.}
\]

Physical t,z derivatives divide the scaled derivative bound by ell.
These are coordinate-amplitude bounds. Conversion to covariant vector or
metric/scalar component norms requires its explicit background factors.
The checker stores sharper exact rational majorants, not only these caps.

## 2. Retarded maps preserve the boundary and are C1

For continuous physical forcing g on D define its odd/even spatial extension
only as integration notation, `g_D(s,q)=sgn(q)g(s,|q|)` and
`g_N(s,q)=g(s,|q|)`. Then

\[
 K_{D/N}g(x,r)=\frac12\int_0^x\int_{r-x+s}^{r+x-s}g_{D/N}(s,q)dq\,ds.
\]

Equivalently use the physical positive-depth kernel with indicator difference
for D and indicator sum for N. Both kernels are nonnegative. Their value
row masses are bounded by `x^2/2` (equal for N). Their first derivatives are

\[
 K_x=\tfrac12\int_0^x[g_{D/N}(s,r+x-s)+g_{D/N}(s,r-x+s)]ds,
\]
\[
 K_r=\tfrac12\int_0^x[g_{D/N}(s,r+x-s)-g_{D/N}(s,r-x+s)]ds.
\]

Each first-derivative absolute row mass is at most x, for **both** modes.
The value of the odd extension at zero is irrelevant. Its possible jump
does not obstruct C1: each endpoint ray crosses the jump at at most one
integration time. Difference quotients of the inner primitive are uniformly
bounded by the forcing supremum, so dominated convergence gives the formulas.
The same argument gives continuity of the ray integrals, including the
characteristic `x=r` and one-sided shell/domain-edge traces.

At the shell,

\[
 K_Dg=0,\quad(K_Dg)_x=0,\quad(K_Dg)_r=\int_0^xg(s,x-s)ds,
\]
\[
 (K_Ng)_r=0,\quad(K_Ng)_x=\int_0^xg(s,x-s)ds.
\]

The Neumann value need not vanish. These identities preserve the requested
rows without imposing zero b,d values or a further harmonic boundary row.

Every effective point obeys `|q|+s<=r+x<=8`; the triangle is backward closed.
No physical coefficient or forcing at negative depth or outside the retained
collar is required. In particular coefficient reflections need not be smooth.

## 3. The complete C1 fixed point

Use the Banach space of triples C1 on the closed triangle, with one-sided
jets at its boundary, and norm

\[
 N(U)=\max_i\{2\|U_i\|_\infty/X,\|U_{i,x}\|_\infty,
                         \|U_{i,r}\|_\infty\}.
\]

The zero Cauchy and mixed boundary rows define a closed subspace. Completeness
follows from uniform convergence of the functions and their first partials
on this convex compact domain. Let K act componentwise as `(KD,KN,KN)`.
For continuous forcing, `N(Kg)<=X||g||`. Therefore the map

\[
 \mathcal T(U)=K(PU_x+QU_r+RU+s)
\]

preserves the subspace and has Lipschitz constant at most q above. Its Picard
sequence `U0=0`, `U(n+1)=T(U(n))` converges in C1 over the **entire** D.
It gives a unique C1 mild solution and satisfies the wave equation in
distributions. Conversely a C1 solution of that distributional wave equation
with these rows has the same retarded representation, hence the same solution.

The quantitative estimate is

\[
 N(U)\le X\|s\|_\infty/(1-q).
\]

Using `q<3/250` and `||s||<8/10^7` gives
`N(U)<6/2470000`, and `||U||<9/2470000`, which imply the displayed caps.
This is an actual continuous solution bound, not only a finite collection
of point jets or a bound on one Picard iterate.

An equivalent direct characteristic construction uses `p=Ux+Ur`,
`m=Ux-Ur` and reflection matrix `J0=diag(-1,1,1)`. The incoming p is the
integral of the forcing along its full incoming ray. The outgoing m integrates
from the initial line when `r>=x`; when `r<x`, it starts from
`m(x-r,0)=J0 p(x-r,0)`. The reflected path has total time length x.
The distributional compatibility
`D_x(p-m)=D_r(p+m)` and zero initial defect imply that integrating
`Ux=(p+m)/2` also gives `Ur=(p-m)/2`. Thus the characteristic construction
produces the same C1 field, with both boundary rows retained.

## 4. Why higher derivatives need trace terms

For physical C2 forcing write `a_+=r+x-s`, `a_-=r-x+s`.
For D let `g1(s,q)=g_r(s,|q|)` and
`g2(s,q)=sgn(q)g_rr(s,|q|)`. Set
`B=1_(x>r)g(x-r,0)`, `B1=1_(x>r)g_x(x-r,0)`. Then

\[
 K_{rr}=\tfrac12\int(g1_+-g1_-)ds-B,\qquad
 K_{xr}=\tfrac12\int(g1_++g1_-)ds+B,
\]
\[
 K_{rrr}=\tfrac12\int(g2_+-g2_-)ds+B1,\qquad
 K_{xrr}=\tfrac12\int(g2_++g2_-)ds-B1.
\]

For N instead let `g1=sgn(q)g_r(s,|q|)`, `g2=g_rr(s,|q|)` and
`Br=1_(x>r)g_r(x-r,0)`. The second formulas have no B terms; the third are

\[
 K_{rrr}=\tfrac12\int(g2_+-g2_-)ds-Br,\qquad
 K_{xrr}=\tfrac12\int(g2_++g2_-)ds+Br.
\]

For either mode, `Kxx=g+Krr`, `Kxxr=g_r+Krrr` and
`Kxxx=g_x+Kxrr`. In particular at the shell D has `Krr=-g`, while N has
`Krrr=-g_r` and `Kxxr=0`. The Neumann third trace must not be dropped.

D joins through second order if `g(0,0)=0`, and through third if also
`g_x(0,0)=0`. N joins automatically through second order, and through
third if `g_r(0,0)=0`. At positive time none of these statements requires
the physical D forcing or the N radial forcing derivative to vanish.
For example `KD(1)` is C1 but has a second-derivative jump at `x=r`;
`KN(r)` is C2 but has a third-derivative jump there. The latter is

\[
 K_N(r)=\begin{cases}x^2r/2,&r\ge x,\\
 x^3/6+xr^2/2-r^3/6,&0\le r\le x.
 \end{cases}
\]

These are explicit counterexamples to an automatic reflection-smoothing claim.

## 5. Constructive C3 convergence, not only uniform higher bounds

Assume P,Q,R,s are physical C2 on D, as supplied by the smooth retained
background and the C3 first master. The source involves the first master
through one derivative, so its two derivatives use at most master order three.
The actual compact packet also gives `s=0` on the initial shell wedge
`x+r<1`. Each Picard iterate is identically zero there, by causal closure.
Consequently every initial shell corner jet needed below vanishes. This
strong flatness assumption is used to make **each** Picard iterate C3; a
general nonflat problem can have compatible cancellations without each
separate iterate being C3, and is not the theorem proved here.

The same causal induction gives `Un=0` when `r>2+x`, because the actual
master source has that support. Its smooth front is inherited under the
retained compatible derivatives. No incoming data on an artificial interior
boundary are introduced.

The trace formulas above show by induction that each iterate is C3.
Define `Vn=partial_x Un` and `Wn=partial_x^2 Un`. The time-independent
D/N propagator obeys, for smooth physical forcing,

\[
 \partial_x Kg=S_{D/N}(g(0,\cdot))+K g_x,
\]

where S denotes the free wave with zero displacement and the stated initial
velocity. This identity follows by integration by parts in the forcing time;
it does not differentiate coefficient parity extensions in r.

Let `L=K(P partial_x+Q partial_r+R)`, the same bounded linear operator
with norm at most q. Then

\[
 V_{n+1}=LV_n+S(s_0)+K R_1(U_n),
\]
\[
 R_1(U)=P_xU_x+Q_xU_r+R_xU+s_x.
\]

Here `V0=0`; every `Vn` for `n>=1` has initial displacement zero and
initial velocity `s0`. The free term is the same for every update.

For the second time derivatives the finite offset is essential:

- `W0=0`.
- `W1` has initial displacement `s0` and velocity `s_x0`.
- **W2 and every later Wn** have initial displacement `s0` and velocity
  `P0 s0+s_x0`.

For `n>=1` one therefore has the common-data recurrence

\[
 W_{n+1}=LW_n+C(s_0)+S(P_0s_0+s_{x0})+K R_2(U_n,V_n),
\]
\[
 R_2(U,V)=2P_xV_x+2Q_xV_r+2R_xV
              +P_{xx}U_x+Q_{xx}U_r+R_{xx}U+s_{xx}.
\]

C is the free wave with the stated displacement and zero velocity.
The Neumann corner displacement derivative is
`D_r W_N(0,0)=D_r s_N(0,0)=0`: the actual source has only the D component,
so its N components vanish identically. The D displacement and velocity
corner values vanish because of the initial flat wedge. Thus these free
waves satisfy the required mixed rows and are C1; the stronger flatness also
makes the exceptional first two iterates harmless.

To prove convergence, put `theta=X/2` and

\[
 \alpha_k=\sup\max_i\sum_j(|\partial_x^kP_{ij}|+
 |\partial_x^kQ_{ij}|)+\theta\sup\max_i\sum_j|\partial_x^kR_{ij}|.
\]

Define `delta U_n=N(U(n+1)-Un)`, and define `delta V_n` and
`delta W_n` in the same way. Differences of successive iterates satisfy

\[
 \delta U_n\le q\delta U_{n-1},\quad
 \delta V_n\le q\delta V_{n-1}+X\alpha_1\delta U_{n-1},
\]
\[
 \delta W_n\le q\delta W_{n-1}
                +X(2\alpha_1\delta V_{n-1}+\alpha_2\delta U_{n-1})
 \quad(n\ge2\text{ for the W difference}).
\]

Here every delta is the weighted C1 norm. Hence the differences are bounded
by summable geometric tails times polynomials of degree at most zero, one
and two respectively. In particular **Vn and Wn converge in C1**, not merely
in value or under a uniform C3 bound. Since `Un->U` in C1 and
`Vn=partial_x Un`, their limit is `Ux`; since `Vn->Ux` in C1 and
`Wn=partial_x Vn`, its limit is `Uxx`.

It follows that `Ux` and `Uxx` are C1. The wave equation recovers the remaining
radial derivatives as continuous classical jets:

\[
 U_{rr}=U_{xx}-F,\quad U_{xrr}=U_{xxx}-F_x,
 \quad U_{rrr}=U_{xxr}-F_r,\qquad F=PU_x+QU_r+RU+s.
\]

These identities hold first in distributions, then classically because their
right sides are continuous. The recovered shell traces retain the D and N
reflection contributions listed above. Thus U is C3 up to the physical shell
and all one-sided edges of D, and solves the displayed wave equations
classically. No smooth odd/even coefficient extension or high-Sobolev solve
has entered the argument.

## 6. Explicit reusable higher-jet estimates

The same proof supplies continuous numerical majorants once the required
coefficient/source derivative bounds are entered. Let `S_k=||partial_x^k s||`
for `k=0,1,2`, and `S_r=||s_r||`. Put

\[
 B=XS_0/(1-q),\quad
 B_V=[2\|s_0\|+X(\alpha_1B+S_1)]/(1-q).
\]

Then `N(U)<=B`, `N(Ux)<=B_V`. Write
`H0=||P0 s0+s_x0||`, `G1=||partial_r s0||` and

\[
 N_{W,free}\le\max\{2\|s_0\|/X+2H0,\ G1+H0\}.
\]

The free-data derivatives here are physical half-line derivatives in scaled r.
Their compatible extensions are C1 because the initial source is flat at the
shell. Then

\[
 B_W=[N_{W,free}+X(2\alpha_1B_V+\alpha_2B+S_2)]/(1-q)
\]

bounds `N(Uxx)`. A bound for all second partials is
`B2=B_V+alpha0 B+S0`. For all third partials use

\[
 F_x^*=L_1B_V+L_0B+\alpha_1B+S_1,
\]
\[
 F_r^*=L_1B2+L_0B+\beta_1B+S_r,\quad
 \beta_1=\sup\max_i\sum_j(|P_{ij,r}|+|Q_{ij,r}|)
                    +\theta\sup\max_i\sum_j|R_{ij,r}|,
\]
\[
 B3=B_W+\max(F_x^*,F_r^*).
\]

These are bounds for derivatives of the actual constructed solution, not
derivatives of a C0 error bound. The companion checker evaluates the C1
bound only; numerical higher-jet envelopes can be obtained independently from
the coefficient/source certificate using these recurrences.

## 7. Approximation identity and parity audit

For the displacement-only free master
`Wdisp=[Qo(r+x)+Qo(r-x)]/2`, let `Po(s)=s Q(|s|)` be the odd extension of
the positive-depth profile `rQ(r)`. In particular **Po is not s Qo(s)**;
the latter is even. For any constant c, the physical half-line field

\[
 p=\frac c4[Po(r+x)+Po(r-x)-r(Qo(r+x)+Qo(r-x))]
\]

has zero Cauchy data for `r>=0`, Dirichlet value zero, and satisfies
`p_xx-p_rr=c partial_r Wdisp`. Indeed the Po free wave solves the homogeneous
equation and `L(rWdisp)=-2 partial_r Wdisp`. This proof uses only the
physical half-line; the displayed expression need not define an odd solution
on an auxiliary negative-depth line.

With `0<=Q<=1`, `||Q'||<=16`, and `X=3`, direct piecewise differentiation
gives `|p|<=1.5|c|`, `|p_x|<=24.25|c|`, `|p_r|<=24.5|c|`.
For the actual first master, the retarded first-derivative formula gives

\[
 \|U_r-\partial_r Wdisp\|
 \le\ell+\|y'\|+X\|\ell^2\mathcal W U\|.
\]

The terms are the initial-velocity rays, causal boundary lift and potential
forcing rays. No derivative of a value-only error certificate is taken.

## 8. Reproduction and finite-check limits

`check_mixed_wave_regularity.py` accepts explicit `--background`, `--master`
and `--output` paths; all output may be outside the source tree. It pins the
whole-collar actual potential/reconstruction report and the exact global
master-jet report. It imports only the Python standard library and stores
no absolute input paths or timestamps in its result.

Its exact monomial kernel integrations are differentiated independently by
tensor Lagrange weights and compared with the D/N formulas through order
three, including shell and characteristic witnesses. Compatible characteristic
joins and both counterexamples are checked. Six deliberate errors are each
rejected by explicit exceptions in normal and optimized Python. These finite
witnesses supplement the analytic Green, contraction and convergence proofs;
they do not replace those all-point arguments or certify a nonlinear problem.
