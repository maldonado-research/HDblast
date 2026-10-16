# First relative-harmonic vector: exact operator and constructive boundary map

16 October 2026. Internal mathematical derivation. This file proves the displayed
operator identities and a constructive, conditional C1 vector-wave theorem on a
specified triangle. It does not solve the nonlinear Einstein-scalar boundary
problem or compute transformed second-order sources. Authentication of the
actual background and master-wave bounds is supplied by the separate coefficient
and solver certificates. Higher regularity is treated separately in `regularity`.

## 1. Conventions and geometric operator

Use

\[
g_0=dy^2+\rho(y)^2\bar g,\qquad
\bar g=-dt^2+\cosh^2t\,\kappa_{S^3},\qquad a=\rho_y/\rho,
\qquad dy=-\rho\,dz.
\]

The shell is `z=0`. The longitudinal first metric has
`h_yy=-4 Psi`, `h_ya=0`, `h_ab=2 Psi (g0)_ab`, with
`Psi=u(y,t)Y` and `Delta_kappa Y=-8Y`. For the relative harmonic vector

\[
C^A=g^{BC}(\Gamma[g]^A{}_{BC}-\Gamma[g_0]^A{}_{BC}),\qquad
DC[h]_A=\nabla^B h_{BA}-\tfrac12\nabla_A\operatorname{tr}_0h,
\]

the trace is `4 Psi` and

\[
C^{\rm long}_{1,y}=-6(\Psi_y+4a\Psi),\qquad C^{\rm long}_{1,a}=0.
\]

Writing `h_wave=h_long+L_xi g0`, covariance gives

\[
DC[\mathcal L_\xi g_0]=(\Box_0+\operatorname{Ric}_0)\xi,
\qquad (\Box_0+\operatorname{Ric}_0)\xi=-C^{\rm long}_1. \tag{1}
\]

For contravariant components `X=xi^y`, `V^a=xi^a`, the exact warped split is

\[
(P\xi)^y=X_{yy}+4aX_y+\rho^{-2}\bar\Box X
-2a\overline{\operatorname{div}}V-(4a_y+8a^2)X, \tag{2}
\]
\[
(P\xi)^a=V^a_{yy}+6aV^a_y+
\rho^{-2}(\bar\Box V^a+3V^a)+2a\rho^{-2}\bar\nabla^aX. \tag{3}
\]

Here the barred wave operator on V is the **rough vector** operator. A direct
derivation, independent of a guessed vector-wave formula, applies DC to

\[
h_{yy}=2X_y,\quad h_{ya}=\partial_aX+\rho^2\partial_yV_a,\quad
h_{ab}=\rho^2(2aX\bar g_{ab}+\bar\nabla_aV_b+\bar\nabla_bV_a),
\]

where V_a is lowered with the barred metric. The trace is
`2X_y+8aX+2 div_bar V`. Its radial divergence before subtracting half the
trace derivative is

\[
2X_{yy}+\rho^{-2}\bar\Box X+(\operatorname{div}V)_y
+8aX_y-8a^2X-2a\operatorname{div}V.
\]

The lower tangential divergence is

\[
\partial_aX_y+\rho^2V_{a,yy}+6a\rho^2V_{a,y}
+6a\partial_aX+\bar\Box V_a+\partial_a\operatorname{div}V+3V_a.
\]

Subtracting the corresponding trace derivatives gives (2) and (3). In
particular no extra tangential warp potential survives.

## 2. Closed scalar angular sector and original shell rows

Take raw coordinate amplitudes

\[
\xi^y=f(y,t)Y,\qquad \xi^t=b(y,t)Y,\qquad
\xi^i=d(y,t)D_\kappa^iY. \tag{4}
\]

The symbol f here is the radial gauge-vector amplitude, not the second spatial
metric amplitude used in earlier initial-data work. It has **not** been divided
by rho. Let `lambda=8`, `H=tanh t`, `S=sech^2 t`. Then

\[
\overline{\operatorname{div}}V=(b_t+3Hb-\lambda d)Y,
\]
\[
(\bar\Box+3)V^t=
[-b_{tt}-3Hb_t+(3+3H^2-\lambda S)b-2\lambda Hd]Y,
\]
\[
(\bar\Box+3)V^i=
[-d_{tt}-5Hd_t+(4-\lambda)Sd+2HSb]D_\kappa^iY.
\]

For example the lower components of `L_V bar g` are
`h_tt=-2b_t Y`, `h_ti=(cosh^2(t)d_t-b)D_iY`, and
`h_ij=2cosh^2(t)[Hb kappa_ij Y+d Hess_ij Y]`.
Their trace is `2b_t+6Hb-2lambda d`.
Using `D^j Hess_ji Y=(2-lambda)D_iY` in their divergence gives exactly the
two displayed vector formulas. The checker also independently verifies the
scalar-gradient identity and the de Sitter boost Killing mode.

The three equations in the original radial coordinate are

\[
\begin{split}
f_{yy}+4af_y+\rho^{-2}(-f_{tt}-3Hf_t-\lambda Sf)
-2a(b_t+3Hb-\lambda d)\qquad\\
-(4a_y+8a^2)f=6(u_y+4au),
\end{split} \tag{5}
\]
\[
b_{yy}+6ab_y+\rho^{-2}[-b_{tt}-3Hb_t+
(3+3H^2-\lambda S)b-2\lambda Hd-2af_t]=0, \tag{6}
\]
\[
d_{yy}+6ad_y+\rho^{-2}[-d_{tt}-5Hd_t+
(4-\lambda)Sd+2HSb+2aSf]=0. \tag{7}
\]

The original shell conditions are

\[
f|_b=0,\qquad b_y|_b=d_y|_b=0. \tag{8}
\]

They are the scalar-sector form of `xi^y=0`, `partial_y xi^a=0`.
Indeed

\[
(\mathcal L_\xi g_0)_{ya}=\partial_aX+\rho^2\bar g_{ab}\partial_yV^b,
\qquad (\partial_y-2a)\xi_a=\rho^2\bar g_{ab}\partial_yV^b.
\]

Thus the shell stays fixed, and the mixed metric boundary components remain
zero at first order. Tangential values of the vector remain free: they
reparametrize the shell and its time label. There are five geometric rows
(one normal Dirichlet and four tangential Neumann), reduced to the three rows
(8) in this invariant scalar sector. No separate boundary condition C=0 or
zero tangential vector value is being added.

The initial data are `f=b=d=f_t=b_t=d_t=0`. Consequently
`f_tt=-6rho^2(u_y+4au)` initially and `b_tt=d_tt=0`, matching the previously
computed initial harmonic shift acceleration.

## 3. Exact inward and scaled coefficient matrices

Because

\[
\partial_y=-\rho^{-1}\partial_z,\qquad
\partial_y^2=\rho^{-2}\partial_z^2+(a/\rho)\partial_z,
\]

the radial drifts in (5) and (6),(7) become respectively `-3rho a` and
`-5rho a` after multiplication by rho squared. For U=(f,b,d)^T,

\[
U_{tt}-U_{zz}=A U_t+B U_z+D U+F, \tag{9}
\]

\[
A=\begin{pmatrix}-3H&-2a\rho^2&0\\-2a&-3H&0\\0&0&-5H\end{pmatrix},
\qquad B=-a\rho\,\operatorname{diag}(3,5,5),
\]
\[
D=\begin{pmatrix}
-\lambda S-\rho^2(4a_y+8a^2)&-6a\rho^2H&2\lambda a\rho^2\\
0&3+3H^2-\lambda S&-2\lambda H\\
2aS&2HS&(4-\lambda)S
\end{pmatrix},
\qquad F=\begin{pmatrix}6\rho u_z-24a\rho^2u\\0\\0\end{pmatrix}. \tag{10}
\]

Let `ell=1/1000000`, `x=t/ell`, `r=z/ell`, with

\[
\mathcal D=\{0\le x\le3,\quad r\ge0,\quad x+r\le8\}.
\]

All amplitudes remain raw coordinate components. The scaled equation is

\[
U_{xx}-U_{rr}=\mathsf A U_x+\mathsf B U_r+\mathsf D U+\mathsf F,
\quad \mathsf A=\ell A,\quad \mathsf B=\ell B,
\quad \mathsf D=\ell^2D,\quad \mathsf F=\ell^2F. \tag{11}
\]

Here `H=tanh(ell x)`, `S=sech^2(ell x)`, and radial coefficients are
evaluated at `z=ell r`. The shell rows are exactly

\[
f(x,0)=0,\qquad b_r(x,0)=d_r(x,0)=0. \tag{12}
\]

For the actual master reconstruction, use a notation distinct from the old
master drift coefficient: define `R=alpha/rho^2`, `h=U_potential'/p`,
`vartheta(t)=cosh(t)^(-3/2)`, and `u=R vartheta w`. The time factor is part
of the retained conjugated-master normalization, even though it equals one
on the initial slice. Since `alpha_z=rho(5a/2-h)alpha`,

\[
R_z=\rho(9a/2-h)R,\qquad
\mathsf F_1=6\ell\rho R\vartheta w_r+
6\ell^2\rho^2R\vartheta(a/2-h)w. \tag{13}
\]

The old shell report's symbol `C=rho(2h-5a)` is **not** R. Formula (13) uses
the actual master w defined by its retarded equation; replacing w by a free
wave comparator would change the forced gauge problem.

## 4. All-point coefficient sizes and domain premise

Suppose the actual branch on `0<=z<=8ell` satisfies

\[
129\le\rho\le131,\qquad |a|\le14/125,
\qquad |p|\le7/50000,\qquad a_y=-\rho^{-2}-p^2/3.
\]

Together with `|H|<=3ell` and `0<S<=1`, exact rational row-sum estimates give

| Matrix | Valid uniform row-sum cap |
|---|---:|
| mathsf A | `3844064009/1000000000000 < 0.003845` |
| mathsf B | `917/12500000 = 0.00007336` |
| mathsf D | `337883156437144543/10400625000000000000000000 < 3.3e-8` |

The standalone checker proves these implications with exact rationals. The
older `PRIOR_ANALYTIC_SHELL_BOX_DELTA_1E3.json` alone proves only depth `8e-10`.
The required `8e-6` depth is supplied by the later independent first-exit
extension, reproduced in `coefficients/check_global_gauge_coefficients.py`.
Its proof starts from the retained actual shell box and checks the full
`8ell` flow remains strictly inside the same outer rectangle. This domain
extension, not the older short-collar prose, is the premise of the application
here. The coefficient certificate also provides sharper caps and the
derivative/source bounds needed for the actual solution.

## 5. Constructive C1 theorem with the actual shell rows

The following elementary theorem applies to arbitrary continuous coefficient
matrices and continuous forcing on D. Assume their row-sum bounds are
`A0,B0,D0`, and let `T=3`. Put

\[
q_*=T(A_0+B_0+TD_0/2)<1. \tag{14}
\]

Then (11), with zero Cauchy data and (12), has a unique C1 solution whose
wave equation holds in distributions. It is constructed by a uniformly
convergent characteristic iteration, and no artificial outer boundary is
required. With the displayed caps even the larger number

\[
q=T(A_0+B_0+TD_0)
=\frac{13581557244358312144543}{1155625000000000000000000}
<3/250 \tag{15}
\]

is below one.

Here is a proof including the boundary trace. Let

\[
p=U_x+U_r,\quad m=U_x-U_r,\quad J=\operatorname{diag}(-1,1,1),
\quad G=\tfrac12(\mathsf A+\mathsf B)p+
\tfrac12(\mathsf A-\mathsf B)m+\mathsf D U+\mathsf F.
\]

Given a continuous triple (U,p,m), define its successor by

\[
p^+(x,r)=\int_0^xG(s,r+x-s)\,ds, \tag{16}
\]
\[
m^+(x,r)=
\begin{cases}
\displaystyle\int_0^xG(s,r-x+s)\,ds,&r\ge x,\\[3pt]
\displaystyle Jp^+(x-r,0)+\int_{x-r}^xG(s,r-x+s)\,ds,&r<x,
\end{cases} \tag{17}
\]
\[
U^+(x,r)=\tfrac12\int_0^x(p^++m^+)(s,r)\,ds. \tag{18}
\]

All paths stay in D. The reflected path in (17) has lengths `(x-r)+r=x`;
J has norm one. Hence, for `M=||G||_infinity`,

\[
|p^+(x,r)|,|m^+(x,r)|\le xM,\qquad |U^+(x,r)|\le x^2M/2.
\]

Use the complete space of continuous triples, with norm

\[
\|(U,p,m)\|_*=
\max\{2\|U\|_\infty/T,\|p\|_\infty,\|m\|_\infty\}.
\]

The source differences obey
`||delta G|| <= (A0+B0+TD0/2)||delta(U,p,m)||_*`; thus (16)-(18)
contract by at most q_*. Starting from zero gives a unique continuous fixed
point and explicit geometric remainder bounds. In particular

\[
\|U_x\|_\infty,\|U_r\|_\infty
\le\frac{T\|\mathsf F\|_\infty}{1-q_*},\qquad
\|U\|_\infty\le\frac{T^2\|\mathsf F\|_\infty}{2(1-q_*)}. \tag{19}
\]

To justify the derivatives, the continuous fixed-point fields satisfy
`p_x-p_r=G` and `m_x+m_r=G` distributionally in the open triangle. The
characteristic seam `r=x` contributes no delta because (17) is continuous
there and `p^+(0,0)=0`. Subtraction gives
`partial_x(p-m)=partial_r(p+m)`. Using the zero initial values and (18),
`U_x=(p+m)/2` and `U_r=(p-m)/2` as distributions. Both right-hand sides are
continuous, so U is C1, with continuous one-sided derivatives up to the shell
and the other edges. The sum of the transport equations gives (11).

At the shell, (17) gives `m=Jp`. The Dirichlet component therefore has
`f_x=0`, and its zero initial value gives `f=0`; the Neumann components have
`b_r=d_r=0`. Conversely, a C1 distributional solution gives these transport
relations and hence the same characteristic integrals, proving uniqueness.
No data are missing at `x+r=8`: all past characteristics used here remain
inside that edge. The theorem is constructive even without differentiable
coefficients; the separate regularity argument is needed before identifying
the gauge equation pointwise at second-derivative level.

## 6. Reflection, higher jets and limits of the conclusion

A full-line parity extension would set `U(x,-r)=J U(x,r)` and must use

\[
\mathsf A(x,-r)=J\mathsf A(x,r)J,\quad
\mathsf B(x,-r)=-J\mathsf B(x,r)J,\quad
\mathsf D(x,-r)=J\mathsf D(x,r)J,\quad
\mathsf F(x,-r)=J\mathsf F(x,r). \tag{20}
\]

These extensions generally jump at the shell. In particular the nonzero
diagonal radial drift becomes odd, as do the Dirichlet-to-Neumann entries of
mathsf A. The Dirichlet forcing need not vanish at the shell after contact.
Thus differentiating a reflected forcing as if it were smooth loses trace
terms. Equations (16)-(18) retain the reflected trace explicitly and do not
make this assumption. The resulting U has a C1 parity extension; its second
normal derivatives may jump across the reflected shell without any loss of
one-sided physical regularity.

For C2/C3 conclusions, use the separate time-differentiated construction with
the actual coefficient/source jets and corner conditions. It keeps the
positive-time boundary forcing terms when recovering radial derivatives from
the equation. The flat initial packet collar supplies the required corner
compatibilities. Flatness initially does not impose vanishing Dirichlet
forcing for all later shell times. With C2 regularity established, (1) holds
classically and `C1[h_wave]=0`; with only the theorem above it holds in the
distributional sense.

A vector tangent to the shell preserves the physical first-order shell laws
by covariance, while allowing tangential and clock transport. This does not
retain a longitudinal second-order forcing formula after gauge conversion.
The transformed metric and scalar, their second-order source terms, the
second-order evolution, a nonlinear remainder, and an observational response
require additional work.

One exact preparation for that later work follows from covariance of the
**unreduced natural** Einstein-scalar residual E. Assume sufficiently smooth
Phi0, h and xi so that every term below is classical; without that premise,
the displayed calculation is only a formal identity. Regard Phi as the
complete metric-scalar pair, put `h'=h+L_xi Phi0` and
`G=2L_xi h+L_xi^2 Phi0`, and use second parameter derivatives (the coefficient
of epsilon squared is one half the second jet). Twice differentiating
`E(exp(epsilon xi)^*(Phi0+epsilon h)) = exp(epsilon xi)^*E(Phi0+epsilon h)`
gives

\[
Q(h',h')=Q(h,h)-E'[G],\qquad Q=E''[\Phi_0], \tag{21}
\]

provided `E(Phi0)=0` and `E'[h]=0`. Without these assumptions the right side
also contains `2L_xi E'[h]+L_xi^2 E(Phi0)`. This identity concerns the natural
unreduced residual. A reduced harmonic residual with fixed reference data is
not automatically the same natural operator. Equation (21) neither computes
its transformed forcing nor supplies a second gauge vector or a second-order
solution.

The quantitative certificate currently controls xi only through C3. For
smooth background and h, this suffices to form the direct second-order
expression `Q(h',h')`, because h' uses one derivative of xi. In contrast, G
contains second derivatives of xi, so a separate classical evaluation or
termwise bound of `E'[G]` can require fourth derivatives of xi. No C4 bound
has been certified here. The present C3 caps therefore do not separately
bound `E'[G]`, and (21) is not yet a shortcut for transformed forcing norms;
that use needs additional regularity estimates or an explicitly proved
derivative cancellation.

## 7. Portable evidence

`check_gauge_operator.py` uses exact sparse Laurent polynomials to check 23
all-point identities, including two independent tensor-divergence reductions,
the de Sitter gradient and Killing checks, the differentiated inward
coordinate conversion, all three matrix rows, source reconstruction, and
boundary parity. It also checks exact rational matrix caps and eight
mutation-specific nonzero/failing witnesses. No numerical sample is used as
a proof of an operator identity.

Normal and optimized runs must produce byte-identical JSON. The separate
control runner invokes each intentional mutation in both modes, requiring
the named failure and absence of output. An existing output path is also
rejected. These tests establish the finite algebra and reported bounds;
the mathematical construction is the argument in sections 5 and 6.
