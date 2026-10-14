# Complete quadrupole-generated l=4 scalar bulk equations

14 October 2026. Independent bulk continuation of the fixed Einstein–scalar
model. This note supplies the complete quadratic residual and the general
linear angular-scalar operator in a stated gauge. It does not establish a
scalar shell boundary estimate, a full second-order contact solution, or an
exact nonlinear family. The accompanying direct-geometry calculation was
performed independently of the other scalar-source derivation.

## 1. Model, harmonic, and derivative convention

The controlling model and first reconstruction are those of
`HDBLAST_REGULAR_CONE_EVOLUTION_20261013/manuscript/REGULAR_CONE_EVOLUTION.md`.
The normalization is

\[
\mathscr E_{AB}=R_{AB}-\phi_{,A}\phi_{,B}-\tfrac23U(\phi)g_{AB},
\qquad \mathscr F=\Box_g\phi-U'(\phi).
\tag{1}
\]

Write

\[
g_0=dy^2+\rho(y)^2[-dt^2+c(t)^2\kappa],\qquad
\phi_0=\phi_0(y),\quad c=\cosh t,\quad h=\tanh t,
\quad a=\rho'/\rho,\quad p=\phi_0'.
\tag{2}
\]

Here \(\kappa\) is the unit round metric on \(S^3\). The actual background
solves (1), in particular

\[
a'=-\rho^{-2}-p^2/3,\qquad
 a^2=\rho^{-2}+p^2/12-U/6,\qquad p'+4ap=U'.
\tag{3}
\]

The source formulas below do **not** require substituting (3) or the
first-order field equations. Their source-compatibility interpretation does.
All potential derivatives in this note are evaluated at \(\phi_0(y)\).

The specified first-order longitudinal representative is

\[
(h_1)_{yy}=-4\Psi,\quad (h_1)_{y\mu}=0,\quad
(h_1)_{\mu\nu}=2\rho^2\Psi\gamma_{\mu\nu},\quad
\Psi=u(y,t)Y,\quad \chi=v(y,t)Y,
\]
\[
Y=n_1^2-n_2^2,\qquad v=-\frac{3(u_y+2au)}p.
\tag{4}
\]

The selected degree-four scalar harmonic is

\[
Q=Y^2-\frac{n_1^2+n_2^2}{2}+\frac1{12},\qquad
\Delta_\kappa Q=-24Q.
\tag{5}
\]

Throughout, \(g_\epsilon=g_0+\epsilon h_1+\epsilon^2h_2/2+\cdots\) and
\(\phi_\epsilon=\phi_0+\epsilon\chi+\epsilon^2w_2/2+\cdots\).
Define \(S=D^2\mathscr E[(h_1,\chi)^2]\) and
\(S_\phi=D^2\mathscr F[(h_1,\chi)^2]\). The second equations are
\(D\mathscr E[h_2,w_2]=-S\) and \(D\mathscr F[h_2,w_2]=-S_\phi\).
In particular, every displayed source is twice the coefficient of
\(\epsilon^2\) in the residual along the **linear** path.

## 2. General scalar gauge and its Regge–Wheeler representative

Let the Lorentzian two-dimensional base metric and areal radius be

\[
q_{ab}dx^a dx^b=dy^2-\rho^2dt^2,\qquad R=\rho c,
\qquad r_a=D_a\log R,
\tag{6}
\]

where \(D_a\) is the connection of \(q\); angular derivatives are written
\(\nabla_i\). For a scalar harmonic of eigenvalue \(-\lambda\), put

\[
\mathbb H_{ij}=\nabla_i\nabla_jQ+\frac\lambda3\kappa_{ij}Q.
\tag{7}
\]

The general scalar metric has seven amplitudes:

\[
h_{ab}=f_{ab}Q,\quad h_{ai}=R f_a\nabla_iQ,\quad
h_{ij}=R^2(K\kappa_{ij}Q+G\mathbb H_{ij}),\quad w_2=wQ.
\tag{8}
\]

For a second-order gauge vector with covariant components
\(\xi_a=X_aQ\), \(\xi_i=R^2L\nabla_iQ\), the convention
\(h_2\mapsto h_2+\mathcal L_\xi g_0\) gives

\[
\begin{aligned}
\delta f_{ab}&=D_aX_b+D_bX_a,&
R\delta f_a&=X_a+R^2D_aL,\\
\delta G&=2L,&
\delta K&=2r^aX_a-\frac{2\lambda}3L,&
\delta w&=X^ap_a.
\end{aligned}\tag{9}
\]

For \(l\ge2\), choose \(L=-G/2\) and
\(X_a=-Rf_a-R^2D_aL\). This gives the complete scalar
Regge–Wheeler gauge

\[
\boxed{f_a=G=0,\qquad
h_{ab}=f_{ab}Q,\quad h_{ij}=R^2K\kappa_{ij}Q,\quad w_2=wQ.}
\tag{10}
\]

No base lapse or shift is deleted. The remaining five amplitudes are the
three components of \(f_{ab}\), \(K\), and \(w\). There is no residual
scalar coordinate freedom preserving (10) for \(l\ge2\). The exceptional
\(l=0,1\) gauges are not covered by this argument.

This is an **order-two** gauge choice with the first fields (4) fixed, so
it does not change \(S\): the extra linear residual is
\(D\mathscr E[\mathcal L_\xi g_0,\mathcal L_\xi\phi_0]
=\mathcal L_\xi\mathscr E_0=0\). A change of first-order gauge must
instead carry along the quadratic conversion
\(2\mathcal L_{\xi_1}h_1+\mathcal L_{\xi_1}^2g_0\), its scalar analogue,
and the changed quadratic source. No source here is claimed universal
under that different operation.

## 3. The full linear scalar operator

The following formulas hold for any \(\lambda=l(l+2)\), \(l\ge2\), on
the warped background (6). Use base indices raised by \(q\), and define

\[
f=q^{ab}f_{ab},\qquad J=f+K,\qquad
A^a{}_{bc}=\tfrac12(D_bf^a{}_c+D_cf^a{}_b-D^af_{bc}),
\]
\[
B^a=r^bf^a{}_b-r^aK-\tfrac12D^aK.
\tag{11}
\]

Write the linear Ricci tensor as
\(\delta R_{ab}=\mathcal R_{ab}Q\),
\(\delta R_{ai}=\mathcal R_a\nabla_iQ\), and
\(\delta R_{ij}=\mathcal R_L\kappa_{ij}Q+\mathcal R_T\mathbb H_{ij}\).
Then

\[
\boxed{\begin{aligned}
\mathcal R_{ab}={}&D_cA^c{}_{ab}-\tfrac12D_aD_bf+3r_cA^c{}_{ab}
 -\tfrac32(D_aD_bK+r_aD_bK+r_bD_aK)
 +\frac\lambda{2R^2}f_{ab},\\
\mathcal R_a={}&\tfrac12\left[D_bf^b{}_a-D_af-2D_aK
                  +r_bf^b{}_a+r_af\right],\\
\mathcal R_L={}&R^2\left[D_aB^a+3r_aB^a-\tfrac12r^aD_aJ\right]
                 +\frac\lambda2K+\frac\lambda6J,\\
\mathcal R_T={}&-\tfrac12J.
\end{aligned}}\tag{12}
\]

The corresponding trace-reversed Einstein coefficients are

\[
\boxed{\begin{aligned}
\mathcal E_{ab}&=\mathcal R_{ab}-p_aD_bw-p_bD_aw
                  -\tfrac23(Uf_{ab}+U'wq_{ab}),\\
\mathcal E_a&=\mathcal R_a-p_aw,\\
\mathcal E_L&=\mathcal R_L-\tfrac23R^2(UK+U'w),\\
\mathcal E_T&=-J/2,
\end{aligned}}\qquad p_a=(p,0).
\tag{13}
\]

The complete scalar-field operator is

\[
\boxed{\begin{aligned}
\mathcal F={}&\left(D^aD_a+3r^aD_a-\frac\lambda{R^2}-U''\right)w
 -f^{ab}D_aD_b\phi_0\\
&-\left[D_af^{ab}+3r_af^{ab}
             -\tfrac12D^b(f+3K)\right]p_b.
\end{aligned}}\tag{14}
\]

These formulas include the mixed base metric component. They contain no
inverse background scalar derivative and remain operator identities even
where \(p=0\); the inherited reconstruction (4) separately needs \(p\ne0\).

A short derivation keeps the connection variation rather than commuting
second derivatives. Its nonzero blocks in the gauge (10) are

\[
\begin{array}{ll}
\delta\Gamma^a{}_{bc}=A^a{}_{bc}Q,&
\delta\Gamma^i{}_{ab}=-f_{ab}\nabla^iQ/(2R^2),\\
\delta\Gamma^a{}_{bi}=\tfrac12f^a{}_b\nabla_iQ,&
\delta\Gamma^i{}_{aj}=\tfrac12D_aK\,\delta^i_jQ,\\
\delta\Gamma^a{}_{ij}=R^2B^a\kappa_{ij}Q,&
\delta\Gamma^i{}_{jk}=\tfrac K2
(\delta^i_j\nabla_kQ+\delta^i_k\nabla_jQ-\kappa_{jk}\nabla^iQ).
\end{array}\tag{15}
\]

Their trace is \(\delta\Gamma^A{}_{AB}=\partial_B[(f+3K)Q]/2\).
Substitution into
\(\delta R_{AB}=\nabla_C\delta\Gamma^C{}_{AB}
-\nabla_B\delta\Gamma^C{}_{AC}\) proves (12). In (14), the apparent
terms proportional to \(K r^ap_a\) cancel between the varied Hessian
contraction and the varied connection contraction.

## 4. Complete l=4 quadratic sources

Set

\[
\mathcal B_2z=-z_{tt}-3hz_t-8c^{-2}z,
\]
\[
\boxed{\begin{aligned}
L={}&2u\mathcal B_2u+2(-u_t^2-4c^{-2}u^2)\\
&-\rho^2\left[4u_y^2+4uu_{yy}+52auu_y
 +24(a'+4a^2)u^2+\tfrac13U''v^2+\tfrac43U'uv\right].
\end{aligned}}\tag{16}
\]

Here \(L\) is an auxiliary source expression, not a metric amplitude.
Decompose the selected harmonic of the residual as

\[
S_{ab}^{(4)}=s_{ab}Q,\quad S_{ai}^{(4)}=s_a\nabla_iQ,\quad
S_{ij}^{(4)}=s_L\kappa_{ij}Q+s_T\mathbb H_{ij},\quad
S_\phi^{(4)}=s_\phi Q.
\tag{17}
\]

The full eight scalar coefficients are

\[
\boxed{\begin{aligned}
s_{yy}={}&-8u_y^2+16uu_{yy}-32auu_y
 +\rho^{-2}(8uu_{tt}+24huu_t-16u_t^2)\\
&-2v_y^2-\tfrac23U''v^2+\tfrac{16}3U'uv,\\
s_{yt}={}&12uu_{yt}-48auu_t-2v_yv_t,\\
s_{tt}={}&24uu_{tt}+12u_t^2-2v_t^2-2L,\\
s_y={}&6uu_y-24au^2-vv_y,\\
s_t={}&18uu_t-12hu^2-vv_t,\\
s_L={}&c^2(2L-24huu_t)-80u^2+\tfrac83v^2,\\
s_T={}&8u^2-\tfrac23v^2,\\
s_\phi={}&12u_yv_y+48puu_y+8u(v_{yy}+4av_y)
 +32u^2(p'+4ap)-4\rho^{-2}u\mathcal B_2v-U'''v^2.
\end{aligned}}\tag{18}
\]

In \(s_{yy}\), angular terms cancel after the degree-four projection;
this cancellation must not be imposed before doing the projection. The
scalar source may use \(p'+4ap=U'\) on the actual background.

Equations (12)–(14), with \(\lambda=24\), and the eight right-hand sides
(18) give the explicit complete sourced scalar system:

\[
\boxed{\mathcal E_{ab}=-s_{ab},\qquad
\mathcal E_a=-s_a,\qquad
\mathcal E_L=-s_L,\qquad
\mathcal E_T=-s_T,\qquad \mathcal F=-s_\phi.}
\tag{19}
\]

This is a closed constrained differential system: no scalar residual
component or forcing amplitude is left unspecified. It is not yet a
proved hyperbolic reduction with a constraint-preserving shell boundary
prescription.

## 5. Independent derivation of the quadratic terms

Along the linear first-order metric path write

\[
N^2=1-4\epsilon\Psi,\quad
q_{\mu\nu}=\rho^2e^{2\omega}\gamma_{\mu\nu},\quad
\omega=\tfrac12\log(1+2\epsilon\Psi),\quad
\nu=\log N=\tfrac12\log(1-4\epsilon\Psi).
\]

The exact radial ADM identities in this zero-shift chart are

\[
\begin{aligned}
R_{yy}={}&-4[a'+a^2+\omega_{yy}+2a\omega_y+\omega_y^2
                   -(a+\omega_y)\nu_y]\\
&-\rho^{-2}e^{2\nu-2\omega}
 [\Box_\gamma\nu+2d\omega\mathbin\cdot d\nu+|d\nu|_\gamma^2],\\
R_{y\mu}={}&-3[\omega_{y\mu}-(a+\omega_y)\nu_\mu],\\
R_{\mu\nu}={}&3\gamma_{\mu\nu}
 -2\omega_{;\mu\nu}-\nu_{;\mu\nu}
 +2\omega_\mu\omega_\nu-\nu_\mu\nu_\nu
 +\omega_\mu\nu_\nu+\omega_\nu\nu_\mu\\
&+\gamma_{\mu\nu}\left[-\Box_\gamma\omega
 -2|d\omega|_\gamma^2-d\omega\mathbin\cdot d\nu
 -\rho^2e^{2\omega-2\nu}
 \{a'+\omega_{yy}-(a+\omega_y)\nu_y+4(a+\omega_y)^2\}\right],\\
\Box_g\phi={}&e^{-2\nu}[\phi_{yy}+(4a+4\omega_y-\nu_y)\phi_y]
 +\rho^{-2}e^{-2\omega}
 [\Box_\gamma\phi+(2d\omega+d\nu)\mathbin\cdot d\phi].
\end{aligned}\tag{20}
\]

The products use the Lorentzian metric \(\gamma\). Expanding
\(\omega=\epsilon\Psi-\epsilon^2\Psi^2\) and
\(\nu=-2\epsilon\Psi-4\epsilon^2\Psi^2\) gives, in particular,

\[
[\epsilon^2]R_{y\mu}=6\Psi\Psi_{y\mu}-24a\Psi\Psi_\mu,
\qquad
[\epsilon^2]R_{\mu\nu}\equiv
12\Psi\Psi_{;\mu\nu}+6\Psi_\mu\Psi_\nu
\pmod{\gamma_{\mu\nu}}.
\tag{21}
\]

The sphere projections used in (18) are

\[
[Y^2]_4=Q,\quad[|\nabla Y|^2]_4=-4Q,
\quad[Y\nabla_iY]_4=\tfrac12\nabla_iQ,
\]
\[
[(\nabla Y\otimes\nabla Y)^{\mathrm{TF}}]_{4,\mathrm{scalar}}
=\tfrac13\mathbb H_{ij}.
\tag{22}
\]

For the last identity the full decomposition is

\[
(\nabla Y\otimes\nabla Y)^{\mathrm{TF}}
=\tfrac12(\operatorname{Hess}F_2)^{\mathrm{TF}}
 +\tfrac13\mathbb H-\tfrac43\mathcal T,
\quad F_2=n_1^2+n_2^2-\tfrac12.
\tag{23}
\]

The accompanying checker builds \(\mathcal T\) directly by restricting
the ambient Weyl tensor with sectional entries
\(C_{1212}=C_{3434}=-1/3\),
\(C_{1313}=C_{1414}=C_{2323}=C_{2424}=1/6\). It verifies all six
components of (23) as exact rational identities in stereographic
coordinates, as well as the scalar eigenvalue and product identities.
The TT part of (21) and the matter term is consequently
\((16u^2+8v^2/3)\mathcal T\), agreeing with the accepted tensor source.

Every mixed angular one-form in (20) is proportional to
\(Y\nabla_iY\), an exact gradient. Equation (23) contains only
scalar-derived and TT angular tensors. Thus the complete quadratic **bulk**
source in this specified first gauge has zero coexact vector projection.
This statement alone says nothing about vector initial corrections or
boundary gauge conversion.

## 6. A definite failure of the unchanged one-master reconstruction

The angular trace-free row of (19) is the algebraic relation

\[
\boxed{J=f+K=16u^2-\tfrac43v^2=:C.}\tag{24}
\]

The inherited longitudinal metric reconstruction, applied to any putative
second-order master amplitude \(z\), has
\(f_{yy}=-4z\), \(f_{tt}=-2\rho^2z\), \(f_{yt}=0\), \(K=2z\).
It therefore has \(J=0\) identically. It cannot solve (19) wherever
\(8u^2-2v^2/3\ne0\). This is a nonzero residual in an explicitly
retained Einstein component, not just a count of unknowns.

For a nonzero smooth compact radial displacement on a regular collar,
the obstruction cannot vanish identically on the slice. Indeed, if
\(v^2=12u^2\) everywhere and (4) holds, then

\[
|u_y|\le\left(2|a|+\frac2{\sqrt3}|p|\right)|u|.
\]

The coefficient is bounded on the compact collar. Integrating from a
radial point outside the support and applying Gronwall forces \(u=0\).
Hence a nontrivial compact displacement has \(C\ne0\) somewhere. For
pure-velocity initial master data, \(u=v=0\) on the initial slice is a
real exception at that time; the statement applies to slices on which
the displacement becomes nonzero. A more general sourced master
construction with additional metric reconstruction terms is not ruled
out by this argument.

To make the retained lapse and shift explicit, write
\(f_{yy}=N\), \(f_{yt}=B\), \(f_{tt}=\rho^2T\). Then

\[
T=N+K-C,
\]
\[
\boxed{\begin{aligned}
N_y-K_y+3aN-\rho^{-2}(B_t+hB)-2pw&=C_y-2s_y,\\
N_t+2K_t&=B_y+2aB-hN-2hK+2hC+2s_t.
\end{aligned}}\tag{25}
\]

These are the two mixed Einstein rows, not permission to discard the
three base rows, the angular trace row, or the scalar equation. They
show explicitly which metric fields a one-function substitution misses.

## 7. Consistency, constraints, and physical degrees of freedom

The exact off-shell identity is

\[
\nabla^A\mathscr E_{AB}-\tfrac12\nabla_B\operatorname{tr}_g\mathscr E
=-\mathscr F\,\phi_{,B}.
\tag{26}
\]

On an on-shell background with on-shell first fields, its second
variation says that both the linear operator and the quadratic source
obey the same background identity. There are no omitted connection
variation terms: those multiply \(\mathscr E_0\) or
\(D\mathscr E[h_1,\chi]\), which vanish under the stated hypotheses.

For any scalar residual coefficients \((e_{ab},e_a,e_L,e_T,e_\phi)\),
put \(\tau=q^{ab}e_{ab}+3e_L/R^2\). The three projected identities are

\[
\boxed{\begin{aligned}
D^ae_{ab}+3r^ae_{ab}-\frac\lambda{R^2}e_b
 -\frac{3r_b}{R^2}e_L-\tfrac12D_b\tau&=-p_be_\phi,
\qquad b=y,t,\\
D^ae_a+3r^ae_a+\frac1{R^2}
 \left[e_L+\left(2-\frac{2\lambda}3\right)e_T\right]
 -\tfrac12\tau&=0.
\end{aligned}}\tag{27}
\]

The coefficient in the second line is \(-14\) at \(\lambda=24\).
Thus (19) is source-consistent when its first-order inputs really solve
the inherited problem. Arbitrary off-shell functions \(u,v\) still give
the correct residual (18), but need not give a compatible linear source.

For time evolution the actual Einstein residual is the trace reversal
\(\mathcal G_{AB}=\mathscr E_{AB}-g_{AB}\operatorname{tr}\mathscr E/2\).
The scalar harmonic Hamiltonian and two momentum constraints are,
up to nonzero background lapse factors,

\[
\tfrac12(e_{yy}+\rho^{-2}e_{tt}+3e_L/R^2)=0,
\qquad e_{ty}=0,\qquad e_t=0.
\tag{28}
\]

There are three independent Bianchi relations among the eight scalar
residual rows, consistent with five gauge-fixed field amplitudes. Before
gauge fixing, the spatial scalar metric has four amplitudes and the
matter field one. The Hamiltonian and two scalar-sector momentum
constraints, with their three coordinate freedoms, leave two physical
configuration degrees of freedom: a scalar-derived gravitational
polarization and the matter-scalar polarization. This is the generic
local canonical count; it is not a boundary characteristic-rank proof
or an assertion that a particular first-order family excites both.
The inherited master describes a restricted linear solution family and
cannot be presumed to span these general forced data.

A complete scalar contact construction still requires an evolution
reduction, a justified choice of independent constraint and boundary
rows, compatible second-order initial data, and a boundary estimate.
Equations (19), (24), and (27) are the concrete bulk inputs to that
calculation. The other scalar harmonics, especially the exceptional
monopole, remain distinct tasks.

## 8. Reproduction and evidence limits

`verify_scalar_bulk.py` has its own metric-jet, inverse-metric,
Christoffel, differentiated-Christoffel, Ricci, and scalar d'Alembertian
implementation. It imports only the Python standard library and SymPy.
It does not import an inherited checker, a serialized symbolic object,
or this note's radial ADM formulas.

The direct five-dimensional tests retain independent symbolic values
for every first and second \((y,t)\) derivative of the perturbations and
for \(\rho,\rho',\rho''\), \(p,p'\). They check all fifteen independent
Ricci components for both the quadratic first-order path and the general
linear RW path, plus both scalar d'Alembertian variations. Angular jets
are taken at the nonspecial rational stereographic point
\((1/3,1/2,2/5)\), with \(\cosh t=5/4\), \(\sinh t=3/4\).
The angular projection identities are separately verified globally as
rational functions. The all-point and general-\(\lambda\) proof of the
linear operator is the connection derivation (15); the direct local-jet
check alone is not described as an all-point computer proof.

Run, with a research-root path containing the four controlling inputs:

```sh
python verify_scalar_bulk.py --output normal.json --model-root RESEARCH_ROOT
python -O verify_scalar_bulk.py --output optimized.json --model-root RESEARCH_ROOT
python verify_scalar_bulk.py --inject-corruption
python -O verify_scalar_bulk.py --inject-corruption
```

The corruption switch changes the actual candidate \(yt\) quadratic
Ricci formula by adding \(u u_{yt}Y^2\). The defining Christoffel
calculation then rejects that equation through an explicit exception.
No check relies on Python `assert`, so optimization cannot remove it.
Receipts report the script SHA-256, model/source SHA-256 values, runtime
versions, the exact checks performed, and the stated local-jet scope.
These checks support the identities; they are neither numerical
contact evolutions nor existence proofs.

The frozen execution receipts contain **52 passing exact checks** in each of
ordinary and optimized Python. The actual corrupted formula exits with code
1 in both modes, with the direct quadratic Ricci `01` residual reported.
`VERIFICATION_RECEIPT.json` records these executions and the artifact hashes.
