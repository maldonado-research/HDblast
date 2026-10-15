# A complete finite-angular second-order coefficient through shell contact

14 October 2026. Conditional coefficient construction, internally reviewed.
This note concerns
a second Taylor coefficient on a controlled positive-cone collar. It does
not assert an exact nonlinear family through contact, an all-angular
nonlinear boundary theorem, external novelty, or a cosmological origin.

## 1. Statement, model and inputs

Use the registered five-dimensional canonical Einstein–scalar model,

\[
 \mathscr E_{AB}=R_{AB}-\phi_A\phi_B-\tfrac23U(\phi)g_{AB},
 \qquad \mathscr S=\Box_g\phi-U'(\phi).
\]

The background on each regular interior is
\(g_0=dy^2+\rho(y)^2[-dt^2+\cosh^2t\,\kappa_{S^3}]\),
\(\phi_0=\phi_0(y)\). Write \(a=\rho_y/\rho\), \(p=\phi_{0,y}\).
The outward normal points toward increasing \(y\). On two identical copies
the common physical shell laws reduce to

\[
 K_{ab}=\frac{\sigma(\phi)}6h_{ab},\qquad
 n\phi=-\frac{\sigma'(\phi)}2.
 \tag{1}
\]

Thus \(a_b=\sigma_0/6\), \(p_b=-\sigma'_0/2\). In the actual model
\(\sigma''=4\phi\), \(\sigma'''=4\). The bulk potential and shell
tension are the previously registered functions; no new physics is added.

Choose the copy-even aligned first master with
\(Y=n_1^2-n_2^2\), \(A=B=\operatorname{diag}(1,-1,0,0)\).
It solves the complete first bulk and shell equations in its inherited
longitudinal representative,

\[
 (h_1)_{yy}=-4\Psi,\quad(h_1)_{y\mu}=0,\quad
 (h_1)_{\mu\nu}=2\rho^2\Psi\gamma_{\mu\nu},
 \quad\Psi=u(y,t)Y,\quad\chi=v(y,t)Y,
 \quad v=-\frac{3(\rho^2u)_y}{\rho^2p}.
\]

The first shell displacement is zero. Initial master displacement and
velocity use the same nonzero smooth radial bump, with a fixed positive
velocity scale, supported in an annulus strictly inside the shell. This
commuting state passes all ten initial compatibility tests. Use the smooth
localized second-order initial constraint correction established by the
earlier localized checkpoint, including its stated localization/background
hypotheses. These inputs vanish near the initial shell corner.

Let \(x=\log(r/b)\), where \(dy=\Omega(r^2)dr\),
\(\rho=r\Omega(r^2)\). Then \(dx=dy/\rho\), shell \(x=0\), and
the cone is \(x=-\infty\). Let the first support outer edge be \(R_0<b\).
Choose a finite controlled interval

\[
 t_c=\log(b/R_0)<t_*<\operatorname{arsinh}(\tan\tau_*),
 \tag{2}
\]

where the inherited first-order contact solution exists on its conformal
slab \(|\tau|\le\tau_*\). A bump sufficiently close to the shell makes
this interval nonempty. Choose \(L\) containing all initial radial support.
The actual background is smooth on \([-L-t_*,0]\).

**Conditional coefficient theorem.** Under these inputs, there is a smooth
complete second coefficient \((k,w)\), with identical copies, satisfying

\[
 D\mathscr E[k,w]=-D^2\mathscr E[(h_1,\chi)^2],\qquad
 D\mathscr S[k,w]=-D^2\mathscr S[(h_1,\chi)^2],
 \tag{3}
\]

all induced-field matches and the second differentiated physical shell
laws, throughout the controlled collar and through \(t_c\). It can use
a finite-angular localized initial constraint correction and zero free
initial tensor coefficient. A boundary-preserving coordinate construction
below fixes a unique reduced evolution for chosen complete initial gauge
data. The physical coefficients can be returned to the first longitudinal
representative with all quadratic gauge terms retained.

The Taylor convention is
\(g_\epsilon=g_0+\epsilon h_1+\epsilon^2k/2\),
\(\phi_\epsilon=\phi_0+\epsilon\chi+\epsilon^2w/2\).
This polynomial satisfies the equations and shell laws through second
order. It is not an exact solution family for nonzero \(\epsilon\).

## 2. Why an unchanged scalar master is insufficient

The independently calculated \(l=4\) quadratic angular scalar-TF residual is

\[
 \left(8u^2-\tfrac23v^2\right)(\operatorname{Hess}F_4)^{\rm TF},
 \quad F_4=Y^2-\tfrac12(n_1^2+n_2^2)+\tfrac1{12}.
\]

For a general scalar RW coefficient
\(k_{ab}=f_{ab}F_4\), \(k_{ai}=0\),
\(k_{ij}=R^2K\kappa_{ij}F_4\), \(R=\rho\cosh t\), the linear
angular TF row is \(-\tfrac12(\operatorname{tr}_q f+K)
(\operatorname{Hess}F_4)^{\rm TF}\). Thus

\[
 \operatorname{tr}_q f+K=16u^2-\tfrac43v^2.
 \tag{4}
\]

The inherited conformal reconstruction has the left side zero identically.
For compact nonzero displacement, its right side cannot vanish everywhere:
the reconstruction and vanishing condition imply
\(|u_y+2au|=(2/\sqrt3)|pu|\), hence a bounded-coefficient
inequality \(|u_y|\le C|u|\). Gronwall from a radial exterior zero point
forces \(u=0\). Pure velocity data can have zero displacement on the initial
slice; that exception does not establish a scalar-only evolution later.
The independently derived temporal mixed row also shows that merely
repairing (4) is insufficient. The construction below keeps all fields.

## 3. Relative harmonic gauge and its first-order construction

Use the background connection as a fixed reference and define the tensorial
wave-map vector

\[
 C^A(g)=g^{BC}\bigl(\Gamma^A_{BC}[g]-\Gamma^A_{BC}[g_0]\bigr).
\]

Reduce the equations by
\(\mathscr E^{\rm red}_{AB}=\mathscr E_{AB}-\nabla_{(A}C_{B)}\),
with symmetrization containing the factor \(1/2\), and keep \(\mathscr S\).
Their linearization at \(g_0\) is normally hyperbolic on the fifteen metric
components and the scalar. The metric principal part is
\(-\tfrac12\Box_0 k_{AB}\); the scalar principal part is \(\Box_0w\).
Matter couplings and connection terms have at most one derivative of the
unknown coefficients.

The first longitudinal solution is generally not in this gauge. With lower
indices, \(DC_A[h]=\nabla_0^Bh_{AB}-\tfrac12\nabla_A^0\operatorname{tr}_0h\)
and \(DC[\mathcal L_\xi g_0]=(\Box_0+\operatorname{Ric}_0)\xi\).
For the longitudinal input,
\(C^{(1)}_y=-6(\Psi_y+4a\Psi)\), \(C^{(1)}_\mu=0\).
Solve the vector wave equation

\[
 (\Box_0+\operatorname{Ric}_0)\xi=-C^{(1)}_{\rm long},\qquad
 \xi_y|_b=0,\quad (\partial_y\xi_a-2a_b\xi_a)|_b=0,
 \tag{5}
\]

with zero initial \(\xi,\xi_t\). Section 6 applies to this finite-angular
vector system. Its source and data vanish near the initial shell corner.
Then \(h_1^w=h_1+\mathcal L_\xi g_0\),
\(\chi^w=\chi+\mathcal L_\xi\phi_0\) solve the first physical equations
and have \(C^{(1)}=0\). Equation (5) makes
\((h_1^w)_{ya}|_b=0\). The vector is tangent to the boundary, so its local
flow preserves the fixed shell and acts identically on the two copies.

Every quadratic source must now use \((h_1^w,\chi^w)\). In particular,
returning from the wave representative uses

\[
 k^{\rm long}=k^w-2\mathcal L_\xi h_1-\mathcal L_\xi^2g_0,
 \qquad w^{\rm long}=w^w-2\mathcal L_\xi\chi-\mathcal L_\xi^2\phi_0.
 \tag{6}
\]

An additional second-order generator, if chosen, has its own linear Lie
term. No untransformed longitudinal quadratic forcing is substituted in
the wave-gauge equations.

For an explicit reduced source, let \(Q_C=D^2C[(h_1^w)^2]\).
The actual second reduced metric equation is

\[
 D\mathscr E[k^w,w^w]-\nabla^0_{(A}DC[k^w]_{B)}
 =-D^2\mathscr E[(h_1^w,\chi^w)^2]
   +\nabla^0_{(A}(Q_C)_{B)}.
 \tag{6a}
\]

The scalar equation is \(D\mathscr S[k^w,w^w]
=-D^2\mathscr S[(h_1^w,\chi^w)^2]\). First-connection variations times
\(C^{(1)}\) vanish. Covariance also gives the exact source recipe
\(-D^2(\mathscr E,\mathscr S)[(q_1^w)^2]
=-D^2(\mathscr E,\mathscr S)[(q_1^{\rm long})^2]
+D(\mathscr E,\mathscr S)[B_g,B_\phi]\), with the completions in (6).
Thus the independently computed longitudinal bulk sources and the known
first gauge vector determine every physical wave-gauge source.

## 4. The exact primal boundary completion

At the shell impose only

\[
 g_{ya}=0\ (4\text{ rows}),\qquad (1)\ (11\text{ rows}),
 \qquad C_y=0\ (1\text{ row}).
 \tag{7}
\]

The four tangential harmonic constraints are not additional primal rows.
Write \(g_{yy}=N^2\), \(g_{ab}=h_{ab}\) at the boundary. Then
\(K_{ab}=h_{ab,y}/(2N)\). Since the reference connection has
\(\widehat\Gamma^y_{ab}=-a_bh_{0,ab}\), the exact conditions are

\[
 h_{ab,y}=\frac\sigma3Nh_{ab},\qquad
 \phi_y=-\frac{\sigma'}2N,\qquad
 N_y=\frac{2\sigma}3N^2-a_bN^3\operatorname{tr}_h h_0.
 \tag{8}
\]

This last row follows from
\(C^y=N_y/N^3-2\sigma/(3N)+a_b\operatorname{tr}_h h_0\).
There are four Dirichlet and twelve independent normal-derivative principal
rows. After using the first two equations in (8), the last contains no
tangential derivative or additional unknown normal derivative.

For explicit second-order sources write
\(N=1+\epsilon\nu_1+\epsilon^2\nu_2/2\),
\(h=h_0+\epsilon h_1^w+\epsilon^2k_{ab}/2\),
\(\phi=\phi_0+\epsilon\chi^w+\epsilon^2w/2\),
and \(A_1=h_0^{-1}h_1^w\). The differentiated rows are

\[
 k_{ab,y}-2a_bk_{ab}
 -(\sigma'_0w/3+2a_b\nu_2)h_{0,ab}=R_{ab},
\]
\[
 w_y+\sigma''_0w/2+\sigma'_0\nu_2/2=R_\phi,
\qquad
 \nu_{2,y}+4a_b\nu_2-a_b\operatorname{tr}_0 k
 -2\sigma'_0w/3=R_N,
 \tag{9}
\]
where

\[
 R_{ab}=(2\sigma'_0\chi^w/3+4a_b\nu_1)h^w_{1,ab}
 +(\sigma''_0(\chi^w)^2/3+2\sigma'_0\chi^w\nu_1/3)h_{0,ab},
\]
\[
 R_\phi=-2(\chi^w)^2-\sigma''_0\chi^w\nu_1,
\]
\[
 R_N=2\sigma''_0(\chi^w)^2/3+8\sigma'_0\chi^w\nu_1/3
 -16a_b\nu_1^2-2a_b\operatorname{tr}(A_1^2)
 +6a_b\nu_1\operatorname{tr}A_1.
 \tag{10}
\]

The metric wave unknown is \(k_{yy}=2\nu_2+2\nu_1^2\), not \(2\nu_2\).
Substituting this identity in (9) gives normal-derivative equations for
every metric/scalar wave unknown, with known quadratic right sides and
lower-order matrix couplings. The boundary sources are all fixed by the
complete transformed first fields. They vanish near the initial corner.

## 5. Finite angular initial constraint and gauge data

The first real quadrupole representation is \(V_2=(1,1)\) of SO(4),
dimension nine. Its symmetric square is

\[
 \operatorname{Sym}^2 V_2=(0,0)\oplus(1,1)\oplus(2,2)
 \oplus(2,0)\oplus(0,2),\qquad 1+9+25+5+5=45.
 \tag{11}
\]

These are scalar degrees zero, two, four and the two tensor degree-two
orientations. The gauge map (5), bulk operators and boundary operators
are SO(4)-equivariant linear maps on tensor/covector sections. Quadratic
sources are equivariant symmetric bilinear images of the same fixed
quadrupole. Consequently they remain in (11). This also proves the absence
of a forced coexact-vector representation for this selected state.

Take the compact-group isotypic projection \(\Pi\) onto (11), acting on
the full canonical tensor/density fields, not a componentwise scalar
Laplacian cutoff. The background linear constraint operator \(P\) commutes
with \(\Pi\). If \(Pq=S\) is the inherited localized correction, then
\(P(\Pi q)=\Pi S=S\). Rotations preserve radii, so this projection preserves
radial support and smoothness. It supplies finite-angular second physical
initial data. The initial tensor metric/momentum directions may be removed
as in the prior tensor note, without changing the sourced constraints.

Convert the physical initial data to coordinate velocities with the chosen
first lapse/shift. Choose second lapse and shift values that are smooth,
localized and in (11), vanishing near the initial shell corner (zero values
are admissible), and
solve algebraically for their five initial time derivatives using

\[
 C^{(2)}=DC[k]+Q_C=0,\qquad Q_C=D^2C[(h_1^w)^2].
 \tag{12}
\]

In coordinates with \(g_{0,tt}=-\rho^2\), \(g_{0,tI}=0\), the unknown
derivative coefficients are \(-1/(2\rho^2)\) for \(k_{tt,t}\) in
\(C_t^{(2)}\), and \(-1/\rho^2\) for \(k_{tI,t}\) in \(C_I^{(2)}\).
Everything else in these rows is already known from physical initial data,
initial metric choices and \(Q_C\). These coefficients are invertible.
This changes no spatial metric, extrinsic curvature or scalar Cauchy data.

The physical constraints are projections of
\(\mathscr H_{AB}=\mathscr E_{AB}-\tfrac12g_{AB}\operatorname{tr}_g\mathscr E\),
not \(\mathscr E_{nn}\) alone. At an initial orthonormal Cauchy frame,
when \(C^{(2)}=0\), a reduced solution has

\[
 \mathscr H^{(2)}_{00}=\tfrac12\partial_0C^{(2)}_0,
 \qquad \mathscr H^{(2)}_{0I}=\tfrac12\partial_0C^{(2)}_I.
 \tag{13}
\]

Thus the sourced physical initial constraints imply
\(\partial_0C^{(2)}=0\). All chosen gauge jets vanish near the initial
shell corner, since the physical corrections and first fields do there.
The reduced wave recursion therefore satisfies every smooth corner
identity of (9).

## 6. Constructive finite-mode wave lemma

After angular decomposition and multiplication by the nonzero principal
coefficients, each system used here has the form, in inward \(z=-x\ge0\),

\[
 U_{tt}-U_{zz}=A(t,z)U_t+B(t,z)U_z+D(t,z)U+F(t,z),
 \tag{14}
\]

with finitely many components and smooth bounded coefficients on the
controlled support. Its boundary is Dirichlet for a subset and Robin for
the rest:

\[
 U_D(0,t)=d(t),\qquad
 -U_{N,z}(0,t)+M(t)U_N(0,t)=g(t).
 \tag{15}
\]

Couplings to prescribed Dirichlet traces can be absorbed in \(g\).
No symmetry or favorable sign of the lower-order matrix \(M\) is assumed.
Take smooth compatible initial data, all sources zero near the initial
corner, and radial support \(z\le L\). This lemma is a finite-mode
construction; it makes no estimate uniform in all angular frequencies.

Define \(p=U_t+U_z\), \(m=U_t-U_z\). Their source is
\(\mathcal F=(A+B)p/2+(A-B)m/2+DU+F\), and

\[
 p_t-p_z=\mathcal F,\qquad m_t+m_z=\mathcal F.
\]

The outgoing value is
\(p_{o}(t)=p_0(t)+\int_0^t\mathcal F(s,t-s)ds\).
For the Robin trace \(q=U_N(0,t)\), (15) gives

\[
 q'=p_{o,N}-Mq+g,\qquad
 m_N(0,t)=p_{o,N}-2Mq+2g.
 \tag{16}
\]

Solve the finite-dimensional trace ODE using its fundamental matrix:
\(q(t)=\Phi(t,0)q(0)+\int_0^t\Phi(t,s)[p_{o,N}(s)+g(s)]ds\).
For Dirichlet components use \(m_D(0,t)=2d'(t)-p_{o,D}(t)\).
Insert the outgoing characteristic formula into this ODE before forming
the incoming characteristic formula. For \(z\ge t\) the latter starts on
the initial slice; for \(z<t\) it starts at the boundary time \(t-z\).
Finally \(U(t,z)=U_0(z)+\tfrac12\int_0^t(p+m)(s,z)ds\).

This defines a Volterra map on \((p,m,U)\). On an interval of length
\(\delta\), the outgoing difference is bounded by
\(C\delta\|(\Delta p,\Delta m,\Delta U)\|_\infty\).
The solved trace difference is bounded by \(C\delta\) times that outgoing
difference. Thus the incoming boundary difference, including the term
\(M\Delta q\), and every characteristic integral are bounded by
\(C\delta\|(\Delta p,\Delta m,\Delta U)\|_\infty\).
The final time integral has the same bound. Choose \(C\delta<1\).
Contraction gives a unique continuous fixed point. The identities imply
\(U_z=(p-m)/2\) and the prescribed boundary trace: their defects have
zero initial value and zero evolution. Differentiation and the compatible
corner recursion give a smooth solution. Successive bounded intervals
cover any specified finite \(t_*\); initial data at each step are the
compatible jets of the preceding solution.

This argument avoids the unsupported energy cancellation of a nonsymmetric
Robin matrix. It is not a gravitational energy, ghost or stability claim.
The characteristic formulas also prove finite propagation: if the source
vanishes for \(z>L+t\), the solution does there. Extending the bounded
background coefficients to the left of the physical support is only a
construction device and cannot affect the resulting coefficient.

## 7. Subsidiary constraints and completion of the physical equations

The exact trace-reversed identity is

\[
 \nabla^A\mathscr E_{AB}-\tfrac12\nabla_B\operatorname{tr}_g\mathscr E
 =-\mathscr S\,\phi_B.
 \tag{17}
\]

For a reduced solution with \(\mathscr S=0\), it gives
\((\Box_g+\operatorname{Ric}_g)C=0\). At second order,
\(C^{(0)}=C^{(1)}=0\) remove all lower-order inhomogeneous terms, hence
\((\Box_0+\operatorname{Ric}_0)C^{(2)}=0\).

The physical boundary equations induce the tangential conditions. Codazzi
in our outward conventions gives

\[
 R_{na}=D_bK^b{}_a-D_aK=-3D_a(\sigma/6)
 =-\tfrac12\sigma' D_a\phi=(n\phi)D_a\phi.
 \tag{18}
\]

Thus \(\mathscr E_{na}=0\), including all quadratic matter and evaluation
terms. On a reduced solution this equals
\(\nabla_{(n}C_{a)}\). Since \(C_n=0\) at the boundary, the projected
subsidiary law is

\[
 (\nabla_nC)_a-K_a{}^bC_b=0.
\]

In boundary-orthogonal coordinates it is
\(\partial_yC_a-2NK_a{}^bC_b=0\). Its second derivative is

\[
 C_y^{(2)}|_b=0,\qquad
 (\partial_yC_a^{(2)}-2a_bC_a^{(2)})|_b=0,
 \tag{19}
\]

because all first subsidiary fields vanish. Section 6 gives uniqueness
for this finite-angular vector system. Its zero initial value and zero
time derivative from (12)–(13) imply \(C^{(2)}=0\). Consequently the
reduced coefficient solves every physical Einstein row, the scalar row,
the sourced initial constraints and their propagation through contact.
The normal Einstein projection follows from these full equations; no extra
primal boundary condition or arbitrary incoming constraint row is appended.
The inherited exact normal-force identity expresses the force residual in
terms of the configuration match, Israel/scalar residuals and the difference
of the two normal bulk residuals. They all vanish through second order here.
The identical-copy construction also makes the differences of outward
extrinsic traces and scalar normal derivatives zero. Thus the second shell
normal-force row follows, rather than supplying an additional primal row.

## 8. Construction recipe, domain and limits

1. Evolve the known aligned first master on its controlled contact slab.
2. Solve (5); form the complete wave-gauge first fields.
3. Project a smooth localized second physical constraint correction onto
   (11), choose free tensor data as stated, and complete initial gauge
   derivatives using (12).
4. Form the full reduced quadratic source by differentiating the displayed
   reduced equations, including \(Q_C\), and form (9)–(10).
5. Expand in the full finite tensor/covector representation space (11),
   without division by the exceptional \(l=0\) eigenvalue. Solve the
   resulting matrix systems by the convergent characteristic construction.
6. Use (17)–(19) to establish all physical equations, and (6) to return to
   the longitudinal first representative if desired. Use identical copies.

This is a constructive specification of the whole coefficient, not only
its angular TT projection. It uses the independently checked scalar sources
and boundary algebra as checks; it does not require an unproved single-master
or two-decoupled-master ansatz. It provides no all-angular uniform estimate,
nonlinear contact existence, convergence of an all-order formal series or
global stability result. The domain stays in the actual positive-cone
collar with support away from the cone at finite global de Sitter time.

The boundary completion follows the mechanism of Fournodavlos–Smulevici,
*The initial boundary value problem in General Relativity: the umbilic case*,
[arXiv:2104.08851v1](https://arxiv.org/abs/2104.08851v1), equations (16)–(18)
and Proposition 2. Their theorem is four-dimensional vacuum with constant
umbilicity. Equations (8), (17)–(19) are derived here for the actual
five-dimensional scalar-dependent shell. This note does not transfer their
theorem or claim that the extension is externally novel. Reviews and
symbolic verification are internal AI-assisted work, not external peer review.
