# Full packet junction sufficiency and the actual shell current

14 October 2026 UTC. Independent local calculation in an authoring directory.
The frozen research capsules are unchanged. No existence theorem for a weaker
completed space, physical completeness theorem, novelty claim, nonlinear
stability result, or cosmological origin claim is made here.

## 1. Exact statement and hypotheses

Fix an actual analytic regular Einstein–scalar background, with

\[
 ds^2=dy^2+\rho(y)^2\gamma_{\mu\nu}dx^\mu dx^\nu,
 \quad a=\rho'/\rho,\quad p=\phi',\quad H=\rho^{-2},
 \quad R_{\mu\nu}(\gamma)=3\gamma_{\mu\nu}.
\]

At its doubled timelike shell \(y=y_b\), both interior normals point toward
increasing \(y\). The background equations and shell junctions are

\[
 a'=-H-p^2/3,\quad p'=U_\phi-4ap,
 \quad\sigma_b=6a_b,\quad\sigma'_b=-2p_b.
\]

Assume \(a_b\ne0,p_b\ne0,B_b=p'_b/p_b+\sigma''_b/2\ne0\).
The registered branch has \(a_b>0,B_b>0\); positivity is not needed for the
algebraic sufficiency proof. All shell coefficients below are background
constants. The bulk reconstruction hypothesis is the exact proved map in
the frozen regular-cone capsule:

\[
 P=r^4F,\quad \Psi=P/\rho^2,\quad
 \chi=-3P'/(\rho^2p),\quad
 h_{yy}=-4\Psi,\quad h_{y\mu}=0,\quad
 h_{\mu\nu}=2\rho^2\Psi\gamma_{\mu\nu}.
\]

Here \(r'=r/\rho\), \(s=r^2=X\cdot X\), \(\rho=r\Omega(s)\), and

\[
 \mathcal LF=\Box_\eta F-2L(s)EF-[8L+(8/3)s\Phi_s^2]F=0,
 \qquad L=\Omega_s/\Omega+2\Phi_{ss}/\Phi_s.
\]

For every smooth solution of this master equation, the complete bulk
Einstein and scalar residuals vanish. In the shell collar it also gives,
as identities for arbitrary smooth functions and without harmonic division,

\[
 \Psi'+2a\Psi=-p\chi/3,
 \qquad p\chi'=p'\chi+3H(\Box_\gamma+4)\Psi-2p^2\Psi.       \tag{1}
\]

Let \(F_\pm=F_e\pm F_o\). Suppose \(F_e,F_o\) are smooth up to the
shell and satisfy respectively

\[
 nF_e-B_b^{-1}\Box_{h_b}F_e+
 [4/\rho_b-4/(\rho_b^2B_b)]F_e=0,
\]
\[
 nF_o+[4/\rho_b+p_b^2/(3a_b)]F_o=0.                       \tag{2}
\]

**Sufficiency theorem.** Set the normal graph displacements to

\[
 \zeta_e=0,\qquad
 \zeta_o=-\Psi_{ob}/a_b=-P_{ob}/(a_b\rho_b^2),
 \qquad\zeta_\pm=\pm\zeta_o.                             \tag{3}
\]

Then the reconstructed bulk fields, common induced metric/scalar, and these
two embeddings satisfy every linearized metric and scalar junction on the
shell. This statement applies to every smooth master solution satisfying
(2), including arbitrary smooth packets and the exceptional constant and
Hessian-degenerate harmonic subspaces. There is no additional harmonic
constraint on this sufficient solution family. Existence of such master
solutions and initial/boundary compatibility are separate PDE questions.

## 2. Direct graph geometry

Use the common shell coordinates \(x^\mu\), and write each graph as
\(y_i=y_b+\zeta_i(x)\). No tangential reidentification of the copies is
introduced; the choice of common chart is a sufficient matching convention.
To first order the normal covector and graph tangents are

\[
 n_A=(1-2\Psi,-D_\mu\zeta),\qquad
 e_\mu^A=(D_\mu\zeta,\delta_\mu^\nu).
\]

Pullback, normalization, and evaluation at the displaced graph give

\[
 \delta h^{\rm ind}_{\mu\nu}
   =2\rho_b^2(\Psi_b+a_b\zeta)\gamma_{\mu\nu},
 \qquad \delta\phi^{\rm ind}=\chi_b+p_b\zeta,              \tag{4}
\]
\[
 \delta K^\mu{}_{\nu}
   =(\Psi'+2a\Psi+a'\zeta)\delta^\mu{}_{\nu}
       -H D^\mu D_\nu\zeta,
 \qquad \delta p_n=\chi'+2p\Psi+p'\zeta.                 \tag{5}
\]

For clarity, before raising the index, the graph curvature is

\[
 \delta K_{\mu\nu}
 =\rho_b^2[\Psi'+4a\Psi+(a'+2a^2)\zeta]\gamma_{\mu\nu}
       -D_\mu D_\nu\zeta.
\]

The \(a'\zeta\) term in (5) comes from evaluating the background curvature
on the graph; omitting it changes the junction. The normal scalar variation
includes both the perturbed lapse \(2p\Psi\) and the evaluation shift
\(p'\zeta\). These are direct geometry identities for arbitrary spacetime
functions, not equations extracted from a single harmonic.

## 3. All matching and junction residuals

Write every field/displacement as an even plus an odd part. Continuity of
(4) is exactly

\[
 \Psi_o+a\zeta_o=0,\qquad\chi_o+p\zeta_o=0.               \tag{6}
\]

Thus \(\zeta_o=-\Psi_o/a\), and scalar continuity is equivalent to

\[
 P_o'=-\frac{p^2}{3a}P_o.                                \tag{7}
\]

Equation (7) is precisely the odd law (2), since
\(P'=r^4(nF+4F/\rho)\). It involves no derivative or eigenvalue division
in a tangential direction. The common induced perturbations are
\(2\rho_b^2\Psi_e\gamma_{\mu\nu}\) and \(\chi_e\) with (3).

The original mixed Israel residual is

\[
 \mathcal I^\mu{}_{\nu}
 =\sum_i(K_i^\mu{}_{\nu}-K_i\delta^\mu{}_{\nu})
       +\sigma\delta^\mu{}_{\nu}.
\]

Its trace-reversed form
\(\mathcal R^\mu{}_{\nu}=\sum_iK_i^\mu{}_{\nu}
 -\sigma\delta^\mu{}_{\nu}/3\) vanishes if and only if
\(\mathcal I\) vanishes; the trace transformation is invertible in four
dimensions. Inserting (4)–(5), continuity, and \(\sigma'=-2p\) gives

\[
 \delta\mathcal R^\mu{}_{\nu}
 =2(\Psi_e'+2a\Psi_e+p\chi_e/3)\delta^\mu{}_{\nu}
       -2H(D^\mu D_\nu+\delta^\mu{}_{\nu})\zeta_e.
\]

The first term vanishes by (1), hence

\[
 \boxed{\delta\mathcal R^\mu{}_{\nu}
       =-2H(D^\mu D_\nu+\delta^\mu{}_{\nu})\zeta_e.}       \tag{8}
\]

This is the complete tensor residual, including trace and trace-free
equations. The choice \(\zeta_e=0\) makes every component zero pointwise.
The summed scalar junction is independently

\[
 \mathcal S_\Sigma=\sum_i n_i\cdot\nabla\phi_i+\sigma'(\phi_{\rm ind}),
\]
\[
 \delta\mathcal S_\Sigma
 =2[\chi_e'+2p\Psi_e+\sigma''\chi_e/2+pB\zeta_e].
\]

Using the second identity (1) gives

\[
 \boxed{\delta\mathcal S_\Sigma
 =-\frac6{\rho_b^2p_b}
     [B_bP_e'-\rho_b^{-2}(\Box_\gamma+4)P_e]
       +2p_bB_b\zeta_e.}                                \tag{9}
\]

With (3), equation (9) vanishes exactly by the even law (2). Equations
(6), (8), and (9) exhaust the linearized matching/junction conditions in
this reconstructed scalar family. The derivation has divided only by
nonzero background coefficients \(a_b,p_b,B_b,\rho_b\).

There is an important distinction between sufficiency and necessity.
For a fully general reconstructed graph, (8) only requires

\[
 (D_\mu D_\nu+\gamma_{\mu\nu})\zeta_e=0.                  \tag{10}
\]

On global unit de Sitter space this has the five ambient-coordinate
solutions \(\zeta_e=c_A n^A(x)\): the value and gradient at one point
determine the solution, so its dimension is at most five, and these five
independent solutions span it. All satisfy
\((\Box_\gamma+4)\zeta_e=0\). In that kernel the scalar condition (9)
retains \(\zeta_e\); the straight even Wentzell law alone need not describe
every allowed embedding representative. This is a finite-dimensional
exceptional-completeness issue, not a missing constraint on solutions of
(2)–(3). Constant displacement has nonzero tensor in (10) and is excluded
unless it is zero. Merely dividing the trace-free Hessian equation would
miss the trace condition. Neither an entire separated mass value nor all
data with that mass value should be removed from the sufficient theorem.

## 4. Normal and tangential consistency do not add a constraint

With common induced fields, the difference of the two bulk Hamiltonian
constraints gives

\[
 3a_b(\delta K_+-\delta K_-)
       =p_b(\delta p_{n,+}-\delta p_{n,-}).                \tag{11}
\]

Indeed each constraint is
\(\tfrac12(K^2-K_{\mu\nu}K^{\mu\nu}-R[h])
 =\tfrac12p_n^2-U\) at first order; its induced-curvature and scalar-value
variations cancel in the difference. Thus the normal shell force equation
is a consequence of the bulk constraints and the junctions. One can check
it directly for the odd graph using (1), (6):

\[
 \delta K_o=-H_b(\Box_\gamma+4)\zeta_o,
 \qquad\delta p_{n,o}
       =-\frac{3a_bH_b}{p_b}(\Box_\gamma+4)\zeta_o.
\]

Consequently
\(\sigma_b\delta K_o+\sigma'_b\delta p_{n,o}
 =6a_b\delta K_o-2p_b\delta p_{n,o}=0\) identically.
There is no extra odd boundary wave equation to impose on (7).

Tangential shell conservation follows from the bulk Codazzi equations:
the common scalar trace \(C\) obeys
\(\sum_i R_{n_i\mu}=2p_bD_\mu C=-\sigma'_bD_\mu C\),
which equals the divergence of Israel with shell stress
\(-\sigma h_{\mu\nu}\). At the level of (8),
\(D_\nu[(D^\nu D_\mu+\delta^\nu{}_{\mu})\zeta_e]
 =D_\mu(\Box_\gamma+4)\zeta_e\); this uses
\(R_{\mu\nu}(\gamma)=3\gamma_{\mu\nu}\) and shows the divergence of the
original Israel residual vanishes. Bulk reconstruction already supplies
the Einstein momentum and scalar constraints even in exceptional harmonic
subspaces. They were not obtained here by dividing those tensors.

## 5. Actual Einstein–scalar packet current and the moving shell corner

Use the exact action and potential convention of the frozen covariant
capsule:

\[
 S=\sum_{i=\pm}\left[\int_{M_i}\epsilon_g
       (R/2-(\nabla\phi)^2/2-U)+\int_\Sigma\epsilon_hK_i\right]
       -\int_\Sigma\epsilon_h\sigma(\phi),
\]
\[
 \theta^A=\sqrt{-g}\{\tfrac12(\nabla_Bh^{AB}-\nabla^Ah)
                  -\nabla^A\phi\,\chi\},
 \quad\omega^A=\delta_1\theta^A(2)-\delta_2\theta^A(1).
\]

The reconstruction satisfies (1) for all smooth packets. The direct
field-space curl on the conformal zero-shift slice therefore yields,
per interior,

\[
 \omega^y=\sqrt{-\gamma}\rho^4
 [\chi_1\chi_2'-\chi_2\chi_1'-2p(\Psi_1\chi_2-\Psi_2\chi_1)]
 =9\sqrt{-\gamma}\alpha
 [((\Box_\gamma+4)P_1)P_2'-((\Box_\gamma+4)P_2)P_1'],      \tag{12}
\]

\[
 \omega^\mu=\sqrt{-\gamma}\rho^2
 [\chi_1D^\mu\chi_2-\chi_2D^\mu\chi_1
       +6(\Psi_1D^\mu\Psi_2-\Psi_2D^\mu\Psi_1)],
 \qquad\alpha=(\rho^2p^2)^{-1}.                          \tag{13}
\]

These formulas refer to the actual bulk covariant current, not the
auxiliary regular-master energy. Their application to this family does
not inherit a generic-harmonic exclusion: the required longitudinal slice
and identities (1) have been directly proved for every reconstructed F.

For two even solutions the boundary law
\(P_b'=H_b(\Box_\gamma+4)P_b/B_b\) makes (12) zero at the shell
**pointwise**, regardless of angular/time dependence or separated masses.
There is no extra straight-even-shell tangential current on this slice.
The Wentzell coefficient is not an additional independent physical shell
kinetic field.

For two odd solutions put \(c_o=(3a_b\rho_b^2)^{-1}\).
The complete graph/shell quadratic derivative action is

\[
 S^{(2)}_{\Sigma,\mathrm{derivative}}
       =-3a_b\rho_b^2\int\sqrt{-\gamma}(D\zeta_o)^2
       =-9c_o\int\sqrt{-\gamma}(DP_{ob})^2.              \tag{14}
\]

Here both gravitational copies and one common shell are included. The
accounting can be derived directly, without separated harmonics. For
one graph, its GHY density minus the slice GHY, after tangential integration
by parts, contributes
\(\rho_b^2[3a_b(D\zeta_i)^2-D\zeta_i\cdot DN_i]\).
The Einstein tangential-divergence graph corner adds
\(+\rho_b^2D\zeta_i\cdot DN_i\). The shifted upper radial integration
domain adds \(+3\rho_b^2D\zeta_i\cdot D\Psi_i\).
The single common shell volume adds
\(-\sigma_b\rho_b^2(D\zeta_o)^2/2\).
Adding the two copies and imposing \(\Psi_o=-a_b\zeta_o\),
\(\sigma_b=6a_b\) gives (14). Omitting the moving-domain term instead
produces the opposite sign. Mixed even/odd domain terms cancel between
copies; the action is invariant under exchanging the interiors.

One exact way to check the graph term is to use a temporal test function
\(Q(t)\) in a flat de Sitter chart, \(\sqrt{-\gamma}=e^{3t}\).
Let \(n=e^{\epsilon N(y)Q(t)}\),
\(H_t=\rho^2e^{2\epsilon\psi(y)Q(t)}\),
\(V=\rho^4e^{4\epsilon\psi Q+3t}\), and
\(f=y_b+\epsilon zQ\). The normalized graph has
\(d=(1-n^2 f_t^2/H_t)^{1/2}\), and exactly

\[
 \sqrt{-h_{\rm graph}}K_{\rm graph}
 =\frac d n\left[\partial_y(V/d)
        +\partial_t(n^2Vf_t/(H_t d))\right]_{y=f}.
\]

Its second variation minus slice GHY includes the graph evaluation shift.
Integration of \(Q(Q_{tt}+3Q_t)\) gives \(-Q_t^2\), yielding the above
coefficient \(\rho_b^2(3a_bz^2-N_bz)(DQ)^2\).
The tensor/scalar contractions determine the covariant derivative
coefficients, so this test does not impose a homogeneous on-shell harmonic.

Polarizing (14) fixes the completed shell current:

\[
 \boxed{\omega_\Sigma^\mu
 =18c_o\sqrt{-\gamma}
   [P_{1b}D^\mu P_{2b}-P_{2b}D^\mu P_{1b}]
 =6a_b\sqrt{-h_b}
   [\zeta_1D_{h_b}^\mu\zeta_2-\zeta_2D_{h_b}^\mu\zeta_1].} \tag{15}
\]

The odd Robin law inserted in (12) gives, without a common spectral value,

\[
 \boxed{\partial_\mu\omega_\Sigma^\mu
       =\omega_+^y(y_b)+\omega_-^y(y_b).}                \tag{16}
\]

The zeroth-order 4 terms cancel from the antisymmetric divergence. A mixed
even/odd bilinear has zero summed bulk current and zero shell current. For
general parity superpositions the even and odd balances add. Tangential
integrations by parts in deriving the action are justified by compact test
variations, or by retaining cut terms. Actual nonzero solutions are not
assumed compact in time. Equation (15) is precisely the retained moving
shell cut contribution; it is not discarded by that integration.

For an ambient-time slab, write \(T=r_b\sinh t\),
\(R_b(T)=\sqrt{r_b^2+T^2}\). If \(\omega^T,\omega^R\) are bulk vector
densities in polar ambient coordinates, the normal boundary flux is

\[
 \omega^R-R_b'(T)\omega^T=\omega^y/R_b(T).
\]

The shell density transforms as \(\omega_\Sigma^T=\omega_\Sigma^t\);
its angular densities acquire \(dt/dT=1/R_b\). Thus (16) is the exact
moving-boundary balance on finite \(T\) slabs. For a smooth full interior
with no additional physical boundary, angular integration gives

\[
 \frac d{dT}\left[\sum_i\int_{R\le R_b(T)}\omega_i^T\,dR\,d\Omega
                  +\int_{S^3}\omega_\Sigma^T\,d\Omega\right]=0. \tag{17}
\]

One may multiply the entire expression by the sign chosen for future
Cauchy orientation. The cone is an internal smooth hypersurface when its
ambient extension is retained. If a subregion ends at the cone, its known
cone flux must be retained as an additional boundary term. No physical
norm, positivity, or completed phase-space identification follows solely
from this antisymmetric conservation law.

## 6. A clean non-gauge witness for the after-contact existence family

For collar-supported data with a finite sum of \(S^3\) scalar harmonics
\(H_l\), \(l\ge2\), there is a direct gauge obstruction independent of
the small-ball symplectic positivity argument in the frozen capsule.
In polar ambient coordinates write the background as

\[
 g=g_{ab}(z)dz^adz^b+w(z)^2q_{ij}d\vartheta^id\vartheta^j,
 \quad z=(T,R),\quad w=\Omega(R^2-T^2)R.
\]

The reconstructed metric has \(h_{ai}=0\) and its angular tensor is pure
trace. Project a putative smooth gauge covector onto one scalar harmonic:
\(\xi_a=f_a(z)H_l\), \(\xi_i=J(z)D_iH_l\). Its Lie derivative has

\[
 (\mathcal L_\xi g)_{ij}^{\rm TF}
       =2J(D_iD_jH_l)^{\rm TF},
 \qquad
 (\mathcal L_\xi g)_{ai}
       =[\partial_aJ-2(w_a/w)J+f_a]D_iH_l.              \tag{18}
\]

On unit \(S^3\), with \(\lambda=l(l+2)\), the integrated Bochner identity
gives

\[
 \int|\operatorname{Hess}^{\rm TF}H_l|^2
  =\tfrac23\lambda(\lambda-3)\int|H_l|^2
  =\tfrac23l(l+2)(l-1)(l+3)\int|H_l|^2>0.              \tag{19}
\]

Therefore zero angular trace-free perturbation forces \(J=0\), and zero
mixed perturbation forces \(f_a=0\). Scalar, divergence-free vector, and
tensor projections are orthogonal on the compact sphere; other smooth
diffeomorphism sectors cannot supply the missing scalar projection.
Thus no nonzero reconstructed \(l\ge2\) state is pure gauge on an open
polar annulus. Any global smooth pure gauge state would restrict to such a
local pure gauge state, so this also excludes it globally wherever the
nonzero collar solution exists. The shell embedding variation must still
be included for any charge computation; this metric obstruction already
rules out a pure gauge pair independently of that issue. The proof does
not classify \(l=0,1\), finite-cut charges, or weaker nonsmooth gauge groups.

## 7. Nonzero physical pairing and a thin-annulus data image

The exact Cauchy pullback from the frozen capsule extends algebraically to
any compact positive-s annulus. Its earlier small-ball restriction was used
to bound its zeroth-order coefficient, not to derive the identity. At
ambient \(T=0\), write smooth master data as \((q,v)=(F,F_T)\), with compact
support away from the shell, and put

\[
 D_r=r\partial_r,\quad m_r=\frac{9r^3}{4\Omega f^2}>0,
 \qquad f=\Phi_s.
\]

For one copy and future Cauchy orientation the actual pairing is

\[
 \sigma((q_1,v_1),(q_2,v_2))=K(q_1,v_2)-K(q_2,v_1),
\]
\[
 K(q,v)=\int dr\,d\Omega_3
 \left[m_r(D_r+4)q(D_r+5)v+\frac{6r^7}{\Omega}qv\right].  \tag{20}
\]

This follows directly from (13): at \(T=0\),
\(\chi=-3(D_r+4)q/(2\Omega^2f)\),
\(\chi_T=-3(D_r+5)v/(2\Omega^2f)\),
\(\Psi=r^2q/\Omega^2\), \(\Psi_T=r^2v/\Omega^2\), and the transformed
physical measure is \(\Omega^3r^3dr\,d\Omega_3\).
The extra unit in \(D_r+5\) comes from
\([\partial_T,E]=\partial_T\). Using
\(r\partial_r\log m_r=3-2sL\), compact-support integration by parts gives

\[
 K(q,q)=\int dr\,d\Omega_3
 \left[m_r\{(D_rq)^2+(2+9sL)q^2\}
                 +\frac{6r^7}{\Omega}q^2\right].        \tag{21}
\]

Choose an annulus \(I=(r_1,r_2)\) strictly inside \(r_b\), with length
\(w_I=r_2-r_1\), and define

\[
 m_* =\min_{\overline I}(m_rr^2)>0,\qquad
 M_* =\max_{\overline I}m_r(-2-9sL)_+<\infty.
\]

For smooth q compactly supported in I, Dirichlet Poincaré on I and the
nonnegative last term in (21) give

\[
 \boxed{K(q,q)\ge
 [\pi^2m_*/w_I^2-M_*]\int_I dr\,d\Omega_3\,q^2.}         \tag{22}
\]

Thus any fixed sufficiently thin annulus satisfies a strictly positive
bracket. Such an annulus exists in every regular shell collar: \(m_*\)
stays bounded away from zero and \(M_*\) stays bounded while its length
tends to zero. This argument uses only radial Poincaré and works for a
fixed finite set of angular harmonics, including the \(l\ge2\) family
in Section 6. Smooth compact data vanish on an open shell neighborhood
initially, so the moving-shell corner initially vanishes.

Same-parity initial pairs \((q,0)\) and \((0,q)\) therefore have doubled
actual physical pairing \(2K(q,q)>0\). Whenever smooth master solutions
of (2) exist for those data, (17) conserves this nonzero total pairing
through shell contact, including the odd moving-shell cut term. This does
not establish a nonzero *local shell* current bilinear: the total pairing
includes the bulk interior, and nonzero trace alone does not fix a local
Wronskian.

On the image of all smooth \((q,v)\) supported in this fixed thin annulus
and with a fixed finite harmonic set, the physical pairing is weakly
nondegenerate. If \((q,v)\) pairs to zero against all such data, testing
first against \((0,q)\) gives \(K(q,q)=0\), hence \(q=0\); testing against
\((v,0)\) then gives \(K(v,v)=0\), hence \(v=0\). Smooth uniqueness and
(17) transport this statement to its solution image on the available
slab. There is no claim about an all-angular completed space or a full
gauge quotient. K need not be symmetric; the antisymmetrized expression
(20) is the physical symplectic form, and positivity of its diagonal
auxiliary bilinear is not a physical energy theorem.

## 8. Verification and limitations

`check_packet_junctions.py` is independently written and imports no scientific
code from the frozen inputs. It calculates all ten covariant and all ten
mixed graph-curvature components, induced matching, and the scalar normal
variation from a direct Christoffel calculation with arbitrary dependence
on the four de Sitter coordinates. It separately checks the two-sided
residuals, normal-force identity, the field-space curl of the completed
radial potential, arbitrary-packet current balances, the exact graph-GHY
second variation including the graph shift, the full moving-shell
coefficient, the compact-annulus physical Cauchy identity, and the
\(S^3\) Hessian identity. Wrong controls detect missing
scalar lapse, an incorrect odd Robin sign, omitted graph evaluation,
omitted gravitational mixing, an invalid common-spectral-value shortcut,
missing moving-domain and Einstein-corner terms, omitted odd current, and
double-counted shell tension. Normal and optimized runs plus an injected
failure are recorded separately. Check counts are algebraic verifications,
not a count of independent physical results.

The proof closes the full linearized-junction sufficiency gap **conditional
on existence of smooth master solutions obeying (2)**. It supplies their
actual shell current balance and a non-gauge witness for finite angular
collar families with \(l\ge2\). It does not prove the unrestricted smooth
master initial-boundary theorem, compatibility in a finite-energy
completion, all perturbation-sector completeness, an admissible global
gauge quotient, or nonlinear dynamics. Prior-art originality remains
unassessed. There is no higher-dimensional Big Bang conclusion.
