# Independent scalar operator and exceptional-harmonic audit

12 October 2026 Pacific; calculation stamp 2026-10-13 UTC. This is a new
scratch audit. It changes no frozen capsule, manuscript, APS material, or
repository file. The script and its JSON receipt accompany this note.

The generic, reflection-even scalar equations and their displaced-shell
boundary condition agree with the archived operator. The exceptional
**harmonic subspaces** can also be handled directly: the static regular
deformation is excluded by the nonzero static junction determinant, while a
relative translation-type displacement is excluded by `p B != 0`. These are
distinct arguments; simply substituting the exceptional mass values into a
formula obtained by dividing by a harmonic tensor would be invalid.

The positive Sturm–Liouville norm is not itself the physical kinetic norm.
A direct quadratic kinetic calculation below gives, for generic
reflection-even separated modes on the admissible longitudinal gauge slice,

\[
 K_{\rm kin}=2\int_0^T\rho^2(|X|^2+6|\Psi|^2)dy
            =18\mathcal E[P]
            =18\Lambda\mathcal N_{SL}[P].                 \tag{A}
\]

The last equality is only at a spectral root with its endpoint term removed
correctly. The first equality is obtained from Einstein–scalar plus
Gibbons–Hawking–York (GHY) and shell terms, with a straight shell in this
sector. It is not an action guessed from the radial ordinary differential
equation. A global symplectic-domain/gauge-completeness theorem, and the
moving-shell kinetic calculation for the reflection-odd sector if that sector
is admitted, are **not** supplied by this local calculation. No conclusion
about the entire theory being free of ghosts follows from positivity of
\(\mathcal N_{SL}\) alone.

## 1. Model, variational boundary data, and scope

Take unit-curvature \(dS_4\), \(R_{\mu\nu}(\gamma)=3\gamma_{\mu\nu}\), with
signature \((-+++)\), and

\[
 ds^2=dy^2+\rho^2\gamma_{\mu\nu}dx^\mu dx^\nu,
 \quad p=\phi',\quad a=\rho'/\rho,\quad \mathcal H=\rho^{-2}.
\]

The metric action must include its variational boundary completion:

\[
 S=\sum_{i=\pm}\left[\int_{M_i}\!\sqrt{-g}
  \left(\frac R2-\frac12(\nabla\phi)^2-U\right)
  +\int_{\Sigma}\!\sqrt{-h}\,K_i\right]
  -\int_{\Sigma}\!\sqrt{-h}\,\sigma(\phi).
                                                        \tag{1}
\]

Both normals at the shell point toward increasing interior radius. Varying
the common inverse induced metric and scalar gives

\[
 \sum_i(K^{(i)}_{\mu\nu}-K_i h_{\mu\nu})
                  +\sigma h_{\mu\nu}=0,
 \qquad \sum_i n_i\cdot\nabla\phi+\sigma'=0.
\]

For identical copies these reduce to
\(K^\mu{}_{\nu}=\sigma\delta^\mu{}_{\nu}/6\) and
\(p=-\sigma'/2\). Thus \(a=\sigma/6\). The bulk variation gives
\(R_{AB}=\partial_A\phi\partial_B\phi+(2/3)Ug_{AB}\),
\(\Box\phi-U'=0\), and hence

\[
 p'=U'-4ap,\qquad a'=-\mathcal H-p^2/3,
 \qquad a^2=\mathcal H+p^2/12-U/6.                         \tag{2}
\]

For the registered exact model,

\[
 \eta=\phi-1,\quad W=\tfrac13+\eta^2+\eta^3/3,
 \quad U=\tfrac12W'^2-\tfrac23W^2,
 \quad\sigma=2W+\delta(1+c+c\eta),
 \quad c=\frac2{1.0357712571566784}-\frac43,
\]

where the decimal is rational. The algebra here uses an arbitrary smooth
potential, not the September binary64 value of \(c\). The frozen theorem
and the prospective bridge supply, at every fixed \(0<\delta\le0.002\),

\[
 p<0\ (0<y\le T),\quad
 B:=p'/p+\sigma''/2>J:=p_\beta+\sigma''/2\ge w/2=1,
 \quad M_0>0,\quad\det D(R,S)=\mathcal H_bM_0>0.          \tag{3}
\]

Only the reflection-even interpretation is assumed until Section 7. A
doubled background by itself must not silently be treated as a restriction
on the perturbations; if the model is a genuine two-sided interface without
a reflection quotient, its odd sector also has to be included.

## 2. Unreduced covariant constraints and the generic operator

Write a scalar perturbation with \(\Box_\gamma Y=\mu^2Y\) as

\[
 \delta g_{yy}=2NY,\quad \delta g_{y\mu}=\rho^2 B_gD_\mu Y,
 \quad\delta g_{\mu\nu}=2\rho^2(\psi\gamma_{\mu\nu}Y
                         +E D_\mu D_\nu Y),\quad
 \delta\phi=\chi Y.
\]

The subscript on \(B_g\) distinguishes the shift perturbation from (3).
Under the convention \(\delta g=\mathcal L_\xi g\),
\(\xi^y=TY,\ \xi^\mu=L D^\mu Y\), the shifts are

\[
 (N,B_g,\psi,E,\chi)\mapsto
 (N+T',B_g+L'+\mathcal HT,\psi+aT,E+L,\chi+pT).
\]

For independent scalar harmonic tensors define
\(S_g=\rho^2(B_g-E')\),
\(\widetilde N=N-S_g'\),
\(\Psi=\psi-aS_g\), \(X=\chi-pS_g\). These equal the corresponding
fields in longitudinal gauge \(B_g=E=0\). Define
\(d=\Psi'-a\widetilde N\) and \(\Lambda=\mu^2+4\).

The radial ADM Hamiltonian and momentum equations follow from

\[
 G_{nn}=\tfrac12(K^2-K^\mu{}_{\nu}K^\nu{}_{\mu}-R^{(4)}),
 \quad G_{n\mu}=D_\nu K^\nu{}_{\mu}-D_\mu K.
\]

In longitudinal gauge \(\delta K^\mu{}_{\nu}=dY\delta^\mu{}_{\nu}\),
and \(\delta R^{(4)}=-6\mathcal H(\Box_\gamma+4)\Psi Y\).
Consequently the unreduced first-order constraints are

\[
 (-3d-pX)D_\mu Y=0,                                    \tag{4}
\]
\[
 12ad+3\Lambda\mathcal H\Psi
          -p(X'-p\widetilde N)+U'X=0.                    \tag{5}
\]

The trace-free tangential Einstein equation is

\[
 -(\widetilde N+2\Psi)
 \left(D_\mu D_\nu Y-\tfrac14\gamma_{\mu\nu}\Box_\gamma Y\right)=0.
                                                        \tag{6}
\]

For \(D Y\ne0\) and nonzero trace-free Hessian, (4) and (6) give
\(d=-pX/3\), \(\widetilde N=-2\Psi\). Substituting into (5) gives

\[
 \Psi'=-2a\Psi-pX/3,\qquad
 pX'=p'X+(3\Lambda\mathcal H-2p^2)\Psi.                  \tag{7}
\]

This is exactly the archived system, with \(Z=\Psi/p\). The division by
\(p\) is legitimate on the open interval for (3), but not at the cone or
at the undeformed constant-scalar limit.

For an explicit check that (7) is not just a necessary subset, the new
script independently forms five-dimensional Christoffel variations,
linearized Ricci, and the linearized scalar wave operator for
\(Y=e^{rt}\) on the flat de Sitter chart. The trace-free coefficient is
\(-(N+2\psi)\), the momentum coefficient is
\(3aN-3\psi'-pX\), and the Hamiltonian is (5). With
\(r^2+3r+\mu^2=0\), every nonzero Einstein component and the scalar
equation vanish after (7). The covariant ADM argument justifies extending
these coefficient equations beyond that concrete harmonic.

## 3. Full displaced-shell junction before harmonic division

Let the shell be \(y=T+\zeta Y\). Its gauge-invariant displacement for a
generic scalar harmonic is \(\xi=\zeta+S_g(T)\). The directly computed
mixed extrinsic curvature in longitudinal gauge is

\[
 \delta K^\mu{}_{\nu}=
  (\Psi'-a\widetilde N+a'\xi)Y\delta^\mu{}_{\nu}
           -\mathcal H\xi D^\mu D_\nu Y.                \tag{8}
\]

The induced scalar is \(X+p\xi\). Applying (2) and the two background
junctions, the complete linear metric and scalar residuals are

\[
 \boxed{\mathcal J^\mu{}_{\nu}=
  (\Psi'-a\widetilde N+pX/3)Y\delta^\mu{}_{\nu}
       -\mathcal H\xi(D^\mu D_\nu Y+\delta^\mu{}_{\nu}Y),}
                                                        \tag{9}
\]
\[
 \boxed{\mathcal J_\phi=X'-p\widetilde N+\sigma''X/2+pB\xi.}
                                                        \tag{10}
\]

Equation (9) is a tensor identity; it has not divided by a trace-free
harmonic. In a general scalar gauge the same statement follows by inserting
\(S_g\) in the fields and displacement. The direct normal calculation uses
\(n_y=1+NY\), \(n_\mu=-\xi D_\mu Y\), and evaluation on the displaced
surface. Its time and spatial components are independently checked in the
script. Keeping the background-shell evaluation shift is essential.

For generic harmonics, the trace-free part of (9) enforces \(\xi=0\).
The remaining metric residual is precisely (4). Equation (10) reduces to

\[
 M=X'+\frac{\sigma''}2X+2p\Psi
  =BX+\frac{3\Lambda\mathcal H}{p}\Psi=0.                \tag{11}
\]

In particular, dropping \(2p\Psi\) is not an approximation to the same
physical boundary condition. The script detects that error, a reversed
bending sign, and a reversed gravity term.

## 4. Which harmonics are actually exceptional

Define \(T_{\mu\nu}[Y]=D_\mu D_\nu Y-
\gamma_{\mu\nu}\Box Y/4\). Commuting covariant derivatives on unit
de Sitter gives

\[
 D^\mu T_{\mu\nu}[Y]=\tfrac34(\mu^2+4)D_\nu Y.         \tag{12}
\]

Thus the degeneracy is classified locally as follows.

* If \(DY=0\), the nonzero harmonic is constant and \(\mu^2=0\).
* If \(T[Y]=0\) but \(DY\ne0\), then \(\mu^2=-4\) and
  \(D_\mu D_\nu Y=-\gamma_{\mu\nu}Y\). This is the
  translation-type subspace, locally the five embedding-coordinate
  harmonics on \(dS_4\).
* All other harmonics satisfy the generic equations, including nonconstant
  solutions at \(\mu^2=0\) and nontranslation solutions at \(\mu^2=-4\).

For the concrete homogeneous realization, the exceptional exponents are
\(r=0,1\), while \(r=-3,-4\) have the same two respective mass values
and nonzero trace-free Hessian. It would be wrong to delete either whole
mass value from the generic spectral problem.

### 4.1 Constant harmonics

Here the momentum equation (4) vanishes identically, and the scalar
longitudinal construction based on independent gradient/Hessian tensors is
not available. Use radial Gaussian normal gauge instead. Regularity fixes
\(\rho(0)=0,\rho'(0)=1,p(0)=0\). The sole physical regular bulk parameter
is \(u=\phi(0)\), together with the endpoint displacement \(\delta T\);
the radial-origin shift is gauge, not a further background parameter.

The frozen regular-cone family therefore spans the admissible constant
linear deformation by \(\delta u\) and \(\delta T\). For
\(R=a-\sigma/6\), \(S=p+\sigma'/2\), the full static junction linearization
is

\[
 \begin{pmatrix}\delta R\\\delta S\end{pmatrix}
 =\begin{pmatrix}j&-\mathcal H_b\\ hJ&pB\end{pmatrix}
   \begin{pmatrix}\delta u\\\delta T\end{pmatrix},
 \quad\det=\mathcal H_bM_0>0.                            \tag{13}
\]

Thus neither a physical constant scalar zero mode nor a static relative
shell displacement survives. This conclusion uses the regular-cone domain
and fixed model parameters \((\delta,c,U,\sigma)\): varying a Lagrangian
parameter is not a perturbative mode. The nonzero determinant does not by
itself address the nonconstant \(\mu^2=0\) sector; that sector belongs to
(7), (11) and the separate nodeless-reference argument.

### 4.2 Translation-type harmonics

For \(D D Y=-\gamma Y\), the nominal \(E\) field only enters through
\(\bar\psi=\psi-E\). Choose Gaussian normal gauge \(N=B_g=0\). Do not
use (6) to impose \(N=-2\psi\): its harmonic tensor is zero.

The unreduced momentum and Hamiltonian constraints become

\[
 \bar\psi'=-p\chi/3,\qquad p\chi'=p'\chi.              \tag{14}
\]

Since \(p\ne0\), their complete local solution on \((0,T]\) is

\[
 \chi=pe_0,\qquad
 \bar\psi=(a+F)e_0-L_0,\qquad F'=\mathcal H.             \tag{15}
\]

It is exactly the residual Gaussian-normal diffeomorphism with
\(T=e_0\), \(L=L_0-e_0F\). Its two integration constants span the
first-order system because the field-value determinant is \(-p\ne0\).
These are solutions of all the linear equations by diffeomorphism
invariance. This proves exhaustion, not merely that one gauge family happens
to be a solution.

In the usual smooth regular-cone interpretation, \(F=-1/y+O(y)\) up to
its arbitrary additive constant. The leading singularities in \(a+F\)
cancel. The smooth translation gauge vector extends as an ordinary
translation of the regular apex; the other constant is fixed by the
smoothness/gauge convention. A singular or nonvanishing boundary symmetry
cannot simply be declared gauge in a different global state space. Here the
quotient is the usual one by smooth diffeomorphisms of the bulk and shell,
with no physical apex fixed independently of the fields.

The Hessian factor in (9) is identically zero. Substituting (15) into (10)
gives

\[
 \mathcal J_\phi=pB(\zeta+e_0).                          \tag{16}
\]

The relative displacement \(\Delta=\zeta+e_0\) is invariant under that
radial translation. Equation (3) forces \(\Delta=0\). The remaining
combined bulk-plus-shell perturbation is pure gauge; there is no physical
relative translation mode in this even sector. The result is conditional
on the stated smooth-diffeomorphism interpretation, not on a claim that
every formal coordinate transformation is an allowed gauge operation.

## 5. Sturm–Liouville conversion and the physical kinetic coefficient

Set

\[
 P=\rho^2\Psi,\quad p_{SL}=\frac1{\rho^2p^2},
 \quad q_{SL}=\frac2{3\rho^2},\quad
 w_{SL}=\frac1{\rho^4p^2},\quad
 \beta=\rho_b^4p_b^2B_b.
\]

Equations (7), (11) give

\[
 P'=-\rho^2pX/3,\quad
 -(p_{SL}P')'+q_{SL}P=\Lambda w_{SL}P,
 \quad p_{SL}P'(T)=\frac\Lambda\beta P(T).                \tag{17}
\]

For \(B>0\), the abstract norm and form are

\[
 \mathcal N_{SL}=\int w_{SL}|P|^2dy+|P(T)|^2/\beta,
 \qquad \mathcal E=\int(p_{SL}|P'|^2+q_{SL}|P|^2)dy.
\]

The equality \(\mathcal E=\Lambda\mathcal N_{SL}\) uses the boundary
condition and vanishing cone Green term. Positive \(\mathcal N_{SL}\)
proves a mathematical self-adjoint/energy statement for this realization;
it does not specify the normalization or even the sign of a physical action
that happens to have the same equation of motion.

To check that missing issue locally, include GHY in (1), use radial ADM
with spacelike normal, zero shift, lapse \(n\), and induced metric \(h\):

\[
 S_{\rm bulk+GHY}=\int\sqrt{-h}\,n\left[
 \frac{R[h]}2+\frac{K^2-K_{\mu\nu}K^{\mu\nu}}2
 -\frac{(\phi')^2}{2n^2}-\frac12h^{\mu\nu}\partial_\mu\phi\partial_\nu\phi-U
 \right].                                               \tag{18}
\]

Put \(h_{\mu\nu}=\rho^2e^{2\psi}\gamma_{\mu\nu}\),
\(n=e^N\), and \(\phi=\phi_0+\chi\). These nonlinear completions
only fix notation for the quadratic expansion. The intrinsic curvature is

\[
 R[h]=\rho^{-2}e^{-2\psi}
        [12-6\Box_\gamma\psi-6(D\psi)^2].
\]

After tangential integration by parts, the terms quadratic in perturbations
with tangential derivatives, for one copy, are exactly

\[
 \mathcal L^{(2)}_{\rm tangential}/\sqrt{-\gamma}
 =\rho^2[3D N\cdot D\psi+3(D\psi)^2-\tfrac12(D\chi)^2].
                                                        \tag{19}
\]

The extrinsic-curvature and normal-scalar terms have no tangential
derivatives in this zero-shift gauge. On the generic solution gauge slice
\(N=-2\Psi,\psi=\Psi,\chi=X\), (19) becomes

\[
 -\frac{\rho^2}{2}[(DX)^2+6(D\Psi)^2].                  \tag{20}
\]

For separated profiles \(X(y)Q(x),\Psi(y)Q(x)\), the doubled coefficient
of \(-\tfrac12\int\sqrt{-\gamma}(DQ)^2\) is the first quantity in (A).
For the even generic mode the physical shell is straight in longitudinal
gauge by (9), so \(-\sqrt{-h}\sigma\) adds no tangential derivatives.
There is no induced shell Einstein term in this model. Fixed-radial shell
GHY is already included in (18). Radial surfaces and constant-time surfaces
are orthogonal in this zero-shift gauge, so no varying corner angle adds
a kinetic term in this calculation.

Substituting \(P'= -\rho^2pX/3\) gives the exact density identity

\[
 2\rho^2(|X|^2+6|\Psi|^2)
     =18(p_{SL}|P'|^2+q_{SL}|P|^2),                      \tag{21}
\]

and hence (A) at eigenroots. The factor \(\Lambda\) is essential. In
particular, a blanket identification of \(\mathcal N_{SL}\) with the
physical kinetic norm would be false. Equation (20) supplies the ordinary
healthy Lorentzian kinetic sign for finite-energy, nonzero generic even
profiles on this gauge slice. It is a pullback of the action's quadratic
kinetic bilinear to constrained solutions; it is not an unconstrained
single-master-field action valid off shell for arbitrary \(P(x,y)\).

The remaining limitations of this kinetic calculation are precise:

* A global proof that every admitted physical state has such a regular
  longitudinal representative, and that the physical covariant symplectic
  flux has the prescribed zero-cone boundary limit, is separate work. For
  the prospective nonpositive modes, the known cone powers make the
  displayed radial integrals and Green terms harmless; this does not by
  itself prove completeness of the full physical state space.
* Exceptional harmonic sectors must use Section 4, where the constraint
  \(N=-2\Psi\) cannot be imposed by the generic argument.
* A reflection-odd mode has a moving shell in longitudinal gauge. Formula
  (20) gives its bulk contribution only; it cannot be promoted to a full
  physical norm without the moving-shell/embedding/GHY contribution.

Thus the audit improves the local kinetic understanding but should not be
advertised as an unrestricted no-ghost theorem.

## 6. What the available signs actually discharge

The nonconstant branch supplies all algebraic nondegeneracy hypotheses
needed for Sections 2–4: \(p\ne0\), \(B>0\), and the nonzero static
determinant. Together with the independently constructed nodeless reference
\(P_0\), the generic spectral problem admits the positive Picone identity
recorded in the prospective bridge. Subject to the separately audited radial
form domain, this excludes generic normalizable \(\mu^2\le0\) modes.

The signs and determinant also discharge the two exceptional *even*
subspaces under regular-cone and smooth-diffeomorphism assumptions:
constant deformations by (13), and translation-type relative displacements
by (16). Neither statement requires a floating-point scan, a lower bound
on the tiny central \(|p|\), or replacing the exact-decimal model with
the September numerical coupling.

This is a scalar separated-mode conclusion. It does not establish a
uniform-in-\(\delta\) physical mass gap, vector-sector stability, nonlinear
attraction, vacuum-decay suppression, or a general well-posed semigroup
estimate. The divided master variable is singular at \(\delta=0\), which
remains outside the argument.

## 7. Two-sided interface caveat and independent odd-sector cross-check

The other independent audit in this work session treats the case in which
the two interiors may fluctuate independently. Its generic odd continuity
conditions are
\(X_o+p\zeta_o=0\), \(\Psi_o+a\zeta_o=0\), while the summed
trace-free junction gives \(\zeta_++\zeta_-=0\). They imply

\[
 X_o=(p/a)\Psi_o,\qquad
 p_{SL}P_o'(T)=-\frac{1}{3a_b\rho_b^2}P_o(T).            \tag{22}
\]

The sign in (22) follows directly from (17) and was independently checked
here. It agrees with the one-side outward-normal residual (9). It is a
different boundary realization, with no positive endpoint mass of the even
problem. The independent odd-sector audit, rather than this note, is
responsible for completeness of those matching conditions and their
exceptional subspaces.

For the odd kinetic norm, write \(c_o=1/(3a_b\rho_b^2)\). The bulk
identity alone gives \(K_{\rm bulk}=18\mathcal E_o
=18\Lambda\int w_{SL}|P_o|^2-18c_o|P_o(T)|^2\).
The shell volume bending term by itself would contribute
\(\sigma\rho_b^2\zeta_o^2=18c_oP_o(T)^2\), using \(\sigma=6a\),
which suggests a cancellation. This is **not** a derived full odd kinetic
norm: moving integration boundaries and the displaced GHY/corner terms
must first be included consistently. A sign inferred from the partial
cancellation would overstate the present work.

## 8. Reproduction and audit boundaries

Run from this directory:

```sh
PYTHONDONTWRITEBYTECODE=1 /workspace/hdblast-research-work/environment/bin/python check_operator_20261013.py
```

The run used SymPy 1.14.0 and passed 32/32 exact checks. The checks include
three deliberate wrong-formula controls. The script has no imports from
the archived derivation. It records SHA-256 hashes of the historical sources
and its own source in `OPERATOR_CHECK_RESULTS_20261013.json`; the complete
console output is `OPERATOR_CHECK_LOG_20261013.txt`.

The count is an inventory of exact local algebra checks, not an independent
proof of every analytic assumption in this note. The theorem about the
regular-cone family and its Jacobi signs is imported from the frozen proof;
the radial form-domain theorem is being audited separately. No archived
numerical spectrum was rerun, and no new interval spectral certificate is
claimed.
