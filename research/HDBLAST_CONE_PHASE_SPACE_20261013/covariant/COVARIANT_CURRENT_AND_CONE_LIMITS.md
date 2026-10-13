# Einstein–scalar symplectic current, shell balance, and two cone limits

13 October 2026. Internal continuation of the frozen regular-shell spectral
checkpoint. This note changes no repository or frozen source. The results
below concern a stated generic scalar solution class. They do **not** prove
that this class exhausts the global physical phase space.

## 1. What is new and what remains conditional

The Einstein–scalar Lee–Wald current can be pulled back explicitly to the
longitudinal scalar solution slice. For a common separated spectral value
`Lambda`, its radial component is **9 Lambda times the reduced Green current
per interior**, hence `18 Lambda` for a doubled parity perturbation. This
is a current identity derived from the action, not an identification of a
radial Hilbert norm with a physical norm by definition.

For different spectral values, or for packets, the common-Lambda formula
must be replaced by (13) below. The resulting current has an exact shell
balance: no shell current for the straight even shell, and the positive
moving-shell current (17) for the odd shell.

The same current has two inequivalent endpoint limits. Smooth perturbations
have zero flux through a shrinking radial surface over **bounded de Sitter
time**. A fixed nonzero segment of the characteristic cone instead requires
`y -> 0`, `t -> +infinity`, `v=y exp(t)` fixed. For actual generic solutions
with the uniform regular characteristic expansion in Section 7, the limiting
flux is the ordinary scalar null symplectic flux, which is not identically
zero as a bilinear on the allowed leading data. No global solution with those
data and the distant shell conditions is constructed here. This establishes
the current and the conditional obstruction to an unrestricted zero-flux
claim; it is not a global existence or completeness theorem.

## 2. Action, orientation, and potential ambiguity

Use signature `(-++++)`, `8 pi G_5=1`, unit `dS_4`, and

\[
 S=\sum_{i=\pm}\left\{\int_{M_i}\epsilon_g
 [R/2-(\nabla\phi)^2/2-U(\phi)]
 +\int_\Sigma\epsilon_h K_i\right\}
 -\int_\Sigma\epsilon_h\sigma(\phi).                 \tag{1}
\]

Each shell normal points toward increasing interior radius. A displayed
`j^y` is a coordinate vector **density**, oriented toward increasing `y`.
The outward flux through the lower boundary `y=epsilon` of an interior
`epsilon<y<T` is therefore `-j^y`. Every flux sign below uses this convention.
Real bilinears are written; complex solutions can be paired by conjugating
one argument after complexification.

For `h_AB=delta g_AB`, `h=g^{AB}h_AB`, `chi=delta phi`, the exact bulk
presymplectic potential is

\[
 \theta^A(g,\phi;h,\chi)=\sqrt{-g}\left[
 \frac12(\nabla_B h^{AB}-\nabla^A h)-\nabla^A\phi\,\chi\right],
 \qquad
 \omega^A(1,2)=\delta_1\theta^A(2)-\delta_2\theta^A(1).
                                                               \tag{2}
\]

The variations commute. In particular the exact scalar part is

\[
 \omega^A_\phi=\sqrt{-g}\big[
 \chi_1\nabla^A\chi_2-\chi_2\nabla^A\chi_1
 +(h_1^{AB}-h_1g^{AB}/2)\partial_B\phi\,\chi_2
 -(h_2^{AB}-h_2g^{AB}/2)\partial_B\phi\,\chi_1\big].       \tag{3}
\]

The gravitational part is obtained by antisymmetrizing the variation of
the first term in (2), including the volume density and raised indices.
Formula (2), rather than a convention-dependent six-index tensor, fixes it
completely. On linearized solutions, `partial_A omega^A=0`.

Adding a boundary Lagrangian changes `theta` by its field variation and
does not change the bulk `omega`. Changing `theta` by a spacetime exterior
derivative changes `omega` by a corner term. Thus charges on finite cuts
require a specified potential/corner convention. Here the bulk convention
is exactly (2); shell results are obtained from (1), including its two GHY
terms and moving embeddings. No counterterm is added at the cone.

## 3. GHY cancellation and a direct pullback of the covariant current

Write `n=exp(N)`, `h_mn=rho^2 exp(2 psi) gamma_mn`, with zero radial shift,
and `phi=phi_0+chi`. This nonlinear family is used only to evaluate the
second variation at the background. Its radial curvature is
`K=4(a+psi')/n`. Direct substitution in (2) gives

\[
 \theta^y=\sqrt{-h}\left[K(\delta N-\delta\psi)
                  -4\delta\psi'/n-\phi'\delta\phi/n\right].
\]

Consequently, throughout this conformal zero-shift family,

\[
 \theta^y+\delta(\sqrt{-h}K)
       =3\sqrt{-h}K\,\delta\psi
        -\sqrt{-h}\phi'\delta\phi/n
       =\pi^{mn}\delta h_{mn}+\pi_\phi\delta\phi,       \tag{4}
\]
\[
 \pi^{mn}=\frac{\sqrt{-h}}2(Kh^{mn}-K^{mn}),
 \qquad \pi_\phi=-\sqrt{-h}\phi'/n.
\]

This proves that the radial canonical and covariant currents agree
**pointwise on this slice**: the field-space curl of the GHY variation in
(4) vanishes. A missing corner improvement is not being assumed away in
the radial formula. For unrestricted variations with nonzero shift, the
usual tangential corner term must be retained.

At the background, antisymmetrizing (4) yields, per copy,

\[
 \frac{\omega^y}{\sqrt{-\gamma}\rho^4}
 =12[(\psi_1'-aN_1)\psi_2-(\psi_2'-aN_2)\psi_1]
  +\chi_1\chi_2'-\chi_2\chi_1'
  +p(N_1\chi_2-N_2\chi_1)
  -4p(\psi_1\chi_2-\psi_2\chi_1).                     \tag{5}
\]

The tangential components follow directly from (2) as well. Before
linearization, the gravitational potential in this family is

\[
 \theta^\mu_g=n\sqrt{-h}\big[
 -D_h^\mu\delta N-3D_h^\mu\delta\psi
 +(\delta\psi-\delta N)D_h^\mu N\big].                \tag{6}
\]

In particular (5) and (6) do not rely on guessing an action for a divided
master variable.

## 4. Generic constraints and the exact scalar current

For the generic scalar harmonic subspace, use the longitudinal
representative with `N=-2 Psi`, `psi=Psi`, `chi=X`. The full bulk equations
give

\[
 \Psi'=-2a\Psi-pX/3,\qquad
 pX'=p'X+3\rho^{-2}(\Box_\gamma+4)\Psi-2p^2\Psi.       \tag{7}
\]

No inverse of `Box+4` is used here. The exceptional constant and
Hessian-degenerate subspaces do not satisfy these generic implications
by the same argument and must be treated separately.

Substitution of (7) into (5), and direct linearization of (6), give the
following **Lee–Wald current**, per interior:

\[
 \boxed{\omega^y=\sqrt{-\gamma}\rho^4
 [X_1X_2'-X_2X_1'-2p(\Psi_1X_2-\Psi_2X_1)]},          \tag{8}
\]
\[
 \boxed{\omega^\mu=\sqrt{-\gamma}\rho^2
 [X_1D^\mu X_2-X_2D^\mu X_1
   +6(\Psi_1D^\mu\Psi_2-\Psi_2D^\mu\Psi_1)]}.          \tag{9}
\]

The cancellation producing the `-2p` term in (8) uses both gravity and
scalar variations; the scalar current alone has a different coefficient.
Pullback of the quadratic action gives exactly (9), explaining the earlier
local kinetic identity while now fixing the full spacetime current.

Set

\[
 P=\rho^2\Psi,\quad
 X=-3P'/(\rho^2p),\quad
 \alpha=(\rho^2p^2)^{-1},\quad q=2/(3\rho^2),\quad
 w=(\rho^4p^2)^{-1},\quad \mathcal D=\Box_\gamma+4.
                                                               \tag{10}
\]

For separated solutions `P_i(y)Y_i(x)`, `Box Y_i=(Lambda_i-4)Y_i`,

\[
 \omega^y=9\sqrt{-\gamma}\alpha
 [\Lambda_1P_1P_2'-\Lambda_2P_1'P_2]Y_1Y_2,           \tag{11}
\]
\[
 \omega^\mu=9\sqrt{-\gamma}
 [\alpha P_1'P_2'+qP_1P_2]
 (Y_1D^\mu Y_2-Y_2D^\mu Y_1).                          \tag{12}
\]

For general superpositions the radial formula is

\[
 \boxed{\omega^y=9\sqrt{-\gamma}\alpha
 [(\mathcal D P_1)P_2'-(\mathcal D P_2)P_1']}.          \tag{13}
\]

Thus `omega^y=9 Lambda sqrt(-gamma) alpha W_y(P_1,P_2)` only when both
radial factors have the same `Lambda`. All expressions double for two
pure-even or two pure-odd perturbations with the same one-side amplitude.
The sum of cross even/odd bulk bilinears cancels.

An independent conservation check uses only
`(alpha P_i')'=(q-Lambda_i w)P_i`:

\[
 \partial_y\{\alpha[\Lambda_1P_1P_2'
                   -\Lambda_2P_1'P_2]\}
 =(\Lambda_1-\Lambda_2)
       (\alpha P_1'P_2'+qP_1P_2).
                                                               \tag{14}
\]

This cancels the divergence of (12). Applying the equal-Lambda expression
to arbitrary packets would fail this check.

## 5. The two shell contributions and their conservation law

On a fixed shell the GHY-completed first variation of each bulk copy is
`pi_i^{mn} delta h_mn + pi_phi,i delta phi`, modulo the ordinary corner term.
The shell potential has no intrinsic derivatives of `h_mn` or `phi`.
Adding `-delta(sqrt(-h) sigma)` gives the junctions

\[
 \sum_i(K^{(i)}_{mn}-K_i h_{mn})+\sigma h_{mn}=0,
 \qquad \sum_i p_i+\sigma'=0.                         \tag{15}
\]

The common induced fields, rather than independent fictitious shell
fields, are varied here. When the embedding moves, the graph GHY terms,
the shifted bulk domain and the Einstein tangential divergence are needed.
The frozen odd-completion calculation supplies their full quadratic
pullback; polarizing it is sufficient to fix the current on the current
generic solution class.

For clarity, the complete accounting behind that pullback is as follows.
For one graph `y=T+zeta(x)`, graph GHY minus slice GHY together with the
Einstein tangential divergence contributes
`+3 a_b rho_b^2 (D zeta)^2`; the shifted radial integration domain
contributes `+3 rho_b^2 D zeta . D psi_b`. Both terms occur for each
of the two interiors. The single common shell contributes
`-sigma_b rho_b^2 (D zeta)^2/2`. Odd continuity is
`psi_b=-a_b zeta`, `X_b=-p_b zeta`, and `sigma_b=6a_b`.
Thus the sum is `-3 a_b rho_b^2 (D zeta)^2`, or (16) after
`zeta=-P_b/(a_b rho_b^2)`. Dropping the moving domain would reverse
this conclusion. These terms also display why merely taking the shell
volume current without the gravitational edge contributions is not a
derivation, although their final values coincide after continuity.

For the even sector the generic shell is straight and
`alpha P_i'(T)=Lambda_i P_i(T)/beta`. Equation (11) then gives
`omega^y_+(T)+omega^y_-(T)=0` even for unequal `Lambda_i`.
There is no additional tangential kinetic shell current on this slice.
The endpoint mass in the *reduced* Hilbert norm is not an extra physical
shell kinetic field in this gauge.

For the odd sector, set `c_o=1/(3 a_b rho_b^2)>0` and use
`alpha P_i'(T)=-c_o P_i(T)`. The complete moving-shell action is

\[
 S^{(2)}_{\Sigma,kin}
   =-9c_o\int\sqrt{-\gamma}(D P_b)^2,                 \tag{16}
\]

where `P_b(x)` is the plus-side trace. It includes both gravitational
copies and one common shell, not two shell tensions. Its current is

\[
 \omega^\mu_\Sigma=18c_o\sqrt{-\gamma}
 [P_{1b}D^\mu P_{2b}-P_{2b}D^\mu P_{1b}],
 \qquad
 \partial_\mu\omega^\mu_\Sigma
       =\omega^y_+(T)+\omega^y_-(T).                   \tag{17}
\]

Consequently
`partial_mu[integral_0^T(omega_+^mu+omega_-^mu)dy+omega_Sigma^mu]`
equals the lower-end `+y` current. This verifies conservation with the
sign of the completed odd kinetic action. At a common spectral root,
the integrated coefficients are exactly the known
`K_e=18 Lambda N_e` and `K_o=18 Lambda N_o`.

The current (17) is the fully completed shell current on this constrained
generic scalar class. It is not advertised as a formula for arbitrary
unreduced shell embeddings or for independent edge degrees of freedom.

## 6. Finite-time shrinking surfaces: a bounded zero-flux statement

Assume the regular cone expansions
`rho=y(1+O(y^2))`, `p=kappa y(1+O(y^2))`, `kappa!=0`, with differentiated
remainders. Introduce global de Sitter time
`gamma=-dt^2+cosh^2(t)dOmega_3^2`.

If the spacetime metric, scalar, perturbations and first derivatives extend
smoothly in regular ambient coordinates near the vertex, the geometric
current in (2) is bounded there. On `y=epsilon` and any fixed compact
`t` interval the proper surface volume is `O(epsilon^4)`. Its integrated
flux therefore tends to zero. The statement permits smooth angular data
on the complete `S^3`, or compact chart support, and uses bounded time
intervals. It is a geometric statement independent of the divided variable.

One can also state a sufficient radial-operator regularity criterion.
Let `x=int_y^T dt/rho(t)` and `F=P/(rho^(3/2)|p|)`. Then

\[
 \alpha W_y(P_1,P_2)=-W_x(F_1,F_2).                    \tag{18}
\]

For functions in the Liouville operator domain `H^2` near infinity,
both `F` and `F_x` tend to zero. This proves the common-spectral-value
Green limit, not a characteristic assertion. For the packet expression
(13), put `G=(D P)/(rho^(3/2)|p|)` and
`s_L=rho(3a/2+p'/p)`. Then

\[
 \frac{\omega^y}{9\sqrt{-\gamma}}
 =-G_1F_{2x}+G_2F_{1x}
     +s_L(G_1F_2-G_2F_1).                              \tag{19}
\]

For example, `F_i,G_i` in `H^2` with uniform tail bounds over a compact
time/angular slab are sufficient for zero integrated radial flux there.
Taking `P` and `D P` in the radial operator domain is a convenient strong
condition. The mere form domain `H^1` does not provide a pointwise trace
of `F_x`; no such trace is asserted for every finite-form-energy packet.
Weak flux/conservation statements for completed energy spaces require
their own trace or approximation argument.

## 7. A different, characteristic limit and its actual physical flux

Near the future cone use

\[
 u=ye^{-t},\quad v=ye^t,\quad y^2=uv,
 \qquad ds^2_{flat}=du\,dv+\frac{(u+v)^2}{4}d\Omega_3^2.
                                                               \tag{20}
\]

The future cone is `u=0`, `v>0`. A compact nonzero cone segment
`v in [v_1,v_2]`, `v_1>0`, is approached with `t=log(v/y)`, not bounded
`t`. The regular background has a smooth continuation in these coordinates:
`rho^2/y^2=1+O(uv)` and `phi=phi_c+kappa uv/2+O((uv)^2)`.

Here is a conditional theorem directly about the physical current.
Suppose two **actual generic linearized solutions** have, uniformly with
the derivatives used below on this compact cone segment,

\[
 P_i=y^4 A_i(v,\Omega)+O(u^3),\quad
 \Psi_i=uv A_i(v,\Omega)+O(u^2),\quad
 X_i=X_{0i}(v,\Omega)+O(u).                             \tag{21}
\]

The remainders are understood in smooth `(u,v,Omega)` coordinates, and
are at least `C^1` for current limits and `C^2` where the equations are
invoked. Since `v` lies in a compact interval away from zero, replacing
the first remainder by `O(y^6)` is equivalent at the level of orders.
The radial momentum constraint fixes

\[
 X_{0i}=-\frac3\kappa(4A_i+v\partial_v A_i).           \tag{22}
\]

This ansatz is compatible with smooth metric components at the
characteristic surface: to leading order the longitudinal perturbation is

\[
 \delta g_i=2\Psi_i g-6\Psi_i dy^2
 =2uv A_i g-\frac32 A_i(v\,du+u\,dv)^2+O(u^2)
                                                               \tag{23}
\]

in the corresponding smooth tensor expansion. The tensor remainder in
(23) is schematic: substituting an `O(u^2)` scalar remainder in `Psi`
can give an `O(u)` coefficient of `du^2`, still smooth. Smoothness of the
leading expression is what is asserted. No smoothness is inferred just
from a bound on the divided variable without its derivatives.

Equation (8), not an equal-Lambda Wronskian shortcut, gives per interior

\[
 \omega^y\longrightarrow
 \frac{v^4}{8}\sqrt{\gamma_{S^3}}
 (X_{01}\partial_vX_{02}-X_{02}\partial_vX_{01}),
 \quad dt=\frac{dv}{v}\quad(y\text{ fixed}).             \tag{24}
\]

Indeed `rho^4 cosh^3t=(v^3/8)y+O(y^3)` and
`X_i'|_t=(v/y)partial_vX_{0i}+O(y)`. The metric-scalar cross term in
(8) is `O(y^4)` after the measure is included and vanishes. For two
copies with the same parity amplitude the integrated `+y` flux is

\[
 \boxed{\ \lim_{y\to0}\int\omega^y_{doubled}\,dt\,d\Omega
 =\frac14\int_{v_1}^{v_2}\!dv\,v^3\int d\Omega\sqrt{\gamma_{S^3}}
 (X_{01}\partial_vX_{02}-X_{02}\partial_vX_{01})\ }.       \tag{25}
\]

This is also the actual cone pullback, not merely a suggestive radial
scaling. A vector density transforms with Jacobian
`|partial(y,t)/partial(u,v)|=1/(2y)`. Hence

\[
 \omega^u=\frac1{2y}\left[\frac{u}{y}\omega^y-u\omega^t\right]
 \longrightarrow\frac{v^3}{8}\sqrt{\gamma_{S^3}}
 (X_{01}\partial_vX_{02}-X_{02}\partial_vX_{01})          \tag{26}
\]

per copy, using (9). It agrees with the scalar Klein–Gordon current on
`ds^2=du dv+(v^2/4)dOmega_3^2`, since
`sqrt(-g) g^{uv}=v^3 sqrt(gamma_S3)/8` on the cone. The outward flux from
the exterior wedge at its lower cone boundary is the negative of (25).
This is an antisymmetric symplectic bilinear, not a positive energy flux.

The bilinear on the right is not identically zero on smooth leading data.
For example take `X_01=f(v)H(Omega)`, `X_02=v f(v)H(Omega)` with a nonzero
smooth compactly supported `f` in `(v_1,v_2)`. Their Wronskian is
`f(v)^2 H(Omega)^2`. The corresponding leading `A_i` solve the regular
first-order equation (22) on this segment. This is an explicit nonzero
**candidate-data** flux, not a claimed solution of the complete shell
initial-boundary problem. The next-order constraints, global existence,
shell matching and admissibility of such data still need proof.

## 8. Gauge transformations, exceptional modes, and cut charges

For a field-independent infinitesimal diffeomorphism `xi`, the bulk
on-shell identity in the convention (2) is

\[
 \boldsymbol\omega(\delta,\mathcal L_\xi)
 =d\boldsymbol k_\xi(\delta),\qquad
 \boldsymbol k_\xi=\delta\boldsymbol Q_\xi
                  -i_\xi\boldsymbol\theta(\delta),
 \qquad \boldsymbol Q_\xi=-\tfrac12\star d\xi^\flat.
                                                               \tag{27}
\]

Matter does not supply an extra derivative-of-`xi` Noether charge in this
minimally coupled scalar model. A field-dependent `xi` requires the usual
`-Q_{delta xi}` adjustment and cannot be covered silently by (27).

A smooth ambient `xi` and smooth perturbation with bounded first
derivatives have bounded `k_xi`; its integral on a shrinking fixed-time
`S^3` is `O(y^3)` and vanishes. This supplies a sufficient local zero-charge
condition at the vertex. It does not prove that the same diffeomorphism
has zero charge on a finite-radius null cut. At a shell the embedding
must be varied together with the bulk metric/scalar, and its contribution
to the total presymplectic form must be included. Diffeomorphisms supported
away from all physical boundaries are degeneracies by (27). Transformations
acting on specified characteristic boundary data or on a retained null
cut must be checked by their total charge, not declared gauge from their
local tensor form.

The frozen constant-sector determinant and the local translation-type
spanning calculation remain valid under their original hypotheses. The
translation-type harmonics correspond locally to smooth ambient apex
translations. Equation (27) validates their absence of a *shrinking vertex*
charge for smooth perturbations. It does not by itself fix a global
characteristic charge or prove that moving the vertex is a gauge operation
in a state space with externally fixed cone data. No exceptional charge
sector is constructed or excluded globally here.

## 9. Regularity, time support, and the next genuinely missing proof

Compactly supported tangential **test variations** justify the local
quadratic action integrations by parts, with vanishing corner variations
near auxiliary initial/final surfaces. Nonzero on-shell hyperbolic
solutions must not be assumed compactly supported in time; the flux laws
above instead apply on finite slabs with their Cauchy-surface terms kept.
On noncompact de Sitter charts, spatial decay or compact chart support
must likewise be stated before discarding lateral flux. Global `S^3`
slices avoid a spatial boundary in the formulas here.

The next missing result is now concrete: choose a physical initial-data
space and an extension or boundary rule at the null cone; prove existence,
regular gauge representatives and shell compatibility in that space;
control the full current (13), (17), and (25) under completion; and compute
the charges (27) of the residual transformations. A radial limit-point
classification addresses only the radial realization. It does not supply
any of these characteristic decisions automatically.

## 10. Exact checks and source boundaries

`check_covariant_current.py` independently checks the local field-space
curl, constrained current, unequal-spectral conservation, both shell
balances, Liouville Green conversion, density Jacobian, smooth leading
metric components and null-current normalization. Deliberate wrong-formula
controls reject an omitted gravitational scalar cross term, a common-Lambda
shortcut for unequal profiles, a missing odd shell current, and a missing
density Jacobian. These are algebra checks, not proofs of global existence,
functional-analytic traces, or physical completeness.

Imported frozen sources are `operator/PHYSICAL_OPERATOR_AUDIT_20261013.md`,
`operator/ODD_KINETIC_COMPLETION_20261013.md`,
`domain/SCALAR_SELF_ADJOINT_DOMAIN_AND_GAP.md`, and the manuscript
`JACOBI_RESPONSE_TO_SPECTRAL_POSITIVITY.md`. In particular the odd graph
action coefficient, registered background signs, and radial domain theorem
are imported with their stated scope. The new current calculation uses
arbitrary smooth `U` and `sigma` and applies at each fixed nonconstant
registered detuning; it does not take the singular `delta=0` reduction.
