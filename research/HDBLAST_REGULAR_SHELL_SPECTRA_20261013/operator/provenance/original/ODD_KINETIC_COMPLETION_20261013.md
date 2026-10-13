# Local completion of the reflection-odd kinetic action

13 October 2026 Pacific and UTC. This bounded follow-up resolves the *local
moving-shell action gap* left open in Section 7 of
`PHYSICAL_OPERATOR_AUDIT_20261013.md`. It does not replace the separate global
gauge, cone-flux, and physical symplectic-domain assumptions in that note.
The earlier report is preserved as the completed first audit stage.

The complete quadratic kinetic coefficient for a generic odd mode is

\[
 \boxed{K_o=18\left(\mathcal E[P_o]
                  +\frac{|P_o(T)|^2}{3a_b\rho_b^2}\right)
       =18\Lambda\int_0^T w_{SL}|P_o|^2dy.}                \tag{1}
\]

The second equality is at a spectral root of the odd Robin problem. This
calculation includes the shifted bulk integration domain, the graph GHY
term, the bulk Einstein tangential-divergence flux, and the shell volume.
Dropping the moving-domain contribution gives the wrong kinetic sign.

## 1. Perturbations and the action calculation being performed

Use the EH+GHY+shell action and outward-normal conventions of the first
audit. For the two interiors put

\[
 g_i=n_i^2dy_i^2+h^{(i)}_{\mu\nu}dx^\mu dx^\nu,
 \quad h^{(i)}_{\mu\nu}
       =\rho^2e^{2\epsilon\psi_i(y_i)Q(x)}\gamma_{\mu\nu},
 \quad n_i=e^{\epsilon N_i(y_i)Q(x)},
 \quad \phi_i=\phi_0+\epsilon X_i(y_i)Q(x),
\]
\[
 y_i=f_i(x)=T+\epsilon\zeta_i Q(x).
\]

These are convenient nonlinear completions of the linear perturbations;
only the quadratic kinetic coefficient is retained. A purely odd mode has
\((\psi_-,N_-,X_-,\zeta_-)=-(\psi_+,N_+,X_+,\zeta_+)\).
Use \((\psi,N,X,\zeta)\) for its plus-side amplitudes. The generic bulk
constraint is \(N=-2\psi\). Induced metric and scalar continuity give

\[
 \psi_b=-a_b\zeta,\qquad X_b=-p_b\zeta.                  \tag{2}
\]

The two induced metrics coincide through the quadratic terms of these
completions: the linear terms vanish by (2), and the quadratic terms are
even under simultaneous sign reversal. The same observation holds for the
induced scalar. Possible second-order corrections to this choice of
completion multiply the first variation of the full action, which vanishes
on the background solution with its junctions. Thus they do not alter this
quadratic on-background action.

We compute the coefficient of
\(-\tfrac12\int\sqrt{-\gamma}\,K_o(DQ)^2\), allowing smooth
compactly supported test functions \(Q\) when integrations by parts are
performed. This is a local kinetic-bilinear calculation; the profiles may
subsequently be chosen to satisfy the separated equations. The use of test
functions does not require them to obey their four-dimensional wave equation
off shell during the calculation.

## 2. The radial ADM decomposition on a moving domain

With a zero shift and spacelike radial normal, the exact gravitational
density decomposes as

\[
 \frac12\sqrt{-g}R
 =n\sqrt{-h}\left[\frac{R[h]}2
                  +\frac{K_y^2-K_{y\mu\nu}K_y^{\mu\nu}}2\right]
   -\partial_y(\sqrt{-h}K_y)
   -\partial_\mu V^\mu,
 \quad V^\mu=\sqrt{-h}\,h^{\mu\nu}\partial_\nu n.       \tag{3}
\]

On a fixed radial boundary, GHY cancels the radial divergence. On
\(y=f(x)\), the cancellation is incomplete, and the tangential divergence
also has a shell contribution. Indeed,

\[
 -\int d^4x\int_0^{f(x)}dy\,\partial_\mu V^\mu
       =\int d^4x\,V^\mu(f(x),x)\partial_\mu f,          \tag{4}
\]

modulo the outer tangential boundary, where the perturbation is taken to
vanish. Therefore the exact expression relevant to the shell is

\[
 S_{EH+GHY}=\int_{y<f(x)}\mathcal L_{ADM}
 +\int d^4x\left[\sqrt{-h_{graph}}K_{graph}
                 -\sqrt{-h}K_y+V^\mu\partial_\mu f\right]_{y=f(x)}.
                                                        \tag{5}
\]

Using (5) avoids replacing the graph GHY term with the coordinate-slice
one. Its signs can be checked already in the simple warped-product formula
\(R_5=R[h]-K_y^2-K_{y\mu\nu}K_y^{\mu\nu}
-2n^{-1}\partial_yK_y-2n^{-1}\Box_h n\).

## 3. Direct normal and graph GHY expansion

Write \(f_\mu=\partial_\mu f\) and
\(D_f=1+n^2h^{\mu\nu}f_\mu f_\nu\). The unit outward normal and graph
measure are

\[
 n_{graph}^y=\frac1{n\sqrt{D_f}},\qquad
 n_{graph}^\mu=-\frac{n h^{\mu\nu}f_\nu}{\sqrt{D_f}},
 \qquad \sqrt{-h_{graph}}=\sqrt{-h}\sqrt{D_f}.
\]

Computing \(K_{graph}=\nabla_A n_{graph}^A\) gives

\[
 \sqrt{-h_{graph}}K_{graph}
 =\frac{\sqrt{D_f}}n\left[
     \partial_y\left(\frac{\sqrt{-h}}{\sqrt{D_f}}\right)
    -\partial_\mu\left(
       \frac{n^2\sqrt{-h}h^{\mu\nu}f_\nu}{\sqrt{D_f}}\right)
                         \right]_{y=f(x)}.              \tag{6}
\]

To second order, subtracting \(\sqrt{-h}K_y\), adding (4), and then
integrating tangentially gives the simple one-side result

\[
 \left(S_{EH+GHY}-\int_{y<f}\mathcal L_{ADM}\right)^{(2)}_{kin}
       =\int\sqrt{-\gamma}\,3a_b\rho_b^2\zeta^2(DQ)^2.
                                                        \tag{7}
\]

For transparency, before adding (4), the integrated graph-minus-coordinate
GHY term has coefficient
\(\rho_b^2(3a_b\zeta^2-N_b\zeta)\).
The tangential Einstein divergence supplies
\(+\rho_b^2N_b\zeta\), cancelling the lapse cross term. Neither term
alone gives (7). The second-order evaluation shift at \(y=f(x)\) is also
included in (6); omitting it changes the coefficient of \(a_b\zeta^2\).

The new script expands (6) independently in a homogeneous chart with
\(\sqrt{-\gamma}=e^{3t}\), retaining \(Q(t)\), \(\dot Q\), and
\(\ddot Q\). It obtains (7) after time integration by parts, without
substituting the odd continuity condition until all geometric contributions
have been combined. The coefficient is independent of \(N_b\) and
\(\psi_b\), as (7) requires.

## 4. Moving bulk domain and shell volume

The radial ADM intrinsic-curvature term is

\[
 \frac12 n\sqrt{-h}R[h],\qquad
 R[h]=\rho^{-2}e^{-2\epsilon\psi Q}
 [12-6\epsilon\psi\Box_\gamma Q
                  -6\epsilon^2\psi^2(DQ)^2].
\]

Its first-order tangential-derivative term is
\(-3\epsilon\sqrt{-\gamma}\rho^2\psi\Box Q\). Taylor expanding the
upper integration limit in (5) therefore produces, for one copy,

\[
 \left(\int_{y<f}\mathcal L_{ADM}
             -\int_0^T\mathcal L_{ADM}\right)^{(2)}_{kin}
 =\int\sqrt{-\gamma}\,3\rho_b^2\zeta\psi_b(DQ)^2.       \tag{8}
\]

The second-order shift acting on the background density has no tangential
derivatives. Normal scalar derivatives, the bulk potential, and radial
extrinsic terms likewise add no further \((DQ)^2\) term to (8).

The induced metric on the graph has the additional term
\(\epsilon^2\zeta^2\partial_\mu Q\partial_\nu Q\). Consequently
the single shell action contributes

\[
 S_{shell,kin}^{(2)}
      =-\frac12\int\sqrt{-\gamma}\,
                     \sigma_b\rho_b^2\zeta^2(DQ)^2.     \tag{9}
\]

At this order the scalar dependence of \(\sigma\) changes only terms
without tangential derivatives: the graph-gradient factor is already
quadratic.

Adding two copies of (7) and (8), then one copy of (9), gives

\[
 S_{boundary,kin}^{(2)}
 =\int\sqrt{-\gamma}\,\rho_b^2
       [(6a_b-\sigma_b/2)\zeta^2+6\zeta\psi_b](DQ)^2.
\]

Finally use \(\sigma_b=6a_b\) and (2):

\[
 S_{boundary,kin}^{(2)}
       =-3\int\sqrt{-\gamma}\,a_b\rho_b^2\zeta^2(DQ)^2,
 \qquad K_{boundary}=6a_b\rho_b^2\zeta^2.                 \tag{10}
\]

The graph and moving-domain terms cancel one another *after* odd induced
metric continuity is imposed. The final result equals the shell-volume
contribution, but this equality is the result of their explicit cancellation;
it cannot be assumed from the shell term alone.

## 5. Norm identification and remaining scope

The bulk result from the first audit is
\(K_{bulk}=18\mathcal E[P_o]\), where
\(P_o=\rho^2\psi\). Continuity gives
\(\zeta=-P_o(T)/(a_b\rho_b^2)\). Hence (10) becomes

\[
 K_{boundary}=\frac{6|P_o(T)|^2}{a_b\rho_b^2}
        =18c_o|P_o(T)|^2,\qquad c_o=\frac1{3a_b\rho_b^2}>0.
\]

For the outward-normal odd Robin condition
\(p_{SL}P_o'(T)=-c_oP_o(T)\), integration by parts gives

\[
 \mathcal E[P_o]+c_o|P_o(T)|^2
        =\Lambda\int_0^T w_{SL}|P_o|^2dy,
\]

provided the cone Green term vanishes. This proves (1) at eigenroots. The
physical kinetic coefficient is a positive bulk form plus a positive
boundary term for nonzero admissible profiles, and it differs from the odd
SL norm by the same essential \(18\Lambda\) factor as in the even sector.

Time corners have not been discarded by an assumption that a graph is
orthogonal to a time boundary. It generally is not. Instead, use test
perturbations that vanish on a neighborhood of the auxiliary initial/final
time surfaces. The graph then agrees with the background there and its
corner angle has zero perturbation. Corner terms consequently do not
contribute to the local kinetic coefficient calculated above. All time
integrations by parts have zero endpoint variation in this test class.

The calculation does **not** establish the global quotient by admissible
diffeomorphisms, physical symplectic completeness, or a choice of states
across the cone/horizon. Those issues remain as stated in the first audit.
Exceptional harmonic subspaces cannot use \(N=-2\psi\) by the generic
trace-free argument and must continue to be classified separately. This
completion supplies the previously missing local quadratic kinetic
normalization for generic odd modes; it is not a general nonlinear stability
or unrestricted no-ghost theorem.

## 6. Independent verification

`check_odd_kinetic_20261013.py` passed 11/11 exact checks, including two
wrong-formula controls (omitted tangential Einstein flux and omitted
moving-domain term). `ODD_KINETIC_CHECK_RESULTS_20261013.json` records the
expressions and source hash. Run it with the same SymPy 1.14.0 environment
as the first audit.

A separate agent independently derived the covariant graph expression and
the cancellation in Sections 3–4, then reviewed the new script. Its review
found no blocking sign or normalization discrepancy. This is internal
independent verification, not external review. No numerical spectrum or
interval estimate was needed for this action identity.
