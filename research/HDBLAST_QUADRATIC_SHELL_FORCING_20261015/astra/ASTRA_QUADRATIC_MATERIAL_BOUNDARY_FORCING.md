# Known quadratic material shell forcing from the physical first variation

15 October 2026 UTC. Independent derivative and exact-moment derivation for
the registered background and the width ell = 10^-6 first tangent.

The new object is the known constitutive part of the second differentiated
shell equations. The previous conditional construction of a complete
finite-angular second coefficient remains in place. This note neither
replaces that construction nor evaluates its unknown second scalar/metric
response. All prior files used are identified by bytes in
`ASTRA_INPUT_PINS.json`; no prior scientific program was run or imported.

## 1. Fixed conventions and geometric variables

The shell has intrinsic spacetime dimension four. On both interiors the
normal points outward toward the shell. There is one common tension, with

\[
 \sigma(\phi)=2\left(\frac13+\eta^2+\frac{\eta^3}3\right)
       +\delta(1+c+c\eta),\qquad \eta=\phi-1,
\]
\[
 \delta=10^{-3},\quad c=96714321425692/161839258930731,
 \quad \sigma''_b=4+4\eta_b,\quad \sigma'''_b=4.
\tag{1}
\]

Use any specified smooth material identification of the physical shell.
Let \(\psi_\epsilon\) be the pulled-back scalar and
\(p_{i,\epsilon}\) the pulled-back scalar \(n_i\phi_\epsilon\), where
the unit normal, embedding and evaluation all vary. For every field use
derivative coefficients

\[
 q_\epsilon=q_0+\epsilon q_1+\tfrac12\epsilon^2q_2+o(\epsilon^2).
\tag{2}
\]

The differentiation can be performed at the level of formal second jets;
it does not assume or establish a new exact solution family.

Thus \(p_{i,2}\) is not an Eulerian partial derivative of the second bulk
scalar alone. It already contains the normal, embedding and quadratic
pullback terms. On the background shell \(y=y_b\),

\[
 \psi_0=\phi_b,\quad p_{+,0}=p_{-,0}=p_b=-\sigma'_b/2,
 \quad (K_i)^a{}_b=(\sigma_b/6)\delta^a_b.
\tag{3}
\]

All these background scalar or mixed-identity coefficients are constant
over the shell, including its initial section. This follows from
\(\phi_0=\phi_0(y)\) and \(y=y_b\); it is stronger than constancy on
only the initial section. The known first physical scalar in the proper
clock is

\[
 \psi_1(\tau,\Omega)=S(t)Y(\Omega),\quad
 t=\tau/\rho_b,\quad Y=n_1^2-n_2^2,
\]
\[
 S(\ell x)=\frac{\mathcal G}{\ell}\cosh^{-3/2}(\ell x)H(x),
 \qquad \mathcal G=\frac3{\rho_b^3p_b},\quad x=t/\ell.
\tag{4}
\]

## 2. The outward-normal sum and the per-copy factor

Differentiate the physical scalar junction itself:

\[
 p_{+,\epsilon}+p_{-,\epsilon}+\sigma'(\psi_\epsilon)=0.
\]

Its first and second rows are

\[
 p_{+,1}+p_{-,1}+\sigma''_b\psi_1=0,
\]
\[
 \boxed{p_{+,2}+p_{-,2}+\sigma''_b\psi_2
              =f_{\Sigma p}:=-4\psi_1^2.}
\tag{5}
\]

Both normals have the same outward convention. A difference of normal
traces would be a different equation. For the retained identical-copy
family only, division of (5) by two gives

\[
 \boxed{p_{i,2}+\frac{\sigma''_b}{2}\psi_2
              =f_p:=-2\psi_1^2.}
\tag{6}
\]

The sum forcing is \(-4\psi_1^2\); the per-copy forcing is
\(-2\psi_1^2\). There is no extra copy of the tension. In an expansion
without the factorial in (2), the coefficient of \(\epsilon^2\) in
the sum junction is half its second derivative. Those conventions cannot
be mixed. Equations (5)–(6) do not set \(\psi_2=0\).

## 3. Tension, stress and extrinsic trace

The scalar tension has the invariant second derivative

\[
 \sigma_2=\sigma'_b\psi_2+q_\sigma,
 \qquad \boxed{q_\sigma=\sigma''_b\psi_1^2.}
\tag{7}
\]

Define the shell stress by \(T_{ab}=-\sigma h_{ab}\). Its mixed form
is exactly \(T^a{}_b=-\sigma\delta^a_b\). Consequently

\[
 (T^a{}_b)_2=-(\sigma'_b\psi_2+q_\sigma)\delta^a_b,
 \quad (\operatorname{tr}T)_2=-4\sigma'_b\psi_2-4q_\sigma.
\tag{8}
\]

The known stress-trace constitutive term is
\(q_{\operatorname{tr}T}=-4\sigma''_b\psi_1^2\). The sign of that
term does not determine the sign of the full second stress trace. In
particular, raising an index with only the background inverse would fail
to compute the invariant mixed trace. The lower-index stress row is

\[
 T_{2,ab}=-\sigma_bh_{2,ab}-2\sigma'_b\psi_1h_{1,ab}
 -(\sigma'_b\psi_2+\sigma''_b\psi_1^2)h_{0,ab}.
\tag{9}
\]

For completeness, the copy-even Israel law gives

\[
 K_{i,2,ab}-\frac{\sigma_b}6h_{2,ab}
 -\frac{\sigma'_b}6\psi_2h_{0,ab}
 =\frac{\sigma'_b}3\psi_1h_{1,ab}
  +\frac{\sigma''_b}6\psi_1^2h_{0,ab}.
\tag{10}
\]

Here \(K_{i,2,ab}\) is the full material geometric derivative. The
mixed extrinsic trace obeys
\((\operatorname{tr}K_i)_2=(2/3)(\sigma'_b\psi_2+q_\sigma)\).
Its known term is \(2q_\sigma/3\) per copy and \(4q_\sigma/3\)
for the two-copy sum. The metric term in (10), or the corresponding
Eulerian quadratic normal terms, must still be supplied when a full
component boundary source is wanted. Squaring (4) alone does not provide
all such terms.

## 4. Material identification and intrinsic proper-clock consistency

Let a parameter-dependent shell reidentification have first generator
\(v\) and second generator \(w\). For a pulled-back field \(q\),

\[
 \widehat q_1=q_1+\mathcal L_vq_0,\qquad
 \widehat q_2=q_2+2\mathcal L_vq_1
            +\mathcal L_v^2q_0+\mathcal L_wq_0.
\tag{11}
\]

This includes the retained intrinsic geodesic proper-clock map, on its
specified compact slab. Since (3) is constant,

\[
 \widehat\psi_1=\psi_1,\quad \widehat p_{i,1}=p_{i,1},\quad
 \widehat\psi_2=\psi_2+2v(\psi_1),\quad
 \widehat p_{i,2}=p_{i,2}+2v(p_{i,1}).
\tag{12}
\]

The left side of (5) changes by
\(2v(p_{+,1}+p_{-,1}+\sigma''_b\psi_1)=0\), using the first
junction and constancy of \(\sigma''_b\). Its known right side is
unchanged. The per-copy statement follows from its first row. Likewise
\(\widehat\sigma_2-\sigma_2=2\sigma'_bv(\psi_1)\), exactly the
change in the response term of (7); \(q_\sigma\) is unchanged.
Because \(\mathcal L_v\operatorname{Id}=0\), the mixed-stress
statement (8) has the same property. This proves consistency with the
proper clock without suppressing the second clock transport.

The first embedding contribution remains essential:
\(\psi_1=\chi_b+p_b\zeta_1\) before a straight-shell representative
is chosen. Constancy of the background shell trace removes tangential
clock/label contributions, not the normal embedding contribution.

## 5. Exact scalar harmonic projections on the unit three-sphere

Set \(J=n_1^2+n_2^2\), \(F_2=J-1/2\), and
\(F_4=Y^2-J/2+1/12\). Then

\[
 \boxed{Y^2=\frac16+\frac12F_2+F_4,}\qquad
 \Delta_{S^3}F_2=-8F_2,\quad \Delta_{S^3}F_4=-24F_4.
\tag{13}
\]

For an independent proof, use coordinates \(X\in\mathbb R^4\) and
\(r^2=|X|^2\). The homogeneous polynomials
\(X_1^2+X_2^2-r^2/2\) and
\((X_1^2-X_2^2)^2-r^2(X_1^2+X_2^2)/2+r^4/12\) are Euclidean
harmonic. The polar-coordinate Laplacian gives the stated eigenvalues.
Distinct eigenvalues yield orthogonality. No degree-two projection
parallel to the original \(Y\) is present: \(F_2\), not \(Y\),
is the degree-two component of this square.

With \(\operatorname{Vol}(S^3)=2\pi^2\), exact sphere moments give

\[
 \int Y^2=\frac{\pi^2}{3},\quad
 \int Y^4=\frac{3\pi^2}{20},\quad
 \|F_2\|_2^2=\frac{\pi^2}{6},\quad
 \|F_4\|_2^2=\frac{19\pi^2}{360}.
\tag{14}
\]

One direct moment formula, also used by the new checker, is
\(\mathbb E\prod_{j=1}^4 n_j^{2a_j}
=\prod_j(1/2)_{a_j}/(2)_{\sum a_j}\); odd moments vanish.

For any known source \(q_a=aS(t)^2Y^2\), with constant \(a\),
the projection coefficients in the bases \((1,F_2,F_4)\) are
\((aS^2/6,aS^2/2,aS^2)\). Their squared angular norms are

| Degree | Squared norm of the projection of \(q_a\) |
| --- | --- |
| 0 | \(a^2S^4\pi^2/18\) |
| 2 | \(a^2S^4\pi^2/24\) |
| 4 | \(19a^2S^4\pi^2/360\) |

The sum is \(3a^2S^4\pi^2/20\). Only the degree-zero part
survives an angular integral, although the source is not purely isotropic.
Here take \(a=-4,-2,\sigma''_b,-4\sigma''_b\) for (5), (6),
(7), and the known stress trace, respectively.

## 6. Source profile and integrated identities

The common nonnegative source profile is exactly

\[
 \psi_1^2=\frac{\mathcal G^2}{\ell^2}
       \cosh^{-3}(\ell x)H(x)^2Y^2.
\tag{15}
\]

It vanishes for \(0\le x\le1\), because the established first
variation does. This is a statement about the known quadratic material
source; it fixes neither a full second-order arrival time nor the exact
nonlinear first arrival. Zeros or a sign reversal of \(H\) do not
reverse the sign of a squared constitutive source.

Use the background proper-volume measure on the fixed proper-clock slab:

\[
 d\nu_0=d\tau\,d\mu_{c_0(\tau)}
 =\rho_b^4\cosh^3t\,dt\,d\Omega.
\]

For the response coefficient already certified for the same tangent,

\[
 \mathcal C_{3\ell}=\int\psi_1^2d\nu_0
 =\frac{\pi^2\rho_b^4}{3}\int_0^{3\ell}\cosh^3tS(t)^2dt
 =\frac{\mathcal K}{\ell}\int_0^3H^2dx,
 \quad \mathcal K=\frac{3\pi^2}{\rho_b^2p_b^2}.
\tag{16}
\]

Hence the exact known-source identities are

\[
 \boxed{\int f_{\Sigma p}d\nu_0=-4\mathcal C_{3\ell},\quad
 \int f_pd\nu_0=-2\mathcal C_{3\ell},\quad
 \int q_\sigma d\nu_0=\sigma''_b\mathcal C_{3\ell},\quad
 \int q_{\operatorname{tr}T}d\nu_0=-4\sigma''_b\mathcal C_{3\ell}.}
\tag{17}
\]

Their absolute integrals have the corresponding positive constants,
since \(\sigma''_b>0\). These are signed integrals, or equivalently
absolute-source integrals with the sign restored. A spacetime squared
\(L^2\) source norm would involve \(S^4\), not just \(\mathcal C\).

The retained bounds \(-9\,10^{-5}\le\eta_b\le-8\,10^{-5}\) and
\(6.06\,10^{11}<\mathcal C_{3\ell}<7.44\,10^{11}\) give the
following exact-decimal outward corollaries:

\[
 -2.976\,10^{12}<\int f_{\Sigma p}d\nu_0<-2.424\,10^{12},
\]
\[
 2.42378184\,10^{12}<\int q_\sigma d\nu_0<2.97576192\,10^{12},
\]
\[
 -1.190304768\,10^{13}<\int q_{\operatorname{tr}T}d\nu_0
 <-9.69512736\,10^{12}.
\tag{18}
\]

These bounds inherit the prior analytic and arithmetic hypotheses.

They are not second variations of an integral of the entire tension.
For example, if the physical volume density on this fixed clock domain is
\(d\nu_\epsilon=(1+\epsilon m_1+\epsilon^2m_2/2+o(\epsilon^2))d\nu_0\),
then

\[
 \left(\int\sigma(\psi_\epsilon)d\nu_\epsilon\right)_2
 =\int\left[\sigma'_b\psi_2+q_\sigma
       +2\sigma'_b\psi_1m_1+\sigma_bm_2\right]d\nu_0.
\tag{19}
\]

Only the squared-difference response used to define \(\mathcal C\)
loses these extra contributions at quadratic order, because its scalar
integrand and first derivative vanish at the background.

## 7. Verification and scope

`check_astra_material_forcing.py` performs fresh symbolic differentiation
of the registered cubic tension, the scalar sum, the stress product, the
inverse-metric cancellation, clock transport, harmonic polynomials, sphere
moments and the response prefactor. It uses explicit failures, so checks
survive `python -O`. Controls reverse one normal, duplicate the tension,
halve the Hessian term, corrupt a projection, omit the unknown response,
drop the second clock transport, and drop the metric-stress product term.
The receipts record actual subprocess exit codes in both Python modes.
`run_astra_material_controls.py` replays the baseline and all seven
corruptions, writes the individual JSON outputs and a consolidated receipt,
and fails if a baseline is rejected or a corruption is accepted.

This is a boundary forcing certificate for the known first variation.
It gives no evaluated complete second-order solution, backreaction sign,
new exact-family existence theorem, effective finite-amplitude radius,
nonlinear remainder estimate or Big Bang evidence. The earlier conditional
complete finite-angular coefficient construction is preserved with its
own stated input hypotheses and domain. These symbolic checks validate
specified formulas, not an entire analytic proof or external peer review.
