# Independent audit of physical shell-response nonvanishing

15 October 2026 UTC. Internal critical derivation. This note separates the
physical scalar observable from a master trace and checks the proposed
nonvanishing argument. It does not infer response from existence alone.

## 1. The trace which must activate

For the copy-even longitudinal master put `P=r^4 F` and `Psi=P/rho^2`.
The full coupled reconstruction gives

\[
 \chi_b=-\frac{3P_{y,b}}{\rho_b^2p_b},\qquad
 B_bP_{y,b}=\rho_b^{-2}(\Box_\gamma+4)P_b,
 \quad B_b=p_{y,b}/p_b+\sigma_b''/2>0.
\]

The background shell scalar is the constant `phi_b`. With zero first shell
normal displacement in this even representative, the physical pulled-back
scalar variation is `Phi_1=chi_b`. A nonzero `P_b` alone does not establish
this observable, because `Phi_1` involves its normal derivative.

If `Phi_1` vanished on all of `[0,T1]`, then `P_y,b=0`, and the dynamic
boundary law would give `(Box_gamma+4)P_b=0`. For the selected degree-two
harmonic, its radial coefficient obeys

\[
 P_{b,tt}+3\tanh(t)P_{b,t}+[8\operatorname{sech}^2t-4]P_b=0.
\]

Its initial value and velocity are zero because the datum is supported away
from the shell. Ordinary ODE uniqueness therefore gives `P_b=0` throughout
that interval, as well as `P_y,b=0`.

This implication uses the actual full Einstein–scalar reconstruction and
its physical shell law. It is not a test-scalar substitution.

## 2. Bulk conjugation and an explicit compatible profile

The registered master equation is

\[
 P_{yy}-(2a+2p_y/p)P_y+\rho^{-2}(\Box_\gamma+4)P-\tfrac23p^2P=0.
\]

For `z=log(b/r)`, `dy=-rho dz`, and degree two this becomes

\[
 P_{zz}-P_{tt}-3\tanh(t)P_t+C(z)P_z
 +[4-8\operatorname{sech}^2t-2\rho^2p^2/3]P=0,
 \quad C=\rho(3a+2p_y/p).
\]

Write `P=alpha(z) cosh(t)^(-3/2) w`, where `alpha'/alpha=-C/2`.
Then

\[
 w_{tt}-w_{zz}+\mathcal W(t,z)w=0,
\]
\[
 \mathcal W=\tfrac{35}4\operatorname{sech}^2t-\tfrac{25}4
       +\tfrac12C'+\tfrac14C^2+\tfrac23\rho^2p^2.
\]

Both constants and signs follow directly from the conjugation. In fact
`C=-3(log rho)_z-2(log|p|)_z`, so one can take the positive weight
`alpha=rho^(3/2)|p|`, up to a positive constant.

Choose `0<d`, a width `ell>0`, and a nonnegative smooth bump `q` supported
exactly in `[d,d+ell]`, positive in its interior and nondecreasing on
`[d,d+ell/2]`. The standard exponential bump has these properties. Define

\[
 G(r^2)=r^{-6}\alpha(\log(b/r))q(\log(b/r)),\qquad K=\nu G,
 \quad\nu>0.
\]

For the registered ambient-time master datum
`F_0=r^2G(r^2)Y`, `(F_T)_0=r^2K(r^2)Y`, the normalized initial fields are
`w_0=q` and `w_t(0,z)=nu r(z)q(z)`. The latter uses
`partial_t=r partial_T` at the initial section; forgetting this radius
factor changes the datum. These profiles retain the same nonzero aligned
displacement/velocity matrices and are smooth on the annulus. The coupled
reconstruction supplies the linearized constraints. Coordinate reflections
preserve `Y=n_1^2-n_2^2` and this tensorial datum, so the retained general
Gamma-invariant exact initial-data theorem applies. A profile chosen in
this explicit way must be declared; arbitrary pre-existing profiles are
not identified with it without an argument.

## 3. Quantitative contradiction for the zero-trace hypothesis

Assume the physical trace vanishes on `[0,T1]`, with `T1>d`. Section 1
then gives `P=P_z=0` at the shell. Thus `w=w_z=0` there. Extend `w` by zero
to `z<0` and extend the smooth bounded potential across zero. The extension
satisfies the ordinary whole-line wave equation: there is no distributional
delta or derivative of a delta because both traces vanish. Smoothness can
also be checked recursively from the equation. Only a compact backward
diamond is needed.

Fix `t=d+eta`, where `0<eta<min(ell/2,T1-d)`. Write `Q=q(t)>0`, bound
`|W|<=M` in that diamond, and take `Cvel>=nu r` on the initial leading
segment. Finite propagation says every nonzero source slice in the
backward diamond has `z+s>=d`, while its upper endpoint is `z<=t-s`.
Its spatial width is therefore at most `eta`. Every initial interval
sampled inside that diamond intersects the datum only in `[d,t]`; monotonicity
bounds the displacement by `Q` and the velocity by `Cvel Q`.

Let `S` be the supremum of the absolute extended solution in the diamond.
Duhamel's formula and its factor one half give

\[
 S\le Q(1+C_{vel}\eta/2)+\beta S,\qquad
 \beta=Mt\eta/2.
\]

The free solution at the observation point has displacement contribution
`Q/2` and nonnegative velocity contribution. Hence, if `beta<1`,

\[
 w(t,0)\ge\frac Q2-
  \beta Q\frac{1+C_{vel}\eta/2}{1-\beta}
 =\frac{Q[1-\beta(3+C_{vel}\eta)]}{2(1-\beta)}.
\]

For sufficiently small positive `eta`,
`Mt eta(3+Cvel eta)<2`, so this lower bound is strictly positive,
contradicting the assumed zero extension. The bump can be flat to every
order at `d`; only `q(d+eta)>0` and its monotonicity are used. No nonzero
value or derivative at the smooth front is assumed.

This establishes physical `Phi_1` activation somewhere in every interval
just after `d`. The positive estimate belongs to the contradiction solution
under the zero-trace hypothesis. It is not a claimed sign or numerical
lower bound for the actual boundary waveform.

For pure displacement or pure velocity one can alternatively use the
registered time-parity/sideways-wave uniqueness argument. That would invoke
a different tangent if the previous state had nonzero proportional velocity;
the positive-profile argument above retains `K=nu G`.

## 4. Comparing the exact-family tangent with the master solution

No unproved uniqueness in a geometric gauge quotient is needed. From the
longitudinal reconstructed solution solve

\[
 (\Box_0+\operatorname{Ric}_0)\xi=-C[h_L],\qquad
 \xi_y|_b=0,\quad(\partial_y-2a_b)\xi_a|_b=0,
\]

with `xi=xi_t=0` initially. The established full linear mixed solver supplies
this boundary-preserving gauge change. It gives
`h_W=h_L+Lie_xi g0`, `chi_W=chi_L+xi^y p`, with `C[h_W]=0` and
`h_W,ya=0`. The physical shell laws remain satisfied. At the shell
`chi_W=chi_L` because `xi^y=0`; equivalently the shell pullback is geometric.
Canonical initial fields and initial metric values agree, while the remaining
five time velocities are precisely the harmonic completion. The exact-family
initial lapse/shift values must match these values. Fixed-gauge linear
uniqueness then identifies its derivative with this harmonic solution.

## 5. A physical clock and the exact nonnegative response

Coordinate-time slices alone do not make a finite-parameter observable
intrinsic. Use future unit-speed timelike geodesics of the induced shell
metric, launched orthogonally from the fixed initial shell three-sphere.
They define slices `S_epsilon(tau)` on a common finite proper-time interval.
At the background the flow is `(tau,Omega)->(t=tau/rho_b,Omega)`, with
nonvanishing spatial Jacobian proportional to `cosh(t)^3`. On a fixed finite
interval, sufficiently small C2 perturbations preserve the smooth
nonintersecting congruence and spacelike slices. Choose a strict time margin
between the observation window and the end of the controlled slab.

Let `Phi_epsilon` be the shell scalar pulled back by this intrinsic clock
flow and `mu_epsilon` the induced slice volume form. Define

\[
 Q_\varepsilon(\tau)=\int_{S_\varepsilon(\tau)}
                   (\Phi_\varepsilon-\phi_b)^2\,d\mu_\varepsilon.
\]

The reference initial section is part of the observable definition. Since
`Phi_0=phi_b` is constant on the whole background shell, its first pulled-back
variation is `Phi_1=chi_b(tau/rho_b,Omega)`. First clock shifts and slice
measure variations do not change the quadratic coefficient. Therefore

\[
 Q_\varepsilon(\tau)=\varepsilon^2 C(\tau)+o(\varepsilon^2),\quad
 C(\tau)=\int\Phi_1^2\,d\mu_0,\quad
 \left.\partial_\varepsilon^2Q_\varepsilon(\tau)\right|_0=2C(\tau).
\]

For `Y=n_1^2-n_2^2`, `integral_S3 Y^2 dmu_kappa=pi^2/3`, so if
`chi_b=v_b(t)Y`,

\[
 C(\rho_b t)=\frac{\pi^2}{3}\rho_b^3\cosh^3(t)v_b(t)^2
 =\frac{3\pi^2\cosh^3(t)}{\rho_b p_b^2}P_{y,b}(t)^2.
\]

A more robust statement uses the fixed clock window

\[
 R_\varepsilon=\int_0^{\rho_bT_1}Q_\varepsilon(\tau)\,d\tau.
\]

Activation and continuity imply `Cwin=integral C(tau)dtau>0`, whence
`R_epsilon=epsilon^2 Cwin+o(epsilon^2)` and
`R_epsilon >= (Cwin/2) epsilon^2` for sufficiently small nonzero parameters.
This lower bound is analytic and existential unless its coefficient and
remainder are separately computed.

## 6. Scope traps retained in the final claim

Exact localized constraint corrections may occupy a larger annulus than the
first-order packet. They can generate order-epsilon-squared fields before
the first-order front. Thus `d` is the linear activation infimum, not an
asserted exact first-response time; the exact squared response can be order
epsilon-to-the-fourth before it. The window coefficient avoids this issue.

The result concerns one deliberately specified coupled, constraint-compatible
profile under the retained exact-family/background hypotheses. It does not
prove nonvanishing for every profile, define a physical gravitational energy,
show observational significance, establish external novelty or demonstrate
a higher-dimensional origin of the Big Bang.
