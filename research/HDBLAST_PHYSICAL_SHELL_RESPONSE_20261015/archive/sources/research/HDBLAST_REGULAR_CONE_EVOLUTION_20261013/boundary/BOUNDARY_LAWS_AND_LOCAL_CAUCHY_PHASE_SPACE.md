# Regular master boundary laws and a bounded local Cauchy phase space

13 October 2026. Independent analytic continuation of the published cone
checkpoint, confined to this authoring directory. This note proves no
after-contact initial-boundary-value theorem or equality of completed physical
spaces. The bulk reconstruction and its full residual identity are independently
derived in the companion `reconstruction/` note; the local evolution
argument below explicitly assumes that identity.

## 1. Fixed background and regular variable

Fix one nonconstant registered background at positive detuning. In an ambient
Cartesian chart use signature (-++++), `s=R^2-T^2`,
`g=Omega(s)^2 eta`, `phi=Phi(s)`, and `f=Phi_s`. Require that the **actual
Einstein–scalar background**, not an arbitrary smooth extension of its
coefficients, has an analytic continuation on `-s_-<s<s_+` with `Omega>0`
and `f!=0`. Here `Omega(0)=1`, `f(0)=kappa/2!=0`. Analyticity of the
registered regular branch with polynomial potential gives a local candidate
continuation; the background equations continue by the analytic identity
principle. Perturbations and their Cauchy data need only be smooth.

The shell remains the timelike hyperboloid
`R_b(T)=sqrt(r_b^2+T^2)`, with `rho_b=r_b Omega_b`, where
`dr/r=dy/rho`. All shell coefficients below are evaluated on the fixed
background. Define

\[
 F=P/r^4,\quad E=X^A\partial_A,\quad
 L=\Omega_s/\Omega+2f_s/f,\quad A=(\Omega f^2)^{-1}>0.
\]

The exact regular master equation and reconstructed fields are

\[
 \mathcal LF=\Box_\eta F-2L EF-[8L+(8/3)sf^2]F=0,
\]
\[
 h_{AB}=2F(s\eta_{AB}-3X_AX_B),\qquad
 \chi=-\frac3{2\Omega^2f}(E+4)F.                         \tag{1}
\]

The sibling direct-index identity is
`delta E_AB=Omega^-2(3X_A X_B-s eta_AB) L F` and
`delta S=-3(E+6)L F/(2 Omega^4 f)`. It supplies sufficiency for arbitrary
smooth packets without harmonic division. It does not establish completeness.
All coefficients and reconstructed fields in (1) are smooth at the cone and
vertex. Extending master coefficients arbitrarily beyond a causal neighborhood
is a device for Cauchy existence; it does not extend the physical background
equations arbitrarily.

## 2. Exact shell equations: Wentzell and Robin

Let `B_b=p'_b/p_b+sigma''_b/2>0`,
`alpha=(rho^2 p^2)^-1`, `beta=rho_b^4 p_b^2 B_b`,
and `c_o=(3a_b rho_b^2)^-1`, with `a_b>0`. The inherited equations are

\[
 \alpha_b P_y=\beta^{-1}(\Box_\gamma+4)P\quad\text{(even)},
 \qquad \alpha_b P_y=-c_oP\quad\text{(odd)}.
\]

Since `P=r^4 F`, `r_y/r=1/rho` and
`P_y=r^4[nF+(4/rho_b)F]` at the shell, where `n=partial_y`
is the outward physical unit normal. Also `alpha_b beta=rho_b^2 B_b`
and `Box_hb=rho_b^-2 Box_gamma`. Thus the exact regularized laws are

\[
 \boxed{nF-B_b^{-1}\Box_{h_b}F+\gamma_eF=0,
 \qquad \gamma_e=4/\rho_b-4/(\rho_b^2B_b),}              \tag{2}
\]
\[
 \boxed{nF+\gamma_oF=0,
 \qquad \gamma_o=4/\rho_b+p_b^2/(3a_b)>0.}              \tag{3}
\]

The positive even coefficient is an auxiliary Wentzell boundary mass. It
must not be confused with an independently introduced physical shell field.
The previous action calculation has no straight-even-shell kinetic term in
the admitted longitudinal slice.

For the canonical master `psi=sqrt(A)F`,
`n log(A^-1/2)=s_b L_b/rho_b`. Equations (2)–(3) become the same equations
for `psi` with Robin coefficients shifted by `s_b L_b/rho_b`.
In the ambient metric eta, even (2) is

\[
 n_\eta\psi-(\Omega_bB_b)^{-1}\Box_{\Sigma,\eta}\psi
       +(\Omega_b\gamma_e+r_bL_b)\psi=0.               \tag{4}
\]

The canonical boundary kinetic coefficient `(Omega_b B_b)^-1` is positive.
No sign is assumed for its shifted potential term. A bounded negative
zeroth-order term does not by itself obstruct finite-slab energy estimates.

## 3. Global de Sitter clock versus ambient time

At the shell `T=r_b sinh(t)`, `R_b=r_b cosh(t)`. Its exact induced metric is

\[
 h_b=\rho_b^2[-dt^2+\cosh^2t\,d\Omega_3^2]
     =\Omega_b^2[-r_b^2R_b^{-2}dT^2+R_b^2d\Omega_3^2].
\]

For a shell trace `z(T,theta)` the tangential operator is

\[
 \Box_{h_b}z=\rho_b^{-2}
 [-R_b^2 z_{TT}-4Tz_T+(r_b^2/R_b^2)\Delta_{S^3}z],      \tag{5}
\]
\[
 n=(T\partial_T+R_b\partial_R)/\rho_b,
 \qquad n_\eta=(T\partial_T+R_b\partial_R)/r_b.
\]

The `4T` friction coefficient is essential; using `3T` omits the second
derivative of the time conversion. On every finite ambient slab the boundary
lapse `rho_b/R_b` is positive and the shell speed `|T|/R_b<1`. These constants
are not uniform as `|T|` tends to infinity. The global de Sitter foliation only
covers the exterior cone wedge; it is not a Cauchy foliation through the cone.

The bulk equation can be obtained from an auxiliary weighted scalar action
with positive physical kinetic weight `Omega^-4 f^-2`. Adding the auxiliary
boundary action with kinetic coefficient
`Omega_b^-4 f_b^-2/B_b` and mass squared
`B_b gamma_e` gives (2). Odd (3) is supplied by the corresponding Robin
boundary potential. This provides the correct sign for a possible scalar
Wentzell energy method. It is neither the Einstein–scalar quadratic action
nor a proof of the full coupled moving-shell problem after contact.

## 4. Compact-slab existence, energy and shell exclusion

Because `A_s/A=-L`, the canonical field satisfies exactly

\[
 (\Box_\eta-V_c)\psi=0,\qquad
 V_c=3L-2sL_s+sL^2+(8/3)sf^2.                          \tag{6}
\]

Smooth real coefficients and normal hyperbolicity give local smooth Cauchy
solutions and propagation speed one. For a compact causal support and
`m>0`, set

\[
 E_m(T)=\tfrac12\int_{\mathbb R^4}
 (|\psi_T|^2+|\nabla\psi|^2+m^2|\psi|^2)\,d^4X.
\]

If `K_m=sup|m^2-V_c|` on that support, integration by parts gives
`E_m'=Re integral (m^2-V_c) conjugate(psi) psi_T`, hence

\[
 |E_m'|\le(K_m/m)E_m,\qquad
 E_m(T)\le E_m(t_0)\exp[(K_m/m)|T-t_0|].              \tag{7}
\]

Commuted estimates control higher regularity. Since (1) differentiates `F`,
finite reconstructed scalar H1 × L2 data require higher master Sobolev order;
smooth compact data satisfy this automatically. The auxiliary energy is a
control norm, not an asserted gravitational Hamiltonian or ghost criterion.

Suppose initial support is inside `R<=R_0` on `T=t_0`. A transparent sufficient
criterion for a slab `|T-t_0|<=Delta` is

\[
 R_0+\Delta<r_b,\qquad
 (|t_0|+\Delta)^2<s_-,\qquad
 (R_0+\Delta)^2<s_+.                                  \tag{8}
\]

The last bound may be replaced by the known original positive-s background
interval if that covers the support. Negative `s` is controlled by the second
bound. These inequalities keep the full causal support in the actual background
domain and a positive distance from the shell. All reconstructed fields and
derivatives vanish near it. Either parity copying and zero embedding
displacement therefore satisfies every homogeneous linearized junction
condition throughout this slab. Compatibility holds to every order because
the data vanish on an open shell neighborhood. This is an actual short-time
bulk-plus-shell solution, not a statement about what happens after contact.

## 5. Exact actual nonzero-flux pair

The integrated manuscript's exact Cauchy-jet construction is preferable to a
flat-comparator approximation. Choose `t_0>0`, `l>=2`, a normalized real S3
harmonic `H_l`, `lambda=l(l+2)`, and a smooth annular radial cutoff equal to
one near `R=t_0`. Use data
`(F_1,F_1T)=(a zeta H_l,0)` and
`(F_2,F_2T)=(0,b zeta H_l)`, `ab>0`, and satisfy (8).
At `T=R=t_0`, the exact master equation supplies the second time derivatives.
With `C=3/kappa`, the resulting cone scalar jets give exactly

\[
 \chi_1\chi_{2v}-\chi_2\chi_{1v}
       =C^2ab(10+\lambda/2)H_l^2>0
       \quad\text{after angular integration}.           \tag{9}
\]

Every `L(0)` term cancels. The actual one-copy Lee–Wald null current is
`v^3[chi_1 chi_2v-chi_2 chi_1v]/8`; parity doubling gives twice this.
Continuity then supplies nonzero integrated flux on a small cone segment.
The relational scalar `Z=delta[(nabla phi)^2]-(K_g)_s chi/f` independently
excludes pure gauge for each state when `l=2` and `|t_0^2 L(0)|<1`.
Angular degree at least two avoids both exceptional harmonic subspaces.
This is actual registered local solution existence, not a claim that these
states have already been evolved through shell contact.

## 6. Injectivity and a bounded symplectic nondegeneracy statement

The unquotiented differential reconstruction is injective: away from the
vertex the matrix `s eta_AB-3X_A X_B` is nonzero, so `h=0` forces `F=0`.
At the vertex `chi=-6F/f(0)` detects the remaining value. This observation
alone does not prove absence of a gauge kernel or nondegenerate symplectic form.

There is an additional exact local Cauchy result. On `T=0`, take smooth compact
master data `q=F`, `p=F_T` supported in a sufficiently small ball, with zero
fields near the shell. Put `D=r partial_r`, `m_r=9 A r^3/4` (this weight is
unrelated to `m` in (7)). The actual one-copy Einstein–scalar symplectic
pullback, with a consistent future-Cauchy orientation, is

\[
 \sigma((q_1,p_1),(q_2,p_2))=K(q_1,p_2)-K(q_2,p_1),
\]
\[
 K(q,p)=\int dr\,d\Omega_3
 \left[m_r(D+4)q(D+5)p+6r^7\Omega^{-1}qp\right].       \tag{10}
\]

Indeed `chi=-3(D+4)q/(2Omega^2 f)`,
`chi_T=-3(D+5)p/(2Omega^2 f)` on this time-symmetric slice, and
`Psi=r^2q/Omega^2`. Transforming the already derived physical time current
gives measure `Omega^3 r^3`. For real compact data integration by parts yields

\[
 K(q,q)=\int dr\,d\Omega_3\left[
 m_r\{(Dq)^2+(2+9sL)q^2\}+6r^7\Omega^{-1}q^2\right]. \tag{11}
\]

Choose the ball so `2+9sL>=1`. Then `K(q,q)>0` for nonzero `q`.
If a data pair is in the kernel of (10) against all compact test pairs,
testing against `(0,q)` first gives `q=0`, then against `(p,0)` gives `p=0`.
Thus the physical pullback is nondegenerate on this explicit smooth compact
master-data image. Doubling a fixed parity multiplies it by two. Ordinary
smooth zero-charge diffeomorphisms are presymplectic degeneracies on a complete
Cauchy surface; none can produce a nonzero state inside this nondegenerate
image. This does not specify a weaker completed gauge group or identify every
physical state with a master datum. No completed Hilbert-space theorem follows.

The cone trace alone has a different kernel: on a finite positive interval,
`(v partial_v+4)F_0=0` admits `F_0=C(theta)v^-4`. Its scalar cone trace is
zero even when its leading metric trace is nonzero. A cone current therefore
cannot be substituted for the full Cauchy form when testing injectivity.
Smooth extension through the vertex excludes this specific divergent trace,
but a finite-segment analysis alone cannot invoke that exclusion.

## 7. Audit and remaining gate

The analytic hypotheses, support criterion and physical current are essential
parts of these statements. `check_boundary_evolution.py` checks the boundary
conversion, normal/induced metric, canonical transformation, jet cancellation
and Cauchy integration coefficient with deliberate wrong-formula controls.
It imports no prior scientific code and proves no PDE theorem mechanically.
The prospective global boundary problem still needs after-contact evolution,
the complete junction/gauge analysis, and control in a chosen physical
completion. The older M465R1 branch near phi=-1 is not used as a theorem about
the present nonconstant registered branch near phi=+1.
