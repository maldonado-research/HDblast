# Full quadratic TT bulk source in the inherited scalar gauge

Date: 2026-10-14. Scope: a bulk differential identity and its forced tensor
equation. No shell-contact existence or nonlinear boundary estimate is proved.

## 1. Actual inherited ansatz and conventions

The input is the spacetime reconstruction in
`research/HDBLAST_REGULAR_CONE_EVOLUTION_20261013/manuscript/REGULAR_CONE_EVOLUTION.md`,
equations (3), (6), and (9), read from
`/workspace/hdblast-research-work/publish-HDblast-20261014/`:

\[
 g_0=dy^2+\rho(y)^2\gamma,\quad
 \gamma=-d\tau^2+\cosh^2\tau\,\sigma_{ab}dn^a dn^b,
 \quad \operatorname{Ric}(\gamma)=3\gamma,
\]
\[
 h_{1,yy}=-4\Psi,\qquad h_{1,y\mu}=0,\qquad
 h_{1,\mu\nu}=2\rho^2\Psi\gamma_{\mu\nu},\qquad
 \chi=-\frac{3(\rho^2\Psi)_y}{\rho^2p},\quad p=\phi_{0,y}.
\]

Use the trace-reversed equations in the inherited normalization,

\[
 \mathscr E_{AB}=R_{AB}-\phi_{,A}\phi_{,B}-\frac23U(\phi)g_{AB}=0.
\]

Write the parameter expansion as

\[
 g_\epsilon=g_0+\epsilon h_1+\frac{\epsilon^2}{2}k_2+O(\epsilon^3),
 \qquad \phi_\epsilon=\phi_0+\epsilon\chi
                  +\frac{\epsilon^2}{2}\chi_2+O(\epsilon^3).
\]

Thus \(D\mathscr E[k_2,\chi_2]+D^2\mathscr E[(h_1,\chi)^2]=0\).
The symbols \([\epsilon^2]\) below refer to the coefficient in the
**linear path** \((g_0+\epsilon h_1,\phi_0+\epsilon\chi)\), so
\(D^2\mathscr E=2[\epsilon^2]\mathscr E\).

\(\Pi_{\rm TT}\) always means the spatial tensor projection on the
closed unit \(S^3\), at fixed \((y,\tau)\). It does not mean a Lorentzian
four-dimensional TT projection on de Sitter.

## 2. Exact geometric reduction, including the radial lapse

For the linear metric path put

\[
 N^2=1-4\epsilon\Psi,\qquad q_{\mu\nu}=\rho^2e^{2\omega}\gamma_{\mu\nu},
 \qquad \omega=\tfrac12\log(1+2\epsilon\Psi),\quad
 \nu=\log N=\tfrac12\log(1-4\epsilon\Psi).
\]

There is no radial shift. The radial extrinsic curvature satisfies

\[
 K_{\mu\nu}=\frac{1}{2N}\partial_y q_{\mu\nu}
            =\frac{a+\omega_y}{N}q_{\mu\nu},\qquad a=\rho_y/\rho.
\]

Every radial-extrinsic-curvature term in \(R_{\mu\nu}\) is therefore
proportional to \(\gamma_{\mu\nu}\). The only non-trace terms are
\(R[q]_{\mu\nu}-N^{-1}D[q]_\mu D[q]_\nu N\). With derivatives on
\(\gamma\), their exact expression modulo \(\gamma\)-trace is

\[
 R_{\mu\nu}\equiv
 -2\omega_{;\mu\nu}-\nu_{;\mu\nu}
 +2\omega_{,\mu}\omega_{,\nu}-\nu_{,\mu}\nu_{,\nu}
 +\omega_{,\mu}\nu_{,\nu}+\omega_{,\nu}\nu_{,\mu}.
 \tag{1}
\]

Since
\(\omega=\epsilon\Psi-\epsilon^2\Psi^2+O(\epsilon^3)\) and
\(\nu=-2\epsilon\Psi-4\epsilon^2\Psi^2+O(\epsilon^3)\),

\[
 [\epsilon^2]R_{\mu\nu}
 \equiv 6(\Psi^2)_{;\mu\nu}-6\Psi_{,\mu}\Psi_{,\nu}
 =12\Psi\Psi_{;\mu\nu}+6\Psi_{,\mu}\Psi_{,\nu}.
 \tag{2}
\]

For angular indices, a \(\gamma\)-Hessian differs from the unit-sphere
Hessian only by a spatial trace term. The sphere TT projection annihilates
both Hessians and trace terms. Also, \(\phi_0\) has no angular derivative;
the only quadratic anisotropic matter term is
\(-\chi_{,a}\chi_{,b}\). The potential term is pure spatial trace.
Consequently the **full** geometric-plus-matter identity is

\[
 \boxed{
 \Pi_{\rm TT}D^2\mathscr E_{ab}[(h_1,\chi)^2]
 =-12\Pi_{\rm TT}(D_a\Psi D_b\Psi)
   -2\Pi_{\rm TT}(D_a\chi D_b\chi).}
 \tag{3}
\]

This quadratic identity needs no first-order PDE, reconstruction
substitution, or background field equation. Those are needed to make
\((h_1,\chi)\) an actual first-order solution and to simplify the linear
tensor operator, respectively. All radial and time derivatives in (1)
have been retained: their spatial anisotropy is either encoded in the
Hessian or is pure trace.

## 3. The specified quadrupole: exact absence of cancellation

Let \(B=\operatorname{diag}(1,-1,0,0)\),
\(Y(n)=n^TBn\), \(\Psi=u(y,\tau)Y\), and
\(\chi=v(y,\tau)Y\). Let \(\mathcal T_{ab}\) be the sphere tensor
induced by the Weyl part of \(B\wedge B\), with sectional components

\[
 C_{1212}=C_{3434}=-1/3,\qquad
 C_{1313}=C_{1414}=C_{2323}=C_{2424}=1/6.
\]

The exact sphere identity supplied and directly proved in the existing
`PRECONTACT_SECOND_ORDER_DYNAMICS.md`, equation (12a), is

\[
 (DY\otimes DY)^{\rm TF}
  =\tfrac12(\operatorname{Hess}F_2)^{\rm TF}
   +\tfrac13(\operatorname{Hess}F_4)^{\rm TF}-\tfrac43\mathcal T,
\]

where \(F_2=n_1^2+n_2^2-1/2\) and
\(F_4=Y^2-(n_1^2+n_2^2)/2+1/12\). Hence
\(\Pi_{\rm TT}(DY\otimes DY)=-4\mathcal T/3\), and (3) gives

\[
 \boxed{
 \Pi_{\rm TT}D^2\mathscr E_{ab}
       =\left(16u^2+\frac83v^2\right)\mathcal T_{ab},\qquad
 \Pi_{\rm TT}D\mathscr E_{ab}[k_2,\chi_2]
       =-\left(16u^2+\frac83v^2\right)\mathcal T_{ab}.}
 \tag{4}
\]

For real perturbations the coefficient is strictly positive wherever
\((u,v)\ne(0,0)\). The geometric and matter contributions reinforce
each other. A zero total forcing by cancellation is therefore impossible
for this quadrupole in the inherited first-order longitudinal gauge.
This is not a claim of a nonzero source on every time slice: for pure
velocity master data with \(F|_{T=0}=0\), formula (8) gives
\(u=v=0\) initially, even when \(F_T\ne0\).

The stated normalization is \(\langle|\mathcal T|^2\rangle_{S^3}=1/12\).
Thus the averaged contraction of \(D^2\mathscr E\) with \(\mathcal T\)
is \(4u^2/3+2v^2/9\). This is the residual, whereas the right-hand side
of the second-order linear equation has its opposite sign.

## 4. Explicit forced tensor equation

Define the tensor amplitude of the second derivative by

\[
 (k_2)^{\rm TT}_{ab}=\rho^2\cosh^2\tau\,H(y,\tau)\mathcal T_{ab}.
 \tag{5}
\]

The \(l=2\) tensor harmonic obeys
\(\Delta_{S^3}\mathcal T=-6\mathcal T\). The linear trace-reversed
Einstein equation in this sector is

\[
 \Pi_{\rm TT}D\mathscr E_{ab}[k_2,\chi_2]
 =-\frac12\rho^2\cosh^2\tau\,\mathfrak L_T H\,\mathcal T_{ab},
\]
\[
 \mathfrak L_T H=
 H_{yy}+4aH_y-\rho^{-2}
  \left(H_{\tau\tau}+3\tanh\tau H_\tau+8\operatorname{sech}^2\tau H\right).
 \tag{6}
\]

The curvature contribution converts the tensor-Laplacian eigenvalue
\(-6\) into \(-8\) in this wave operator. From (4),

\[
 \boxed{\mathfrak L_T H
 =\frac{16}{3\rho^2\cosh^2\tau}\left(6u^2+v^2\right).}
 \tag{7}
\]

This equation includes all scalar-scalar bulk forcing in this tensor
channel. Other second-order scalar and vector amplitudes do not enter
its linear left side, by the spherical background symmetry.

For an ambient form that continues across the cone away from the
polar origin, let

\[
 T=r\sinh\tau,\quad R=r\cosh\tau,\quad s=R^2-T^2=r^2,
 \quad dy=\Omega(s)dr,\quad \rho=r\Omega(s).
\]

Write the inherited regular master as \(F=\mathcal U(T,R)Y\). Then

\[
 u=\frac{s\mathcal U}{\Omega^2},\qquad
 v=-\frac{3}{2\Omega^2 f}(T\partial_T+R\partial_R+4)\mathcal U,
 \qquad f=\Phi_s.
 \tag{8}
\]

With the same amplitude \(H(T,R)\), equation (7) becomes

\[
 \boxed{
 \left[-\partial_T^2+\partial_R^2+\frac3R\partial_R
 +6\frac{\Omega_s}{\Omega}(T\partial_T+R\partial_R)-\frac8{R^2}\right]H
 =\frac{16}{3R^2}(6u^2+v^2).}
 \tag{9}
\]

Equivalently its right side is

\[
 \frac{16}{3R^2\Omega^4}
 \left[6s^2\mathcal U^2+
 \frac9{4f^2}\big((T\partial_T+R\partial_R+4)\mathcal U\big)^2\right].
 \tag{10}
\]

The derivation using radial de Sitter slices initially applies on
\(s>0\). Equation (9), together with the ambient reconstruction, has
smooth coefficients across \(s=0\) for \(R>0\), \(\Omega>0\), and
\(f\ne0\), so the differential identity continues there. The ordinary
\(R=0\) polar-coordinate regularity conditions still have to be imposed;
equation (9) is not a claim that arbitrary radial amplitudes are regular
at the origin.

For smooth quadrupole data the apparent origin poles can be removed
explicitly. Write \(z=R^2\), \(\mathcal U=R^2U(T,z)\), and
\(H=R^2J(T,z)\), with \(U,J\) smooth at \(z=0\). Then the ambient
spatial tensor is
\((k_2)^{\rm TT}_{ij}=\Omega^2J\,C_{ikjl}X^kX^l\), manifestly
smooth. Equation (9) becomes

\[
 \begin{split}
 &-J_{TT}+4zJ_{zz}+16J_z
 +6\frac{\Omega_s}{\Omega}(TJ_T+2zJ_z+2J)\\
 &\qquad=\frac{16}{3\Omega^4}
 \left[6s^2U^2+\frac9{4f^2}(TU_T+2zU_z+6U)^2\right],
 \qquad s=z-T^2.
 \end{split}\tag{11}
\]

This removes all explicit inverse powers of \(R\); the \(z\) equation
has the usual radial-coordinate degeneracy and does not replace a
smoothness or existence argument.

If one instead denotes the coefficient of \(\epsilon^2\) in the metric
by \(\widehat k_2\), then \(\widehat k_2=k_2/2\); its amplitude and the
right side of (7)/(9) are both one half of those displayed above.

## 5. Gauge scope and remaining contact problem

The calculation fixes the inherited first-order longitudinal
reconstruction. An order-two coordinate change alone adds
\(\mathcal L_{\xi_2}g_0\) to \(k_2\), whose spatial TT part vanishes.
A change of **first-order** gauge also changes the second-order tensor
by the TT part of
\(2\mathcal L_{\xi_1}h_1+\mathcal L_{\xi_1}^2g_0\), and transforms the
quadratic source consistently. Thus (4) is a complete source statement
in the specified first-order gauge, not an assertion that every raw
second-order tensor amplitude or quadratic source is gauge invariant.

The next contact problem must combine (7) or (9) with the complete
second-order tensor junction conditions, including the shell embedding
and the prescribed matter law, and with compatible second-order Cauchy
data. This note supplies neither those junction conditions nor an
energy estimate or a shell-contact existence theorem. A nonzero source
also does not by itself prove a nonzero outgoing radiation observable.

## 6. Independent checks

`verify_quadratic_tt.py` computes the full five-dimensional Christoffel
expansion along the linear scalar metric path in a flat de Sitter
chart. It checks all six off-diagonal tangential components and all
three independent diagonal differences against (2). It also checks
the general tangential-dimension formula, the Taylor factor, and the
exact radial-to-ambient transformation of the tensor operator. It
imports no inherited scientific code. Its stated wrong-formula
controls detect deletion of the geometric source, reversal of the
matter sign, a missing Taylor factor, a wrong tensor eigenvalue, and a
wrong conformal drift.

`verify_tensor_operator.py` independently uses a global de Sitter chart
with Hopf coordinates on \(S^3\). In this chart the actual tensor above
is diagonal:

\[
 \mathcal T_{\theta\theta}=1/6,\quad
 \mathcal T_{\varphi\varphi}=\sin^2\theta(1/6-\sin^2\theta/2),\quad
 \mathcal T_{\psi\psi}=\cos^2\theta(1/6-\cos^2\theta/2).
\]

It checks the trace and divergence, then directly linearizes the full
five-dimensional Ricci tensor and checks all 15 independent
trace-reversed Einstein components against (6). The only background
substitution in this operator check is
\(\tfrac23U=3/\rho^2-\rho_{yy}/\rho-3a^2\).

These are exact symbolic identity checks, not a numerical evolution,
and they do not verify a shell boundary condition.

`verify_sphere_projection.py` separately recomputes the complete
Hessian decomposition in Hopf coordinates, checks every independent
off-diagonal component and diagonal trace difference, and verifies the
normalization and projected coefficient by exact sphere integration.
With \(x=\sin^2\theta\), the normalized Hopf measure makes \(x\)
uniform on \([0,1]\); the azimuthally averaged stress contraction is
\(2x/3-4x^2/3\), and the tensor norm is
\(1/36+(1/6-x/2)^2+(1/6-(1-x)/2)^2\). Their integrals are
\(-1/9\) and \(1/12\), respectively.

`verify_regular_amplitude.py` checks the two differential substitutions
used in (11). The pinned runtime is Python 3.12.14 with SymPy 1.14.0.
All four checkers pass under ordinary Python and under `python -O`.
The deliberately false matter-only claim fails with exit code 1 in both
modes. `VERIFICATION_RECEIPT.json` records counts, versions, commands,
and SHA-256 hashes of these local artifacts and the inherited inputs.
