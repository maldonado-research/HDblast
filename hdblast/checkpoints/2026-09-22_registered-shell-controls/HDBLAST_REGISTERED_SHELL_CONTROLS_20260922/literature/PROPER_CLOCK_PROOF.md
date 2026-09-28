# Exact local proper-clock coordinate check

22 September 2026. This note verifies coordinate identities and the registered pure-tension junction conditions. It does not evolve the Einstein–scalar system, prove global chart coverage, or continue the numerical solution through its late-time coordinate failure.

## Premises and notation

Start with

\[
ds^2=e^{2B(t,z)}(-dt^2+dz^2)+e^{2A(t,z)}d\mathbf x^2,
\qquad \phi=\phi(t,z),\qquad z_b=0.
\]

Take the normal pointing toward increasing \(z\): \(n=e^{-B}\partial_z\). The supplied Chat 13 equation derivation specifies this orientation with bulk on \(z<0\) and the boundary conditions

\[
e^{-B}A_z=e^{-B}B_z=\sigma(\phi_b)/6,
\qquad e^{-B}\phi_z=-\sigma_{,\phi}(\phi_b)/2.
\]

Use one common \(C^3\) function \(F\) on both null coordinates, with finite strictly positive \(F'\) throughout the mapped neighborhood. Assume the original fields are sufficiently smooth for the derivatives used. These assumptions give a local orientation-preserving diffeomorphism. The condition must hold off the brane as well as on it.

Write

\[
u=t-z,\quad v=t+z,\quad U=F(u),\quad V=F(v),
\quad T=(U+V)/2,\quad Z=(V-U)/2.
\]

Define \(p=F'(u)>0\), \(q=F'(v)>0\). In the verification script, the symbol `q` is this Jacobian factor, **not** the particle-production scale in the literature report.

## Metric and shell position

The Jacobian and inverse are

\[
J=\frac12\begin{pmatrix}p+q&q-p\\q-p&p+q\end{pmatrix},
\quad \det J=pq,
\quad J^{-1}=\frac1{2pq}\begin{pmatrix}p+q&p-q\\p-q&p+q\end{pmatrix}.
\]

Exact multiplication gives \(J^{-1}J=JJ^{-1}=I\), and

\[
(J^{-1})^\mathsf T\operatorname{diag}(-1,1)J^{-1}
=\frac1{pq}\operatorname{diag}(-1,1).
\]

Consequently

\[
\widetilde B=B-\tfrac12\log(pq),\qquad
\widetilde A(T,Z)=A(t,z),\qquad
\widetilde\phi(T,Z)=\phi(t,z).
\]

On the shell \(u=v=t\), so \(U=V=F(t)\), \(Z=0\), \(T=F(t)\). If \(L=F'(t)>0\), then \(p=q=L\), and the shell inverse Jacobian is \(L^{-1}I\). The same function on both null directions is essential to this fixed-shell result. No second-brane condition is imposed here.

## Derivatives and invariants at the shell

Let \(w=F''(t)/F'(t)\), and evaluate all quantities below on \(z=0\). For a scalar \(f=A\) or \(\phi\),

\[
f_T=f_t/L,\qquad f_Z=f_z/L,\qquad
f_{TT}=(f_{tt}-wf_t)/L^2.
\]

The logarithmic shift in \(\widetilde B\) is not a scalar transformation. With \(r=F''(u)/F'(u)\), \(s=F''(v)/F'(v)\), its original-coordinate derivatives are

\[
\partial_t\widetilde B=B_t-(r+s)/2,\qquad
\partial_z\widetilde B=B_z+(r-s)/2.
\]

On the shell \(r=s=w\); applying the inverse Jacobian gives

\[
\widetilde B_T=(B_t-w)/L,\qquad \widetilde B_Z=B_z/L,
\qquad e^{-\widetilde B}=L e^{-B}.
\]

It follows immediately that

\[
d\tau=e^{\widetilde B}dT=e^Bdt,
\qquad
H_b=e^{-\widetilde B}\widetilde A_T=e^{-B}A_t,
\qquad
\frac{d\phi_b}{d\tau}=e^{-\widetilde B}\widetilde\phi_T=e^{-B}\phi_t.
\]

The proper derivatives also agree:

\[
\frac{dH_b}{d\tau}
=e^{-2\widetilde B}(\widetilde A_{TT}-\widetilde B_T\widetilde A_T)
=e^{-2B}(A_{tt}-B_tA_t),
\]

\[
\frac{d^2\phi_b}{d\tau^2}
=e^{-2\widetilde B}(\widetilde\phi_{TT}-\widetilde B_T\widetilde\phi_T)
=e^{-2B}(\phi_{tt}-B_t\phi_t).
\]

These are invariants along the same brane worldline, not statements that a field's coordinate derivative is invariant.

## Normal and junction conditions

On the shell,

\[
\widetilde n=e^{-\widetilde B}\partial_Z
=Le^{-B}\left(L^{-1}\partial_z\right)=n.
\]

Thus the scalar normal derivative and mixed extrinsic-curvature eigenvalues agree:

\[
\widetilde n(\phi)=n(\phi),\quad
\widetilde K^i{}_j=e^{-\widetilde B}\widetilde A_Z\delta^i_j
=e^{-B}A_z\delta^i_j,\quad
\widetilde K^T{}_T=e^{-\widetilde B}\widetilde B_Z=e^{-B}B_z.
\]

The trace \(K=e^{-B}(B_z+3A_z)\) is unchanged. Covariant time components transform rather than remaining numerically identical: \(\widetilde h_{TT}=h_{tt}/L^2\) and \(\widetilde K_{TT}=K_{tt}/L^2\).

Since \(\phi_b\) is unchanged, \(\sigma(\phi_b)\) and its field derivative are unchanged. Each normal-form junction residual is therefore exactly invariant, including the registered signs and factors \(1/6\), \(-1/2\). Each coordinate Neumann residual transforms by the nonzero multiplier \(1/L\). Its zero set is unchanged.

## Proper-clock specialization

Locally choose

\[
F'(t)=e^{B_b(t)},\qquad
F(t)=F(t_0)+\int_{t_0}^{t}e^{B_b(s)}ds.
\]

Then \(T=\tau+\mathrm{constant}\), \(\widetilde B_b=0\), and \(w=B_{b,t}\), so \(\widetilde B_{b,T}=0\) and \(H_b=\widetilde A_{b,T}\). The normal derivative \(\widetilde B_Z\) need not vanish. An arbitrary positive \(F'\) does not set the brane lapse to one; the script contains a nonidentity control for that distinction.

A finite proper-time range covered by an old chart remains a finite range under this relabeling. No coordinate transformation manufactures later solution data. A new evolution would need consistent Cauchy data on a new spacelike slice; the image of \(t=\mathrm{constant}\) generally has spatially varying \(T\). A practical gauge implementation additionally needs future gauge information, boundary treatment and numerical stability. None of these is certified here.

## Canonical field dimensions versus registered normalization

The provided `derive_evolution_equations.py` explicitly defines

\[
R_{AB}-\partial_A\phi\partial_B\phi-\frac23U(\phi)g_{AB}=0,
\qquad \Box\phi=U_{,\phi}.
\]

These equations are consistent with the bulk normalization

\[
S_{\rm bulk}=\frac1{2\kappa_5^2}\int d^5x\sqrt{-g}
[R-(\partial\phi)^2-2U(\phi)].
\]

With physical coordinates, \([\kappa_5]=\mathrm{mass}^{-3/2}\), \(\phi\) is dimensionless, \([U]=\mathrm{mass}^2\), and the canonically normalized dimensionful scalar is \(\Phi=\phi/\kappa_5\), with canonical potential \(V(\Phi)=U(\kappa_5\Phi)/\kappa_5^2\). This gives \([\Phi]=\mathrm{mass}^{3/2}\) and \([V]=\mathrm{mass}^5\).

For a canonical four-dimensional brane field \([\chi]=\mathrm{mass}\), the proposed interaction \(g_b^2(\Phi_b-\Phi_*)^2\chi^2/2\) requires \([g_b]=\mathrm{mass}^{-1/2}\). Equivalently, using the dimensionless registered field gives the dimensionful coefficient \(\bar g=g_b/\kappa_5\), with \([\bar g]=\mathrm{mass}\). A prospective crossing scale is \(q_{\rm prod}=\bar g|d\phi_b/d\tau_{\rm physical}|\), of dimension mass squared.

This is conditional dimensional restoration, not a measurement or choice of \(\kappa_5\), a supplied coupling, or the numerical coordinate-to-physical-length calibration. None of those numerical values has been inferred or imposed here. A code-time derivative cannot be used in a production formula without those calibrations.

Local convention source, read-only:
`/Users/ricardomaldonado/Documents/D-Blast 3/untitled folder 150/HDBLAST_CHAT13_NONLINEAR_ROLLOFF_20260921/derive_evolution_equations.py`

SHA-256: `87082ca89e06f63b1cc588df9a054388b1b1dcbb22d11f00d51ef42a4d52a6fc`.

The coordinate freedom has established prior art in [BraneCode, appendix A](https://arxiv.org/abs/hep-th/0309001); the identities above specialize it to the stated single-shell conditions. No novelty claim is made.

## Reproduction and scope

Run from this directory:

```text
python3 -I -S -B -O verify_proper_clock.py
```

The script writes `PROPER_CLOCK_EXACT_CHECKS.json`. It checks 54 exact statements: rational Laurent-polynomial identities, two explicit nonidentity controls, and four exact mass-dimension relations. Rational coefficients cancel exactly; no floating-point tolerance is used. Guards raise exceptions explicitly and are unaffected by optimized Python execution. The proof assumptions, particularly positive Jacobian and smoothness, are premises, not quantities the algebra can certify for an evolved numerical solution.
