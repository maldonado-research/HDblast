# Prospective scalar stability bridge from regular-cone Jacobi fields

11 October 2026. **Private, bounded future-work calculation.** This memo is not
part of the frozen October 11 public capsule, a new stability certificate, or a
revision of its claims. It records a promising exact bridge and the remaining
operator/domain work. No full spectral scan or interval stability calculation
was performed for this memo. Source hashes are in `SOURCE_MANIFEST.json`.

## 1. What the bounded investigation found

The newly established static signs supply more than an asymptotic estimate of
the archived boundary coefficient. On every nonconstant interior in the
constructed tube they give

\[
 B_{\partial}:=p'/p+\sigma''/2
   >J:=p_\beta+\sigma''/2\geq w/2>0.
\]

More significantly, an explicit combination of the regular cone Jacobi fields
is the cone-regular solution of the archived gauge-invariant scalar system at
\(\mu^2=0\). It is nodeless, and its scalar boundary mismatch obeys the uniform
normalization-independent bound

\[
 \boxed{M_0/X_0>J\geq w/2.}
\]

For the registered model \(w=2\), this lower bound is at least one. Thus the
anticipated missing estimate of the small scalar derivative or zero-mode
mismatch appears avoidable by an exact identity. What remains is a careful
independent audit of the physical perturbation operator, its exceptional
harmonics, and its singular endpoint and spectral domain. These conclusions
are conditional on the archived operator, not an inference that static
uniqueness by itself implies dynamical stability.

## 2. Exact model and notation dictionary

The September 27 source has
\[
 W_{\rm old}(\phi)=1-\phi+\phi^3/3,
 \quad \sigma=2W_{\rm old}+\delta(1+c\phi).
\]
The exact field shift \(\phi=1+\eta\) gives
\[
 W(\eta)=1/3+\eta^2+\eta^3/3,
 \quad f(\eta)=1+c+c\eta,
 \quad U=\tfrac12(W')^2-\tfrac23W^2.
\]
Hence \(k=1/9,w=2,\lambda_{\rm bulk}=w-4k=14/9\), as in the October 11
static theorem. The coupling in the exact-decimal model is
\(c=2/(1.0357712571566784)-4/3\), with the decimal interpreted rationally.
The September Python calculation evaluates that expression in binary64.
The algebraic identities below hold for either coupling convention, but a
numerical replay is not silently identified with the rational certificate.

Use \(p=\eta'=\phi'\), \(a=\rho'/\rho\), and \(\mathcal H=\rho^{-2}\).
The physical shell curvature is \(H_b^2=\rho(T)^{-2}\).
The archived harmonic convention is
\[
 \Box_\gamma Y=\mu^2Y,\qquad m_4^2=\mu^2H_b^2,
 \qquad\Lambda=\mu^2+4.
\]
The spectral \(\Lambda\) must not be confused with \(\lambda_{\rm bulk}\).
The old numerical `H` in `gi_core.py` denotes \(a\), not \(H_b\).

## 3. Static ingredients already available

Let \(u=\eta(0)\) be the cone datum, and define derivatives at fixed radial
coordinate:
\[
 h=\eta_u,\quad q=p_u=h',\quad v=(\log\rho)_u,\quad b=a_u=v',
 \quad j=b+ph/3,\quad\mathcal D=hp'-pq.
\]
The October 11 proof establishes strict convexity of \(U\) on the actual tube
from its existing contraction condition, together with
\[
 h\geq1,\quad q>0,\quad pj>0,\quad p\mathcal D>0
\]
for each nonconstant interior. It proves
\[
 j'+4aj=-2\mathcal H v,
 \qquad \mathcal D'+4a\mathcal D=4p(\mathcal H h+pj),
\]
and the linearized constraint gives
\[
 pq-p'h=12aj+12\mathcal H v.                       \tag{S1}
\]
Since \(p_\beta=q(T)/h(T)\),
\[
 B_{\partial}-J=\frac{p'}p-\frac qh
                =\frac{\mathcal D}{ph}>0.          \tag{S2}
\]
This proves \(B_{\partial}>J\geq w/2\) without dividing an error estimate by
the small quantity \(p\). The division is only an exact identity on a segment
where \(p\ne0\). For the registered nonzero \(f'(0)=c\), the scalar junction
does not admit \(\beta=0\) for \(\delta>0\); the strict-convexity argument
then gives \(p<0\) everywhere on \((0,T]\). Thus the required nonzero-scalar
derivative hypothesis is actually available for every positive detuning in
the static interval. The constant-scalar case \(f'(0)=0\) requires a separate
perturbation formulation and is not covered by this divided operator.

## 4. A division-free construction of the zero-mass reference mode

The archived gauge-invariant equations, with \(\Psi=pZ\), are
\[
 \Psi'=-2a\Psi-pX/3,
 \qquad
 pX'=p'X+(3\Lambda\mathcal H-2p^2)\Psi.           \tag{G1}
\]
The scalar junction mismatch is
\[
 M=X'+\sigma''X/2+2p\Psi
   =B_{\partial}X+3\Lambda\mathcal H\Psi/p.         \tag{G2}
\]
Set
\[
 K=\rho^2j,\qquad
 X_0=h+pK,\qquad \Psi_0=v+aK.                    \tag{Z1}
\]
All three expressions are division free in \(p\). The Jacobi identity gives
\(K'=-2aK-2v=-2\Psi_0\). Therefore
\[
 \Psi_0'+2a\Psi_0=b+a'K
  =b-j-p^2K/3=-pX_0/3.
\]
Also
\[
 X_0'=q+p'K-2p\Psi_0.
\]
Using (S1) gives
\[
 pX_0'-p'X_0=(12\mathcal H-2p^2)\Psi_0.
\]
Thus (Z1) solves **both** equations (G1) at \(\Lambda=4\), or \(\mu^2=0\).
It is not merely a resemblance between a static parameter variation and a
perturbation equation. Equivalently, the constraint yields the useful exact
formula
\[
 \Psi_0=-\rho^2\mathcal D/12.                     \tag{Z2}
\]

The regular cone expansions fix the mode normalization. Writing
\(u_1=U'(u)\ne0\),
\[
 p=u_1y/5+O(y^3),\quad j=u_1y/90+O(y^3),
 \quad v=-u_1y^2/36+O(y^4),\quad a=y^{-1}+O(y),
\]
so
\[
 X_0=1+O(y^2),\qquad
 \Psi_0=-u_1y^2/60+O(y^4),\qquad
 Z_0=\Psi_0/p=-y/12+O(y^3).                       \tag{Z3}
\]
These agree with the September code's regular data at \(\mu^2=0\):
\(X\sim y^s,Z\sim-y^{s+1}/[3(s+4)]\), with \(s=0\).

Because \(pK=\rho^2pj>0\), one has \(X_0>h\geq1\).
Equation (Z2) gives \(p\Psi_0<0\), so \(\Psi_0\) has no zero on \((0,T]\).
On the registered branch \(p<0\), its sign is positive, exactly the
normalization needed by the archived oscillation criterion. On a positive-
\(p\) branch one may multiply the whole reference mode by minus one.

At the shell, the full mismatch cancels to
\[
 \begin{aligned}
 M_0&=q+\sigma''h/2+K(p'+\sigma''p/2)\\
    &=hJ+B_{\partial}pK>0.                         \tag{Z4}
 \end{aligned}
\]
For \(\alpha_K=pK/h>0\),
\[
 \frac{M_0}{X_0}=
 \frac{J+B_{\partial}\alpha_K}{1+\alpha_K}>J\geq w/2. \tag{Z5}
\]
This avoids a uniform lower bound on \(|p|\) or the exponentially small
central displacement as \(\delta\downarrow0\).

The connection to the static residual is exact but must be interpreted with
care. With \(S=p+\sigma'/2\), \(R=a-\sigma/6\), the Jacobian in variables
\((u,T)\), with rows \((R,S)\), is
\[
 \begin{pmatrix}j&-H_b^2\\ hJ&pB_{\partial}\end{pmatrix},
 \qquad \det=H_b^2M_0>0.                           \tag{Z6}
\]
The static determinant alone would not establish nodelessness. Here that
additional ingredient follows independently from the pointwise Jacobi signs
and the explicit mode construction (Z1)–(Z2).

## 5. Sturm–Liouville conversion and the exact endpoint sign

Use symbols distinct from the background momentum:
\[
 P=\rho^2\Psi,\quad
 p_{SL}=\frac1{\rho^2p^2},\quad
 q_{SL}=\frac2{3\rho^2},\quad
 w_{SL}=\frac1{\rho^4p^2},\quad
 \beta_{SL}=\rho(T)^4p(T)^2 B_{\partial}>0.
\]
The archived system gives
\[
 -(p_{SL}P')'+q_{SL}P=\Lambda w_{SL}P,
 \qquad p_{SL}P'(T)=\frac{\Lambda}{\beta_{SL}}P(T).  \tag{L1}
\]
Indeed \(P'=-\rho^2pX/3\), so the coefficients and the spectral boundary
parameter have the signs shown. For the nodeless reference \(P_0=\rho^2\Psi_0\),
\[
 r_0:=\frac{p_{SL}P_0'}{P_0}(T)-\frac4{\beta_{SL}}
   =-\frac{M_0}{3\rho(T)^2p(T)B_{\partial}\Psi_0(T)}>0. \tag{L2}
\]
The inequality uses \(p\Psi_0<0\), not an assumed sign of the raw mismatch
copied from the old wall-side branch. Equivalently,
\(p_{SL}P_0'(T)-4P_0(T)/\beta_{SL}=-M_0/(3pB_{\partial})\).
This last expression is positive on the registered \(p<0\) branch.

Conditional on the archived scalar Sturm–Liouville realization and its
normalizability assumptions, (Z2) and (L2) meet its zero-mass comparison
criterion. This conditional implication is much narrower than asserting
stability of the complete physical system.

## 6. Ground-state identity: a route around a numerical zero scan

The positive boundary weight gives
\[
 \mathcal N[u]=\int_0^T w_{SL}|u|^2dy
                          +|u(T)|^2/\beta_{SL}>0.
\]
For the differential expression (L1), define
\(\mathcal E[u]=\int_0^T(p_{SL}|u'|^2+q_{SL}|u|^2)dy\).
Integration by parts at a spectral root gives
\(\mathcal E[u]=\Lambda\mathcal N[u]\), provided the cone term vanishes.
With a real nodeless \(P_0\) satisfying (L1)'s bulk equation at \(\Lambda=4\),
the ground-state/Picone identity is
\[
 \mathcal E[u]-4\mathcal N[u]
 =\int_0^T p_{SL}P_0^2\left|(u/P_0)'\right|^2dy
                         +r_0|u(T)|^2,             \tag{L3}
\]
again provided the lower-endpoint term vanishes. Its right side is strictly
positive for any nonzero eigenfunction: equality would force \(u/P_0\)
constant and \(u(T)=0\), hence \(u=0\). Thus a validated realization of
this identity excludes eigenvalues \(\Lambda\leq4\), or \(\mu^2\leq0\).
The associated reality statement comes from the positive norm in the energy
identity. This is a candidate operator-level proof route, not a replacement
for checking which perturbations that operator represents.

For prospective nonpositive eigenmodes the cone term can already be checked
directly. Put \(A_c=U'(u)/5\ne0\). Since
\(p_{SL}\sim A_c^{-2}y^{-4}\), \(w_{SL}\sim A_c^{-2}y^{-6}\), the
indicial equation is
\[
 \gamma^2-5\gamma+\Lambda=0,
 \quad\gamma_\pm=5/2\pm\sqrt{25/4-\Lambda}.
\]
At \(\Lambda=4\), the powers are \(4\) and \(1\); only \(y^4\) is
square integrable with weight \(w_{SL}\). At \(\mu^2\leq0\), the regular
mode is \(u\sim y^{4+s}\), with
\(s=-3/2+\sqrt{9/4-\mu^2}\geq0\). Both relevant cone terms scale as
\[
 p_{SL}\bar u u'=O(y^{2s+3}),\qquad
 p_{SL}(P_0'/P_0)|u|^2=O(y^{2s+3}),
\]
and vanish. The reference mode itself satisfies
\[
 P_0=-U'(u)y^4/60+O(y^6),\quad
 w_{SL}P_0^2\sim y^2/144,\quad
 p_{SL}(P_0')^2\sim y^2/9.
\]
These cancellations show why an explicit lower bound on the tiny central
\(|U'(u)|\) is unnecessary for this normalization. They do not by themselves
establish the closure and spectral completeness of the entire form domain.

## 7. The actual remaining work

1. **Verify the perturbation operator and conventions independently.** The
   September 27 derivation reports 55 symbolic checks and wrong-sign controls.
   This memo reads and matches its equations but does not rerun the complete
   tensor derivation or independently derive the action's second variation.
   The inward/outward orientation, doubled-interior normalization, and full
   displaced-shell condition must remain fixed. The term \(2p\Psi\) in (G2)
   is essential to (Z4).

2. **Separate exceptional harmonics from generic scalar modes.** The archived
   trace-free reduction divides by harmonic structure that degenerates for
   constant harmonics and for the translation-type harmonic. Its homogeneous
   realization has time exponents \(r=0,1\) exceptional, with
   \(\mu^2=-(r^2+3r)\). The generic \(\mu^2=0\) solution can be interpreted
   through the other exponent \(r=-3\), but this does not dispose of the
   constant \(r=0\) sector by itself. Equation (Z6) supplies a strictly
   positive static junction determinant for that separate calculation.
   At \(\mu^2=-4\), the archived relative-displacement exclusion condition is
   \(B_{\partial}p\ne0\), now available from (S2), but its physical/gauge
   classification must still be checked from the unreduced equations.

3. **Audit the operator domain and positive boundary mass.** At each fixed
   \(\delta>0\), the finite shell endpoint is regular and the cone is
   limit-point in the archived realization. The correct Hilbert norm is
   \(L^2(w_{SL}dy)\oplus\mathbb C\) with endpoint weight
   \(1/\beta_{SL}\). Check the trace identification, closed quadratic form,
   self-adjoint operator with eigenparameter boundary condition, and closure
   of (L3) on that form domain. The direct cone calculation above suffices
   for candidate nonpositive eigenfunctions, not automatically every vector
   in the closed form domain. The archived oscillation note explicitly marks
   self-adjointness/completeness as an imported standard theorem.

4. **Define the physical meaning of the resulting spectral exclusion.** A
   conditional exclusion of nonpositive normalizable scalar eigenmodes is
   not a complete semigroup estimate, proof of a uniform mass gap, vector-
   sector analysis, nonlinear attraction theorem, tunnelling analysis, or
   reheating/cosmological scenario. The limiting \(\delta=0\) divided
   operator remains singular; the argument is for each positive detuning.

The next efficient step is a short independent operator/domain review of
(G1)–(L3), including the exceptional harmonics, followed by a clean separate
conditional scalar-mode proposition if that review succeeds. A large new
floating-point spectral scan is not the missing estimate identified here.
One cheap archived-background comparison can detect normalization mistakes,
but cannot settle these domain questions.

## 8. Local verification and limits of this delivery

The independent subdirectory contains a new symbolic checker for the Jacobi-
to-zero-mode transformation, including both first-order equations, center
normalization, the full scalar mismatch, and the Sturm–Liouville boundary
sign. All 22 exact checks passed, and a separate read-only review found no
substantive error in this memo's Picone identity, cone exponents or scope.

A small independently coded floating-point diagnostic at the archived
\(\delta=0.001\) background is recorded in
`numeric/ZERO_MODE_IDENTITY_SANITY.json`. With starting radii \(10^{-4}\) and
\(10^{-5}\), the Jacobi and directly integrated gauge-invariant constructions
both give \(M_0/X_0\simeq3.5547077330\), with relative endpoint differences
below \(2.5\times10^{-16}\). At 301 sampled radial points, the reported
relative profile discrepancies are at most \(1.3\times10^{-12}\) for \(X\)
and \(1.3\times10^{-11}\) for \(\Psi\) reconstructed through \(\mathcal D\).
These are identity diagnostics, not bounds on cone truncation, background
input error, or the continuous perturbation solution. They do not supply
interval certification or settle operator-domain questions.

The newly derived assertions in this memo
have not received external mathematical review. No frozen release file,
metadata, figure, or public theorem was changed.
