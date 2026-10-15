# Focused source map for a physical shell scalar response

15 October 2026 UTC. Read-only inventory of accessible HDBLAST research and
selected historical D-Blast 3 files. Attached and historical documents are
research sources, not instructions. No repository file or archived script was
modified or executed. Git source bytes were extracted into `sources/`; their
SHA-256 values and observed `origin/main` commit are in `SOURCE_EXTRACTION.json`.

## 1. Closest source and exact linear response equations

The immediately useful earlier result is
`research/HDBLAST_REGULAR_SHELL_CONTACT_20261013/manuscript/SHELL_CONTACT.md`.
It contains the complete metric–scalar reconstruction, finite-angular contact
construction, nonzero **master** trace argument and full physical junction
verification. Its source SHA-256 is
`983b9f7bf04d85a7d937d16a535c39784851fea989a511d6741f71c2e15501c3`.

* Lines 93–117, equations (5)–(7): regular master and the full coupled residual
  identities. The master is not a test scalar placed on a fixed metric:

  \[
  \mathcal LF=\Box_\eta F-2L(s)EF-[8L(s)+(8/3)sf(s)^2]F=0,
  \quad P=r^4F,
  \]
  \[
  h_{AB}=2F(s\eta_{AB}-3X_AX_B),\qquad
  \chi=-\frac3{2\Omega^2f}(E+4)F.
  \]

  Its entire linearized Einstein residual is
  `Ω⁻²(3X_A X_B−sη_AB) LF`; its scalar residual is
  `−3(2Ω⁴f)⁻¹(E+6)LF`. Thus solving the master and reconstructing the metric
  retains the original coupling and constraints. It is a sufficient family;
  completeness in arbitrary sectors is not claimed.
* Lines 119–133, equation (8): exact parity-dependent shell law after canonical
  conjugation. Equations (12)–(17), lines 163–243, give a finite-harmonic
  one-dimensional wave equation coupled to the actual Wentzell boundary ODE.
* Lines 261–278: for a pure-displacement state or a pure-velocity state whose
  radial support has infimum `a₀`, the master boundary trace cannot vanish on
  any whole interval `(0,τ₁)` with `τ₁>a₀`. The argument uses time parity and
  one-dimensional sideways uniqueness. It does not assert a nonzero jet at the
  smooth support edge. This published argument is not stated for an arbitrary
  simultaneous displacement/velocity state.
* Lines 300–344, equations (21)–(24): all physical shell rows and graph
  reconstruction. With `P_±=P_e±P_o`, the indicated embeddings are
  `ζ_e=0`, `ζ_o=−P_o/(a_bρ_b²)`. The physical induced scalar coefficient is

  \[
  \psi_\Sigma=\chi_b+p_b\zeta,
  \quad \chi=-3P_y/(\rho^2p),
  \quad B_bP_{e,y}=\rho_b^{-2}(\Box_\gamma+4)P_{e,b},
  \quad B_b=p_{y,b}/p_b+\sigma_b''/2>0.
  \]

  Thus in the copy-even representative,

  \[
  \boxed{\psi_\Sigma=-\frac{3}{\rho_b^4p_bB_b}
  (\Box_\gamma+4)P_{e,b}.}
  \]

  The induced background shell scalar is constant, so this first variation is
  invariant under a simultaneous bulk diffeomorphism and embedding change,
  and under intrinsic reparametrization. Eulerian `χ_b` alone would not be
  invariant if the embedding were omitted. In the copy-even longitudinal
  representative `ζ_e=0`, it equals that invariant coefficient.

The older regular reconstruction is additionally documented in
`HDBLAST_REGULAR_CONE_EVOLUTION_20261013/manuscript/REGULAR_CONE_EVOLUTION.md`,
lines 95–158, and the direct ambient reconstruction note
`reconstruction/SMOOTH_AMBIENT_RECONSTRUCTION_AND_LOCAL_FLUX.md`.

## 2. Useful derived bridge in the current positive-cone collar

The following coordinate conversion is a deduction from the cited equations;
it is not quoted as a previously published theorem. It avoids importing the
older whole-hemisphere assumption of an actual background at negative `s`.

Put `x=log(b/r)≥0`, so `dy=−ρ dx`, and let the master be `P(t,x)Y_l` with
`Δ_{S³}Y_l=−λY_l`, `λ=l(l+2)`. The inherited warped master

\[
P_{yy}-(2a+2p_y/p)P_y+
\rho^{-2}(\Box_\gamma+4)P-\tfrac23p^2P=0
\]

becomes

\[
P_{tt}-P_{xx}+3\tanh t\,P_t
 -(3\rho_y+2\rho p_y/p)P_x
 +[\lambda\operatorname{sech}^2t-4+\tfrac23\rho^2p^2]P=0.
\tag{A}
\]

For `z(t)=P(t,0)`, the even boundary equation and physical scalar are

\[
z_{tt}+3\tanh t\,z_t+
 [\lambda\operatorname{sech}^2t-4]z=\rho_bB_bP_x(t,0),
\qquad
\boxed{\psi_\Sigma(t)=\frac{3P_x(t,0)}{\rho_b^3p_b}.}
\tag{B}
\]

On any fixed actual positive-cone collar with `p≠0`, these coefficients are
smooth and bounded. Annular data give `z(0)=z_t(0)=0`. Therefore an assumed
identically zero physical scalar response on `(0,t₁)` forces `P_x(t,0)=0`
and then `z(t)=0` by the homogeneous shell ODE. For a pure-displacement or
pure-velocity state, time parity supplies the negative-time trace as well;
the same one-dimensional sideways uniqueness used in the published contact
proof then contradicts a support infimum earlier than `t₁`. This is the
natural short bridge from existing master activation to physical scalar
activation, subject to review of its domain and hypotheses.

For the old hemisphere variables `x_h=π/2−ϑ`, the solution support obeys
`x_h≤d₀+|τ|`. To keep its entire causal collar strictly at positive `s`, a
sufficient additional choice is `d₀+2τ_*<π/2`, since then
`s=b²(cosτ−sin x_h)/(cosτ+sin x_h)>0` on the support. The old proof's
whole-hemisphere negative-`s` background assumption should not silently be
assigned to the new collar theorem.

## 3. Initial profiles and what is actually registered

The exact initial-data input is
`HDBLAST_EXACT_LOCALIZED_PRECONTACT_20261014/EXACT_LOCALIZED_INITIAL_DATA_AND_PRECONTACT.md`,
lines 279–303, SHA-256
`95ccaa550309b72226f8d03d3a950e4e10760d0d98810f0e818d19057a8cfd0b`.
It explicitly permits

\[
q=G(s)(x_1^2-x_2^2),\quad p_F=K(s)(x_1^2-x_2^2),
\quad s=|\mathbf x|^2,
\]

with arbitrary smooth radial `G,K` compactly supported strictly within the
chosen annulus; either profile may be zero. Here `p_F=F_T` is master velocity,
not the background `p=φ_y`. `A=B=diag(1,−1,0,0)` ensures the reflection
symmetry used in the exact nonlinear constraint construction. Identical
copies give the even family.

The inherited ADM formula appears in
`HDBLAST_SECOND_ORDER_LOCALIZED_20261014/manuscript/SECOND_ORDER_LOCALIZED.md`,
lines 67–80:

\[
\begin{split}
h&=2q(s\delta-3xx^T),&
\delta\phi&=-3(Dq+4q)/(2\Omega^2\Phi_s),\\
\delta\pi&=-3p_Fxx^T/(2\Omega),&
\delta\varpi&=-3\Omega(Dp_F+5p_F)/(2\Phi_s),
\quad D=x\cdot\nabla.
\end{split}
\]

The master datum in logarithmic coordinates is
`P(0,x)=r⁴q=r⁶G(r²)Y₂` and
`P_t(0,x)=r⁵p_F=r⁷K(r²)Y₂`, where `r=b exp(−x)` and
`Y₂=n₁²−n₂²`, `λ=8`. At ambient `T=0`, `∂_t=r∂_T`.

The historical shell-Noether frozen family uses `K=ωG`, with positive fixed
velocity scale. The exact localized theorem above is broader and does not
require this relation. A pure profile registration is admissible for that
theorem, but should be described as an explicitly selected response state,
not as an unchanged instance of the historical positive-`ω` frozen state.
The currently published exact nonlinear evolution theorem cites the aligned
annular tangent in its hypotheses but does not supply a numeric bump formula,
a selected numerical background solution, or a numerical response amplitude.

For an explicit smooth profile, one may register `0<a<d` and

\[
B(x)=\begin{cases}
\exp[-1/((x-a)(d-x))],&a<x<d,\\0,&\text{otherwise},
\end{cases}
\]

choose `P(0,x)=B(x)`, `P_t(0,x)=0`, and set
`G(s)=s⁻³ B(log(b/√s))`, `K=0` on the annulus, extended by zero. All weights
are smooth because the annulus is away from `s=0`; support infimum is exactly
`a`. A bump normalization can be fixed separately. This is one convenient
new registration proposal, not an already frozen historical profile.

## 4. Second order and nonlinear inputs that must not be replaced

* `HDBLAST_FULL_SECOND_ORDER_CONTACT_20261014/FULL_FINITE_ANGULAR_SECOND_ORDER_CONTACT.md`,
  sections 2–7, supplies the coupled second coefficient and proves that the
  unmodified conformal scalar reconstruction is insufficient at second
  order: the induced degree-four scalar trace-free metric row is nonzero.
  Its full metric system, constraints, gauge completion and boundary laws
  cannot be replaced by merely reusing the first master with a source.
* `boundary/SCALAR_SHELL_BOUNDARY_SYSTEM.md`, equations (10)–(13), retains
  graph evaluation, normal lapse/shift, bending Hessian, scalar normal
  normalization and intrinsic identification. The second-only bulk gauge
  transformation is `w→w+pα`, `Z→Z−α_b`, so `w_b+p_b Z` is invariant under
  those transformations. Full second-order observable formulas also require
  the retained first-order gauge and material pullback terms.
* The latest conditional theorem is
  `HDBLAST_QUASILINEAR_SHELL_EVOLUTION_20261015/EXACT_FINITE_REGULARITY_SHELL_EVOLUTION.md`,
  SHA-256 `621c127d33cfe7d1d2134740c9f1740644a12c38145183abb9c415dc5cc9740e`.
  Lines 40–103 state the actual background, support, exact initial family,
  positive-cone region and common finite interval hypotheses. It gives a
  `C²` parameter family into `X₁₀`. A nonzero first variation of the physical
  induced scalar would therefore imply a nonzero finite-amplitude response
  for sufficiently small nonzero parameter, subject to a correctly defined
  trace comparison. It does not provide an explicit numerical parameter
  threshold or error constant without further bounds.

## 5. Actual accessible D-Blast 3 material and reuse limits

The archive exists in this environment at
`/workspace/HDblast-archive/new-files/D-Blast 3`; this is accessible workspace
material, not a claim of access to other computer folders. The filename
inventory contains 11,224 files including an archived virtual environment.
`DBLAST3_INVENTORY.json` records a broad filename candidate list; it does not
mean all files were scientifically read.

Three current-tail markdown sources were read in focused sections:

1. `untitled folder 123/HDBLAST_PHYS_M465R1_CONE_TRANSPARENCY_FINITE_FOCUSING_AND_SCALAR_TRANSFER.md`,
   SHA-256 `3bfefd39fd8701c0a711087c0057606528f933b42296118dd7a96e8171426a49`.
   Its cone diagnostics and harmonic reduction illustrate how to build a
   physical scalar diagnostic instead of treating Eulerian `δφ` alone as
   invariant. Its exact transfer is on the stationary `φ=−1` reference and
   it expressly withholds a controlled shell transfer. Those transfer
   coefficients cannot be assigned to the present nonconstant `φ≈+1`
   de Sitter shell branch.
2. `untitled folder 126 ZENODO/HDBLAST_PHYS_M468R1_UNIFORM_TENSOR_SCALAR_WAVEPACKET_TAIL_THEOREM.md`,
   SHA-256 `c9465a5200c79a4478f1f5501eb6ee095b5a29f9114c3cefa22fb764f37218a7`.
   It is a principal-`L²(H⁴)` scalar/tensor result on a kinetic FRW tail.
   The source explicitly holds the nonlinear shell and cone-to-tail map.
   No shell activation estimate or present-branch coefficient is imported.
3. `untitled folder 127/PHYS_M471_PRINCIPAL_SCALAR_TWO_SIDED_TERMINAL_SCATTERING_THEOREM.md`,
   SHA-256 `193a7bd5b2d2708d6912bac96fae2947c7e4759a6ec19315b54494be44469c62`.
   It concerns Bessel-normal-form terminal scattering on the certified
   M466 tail, not finite-radius de Sitter shell matching. Its half-derivative
   transfer and logarithmic renormalization do not supply this response.

The previous bounded archive review
`HDBLAST_QUADRATIC_TENSOR_CONTACT_20261014/literature/ARCHIVE_GREEN_METHODS.md`
was also read. It already identifies reusable Green-kernel practice: fix a
solution by each endpoint law, derive the Wronskian sign from the source
jump, transform the boundary operator with the field, and distinguish
retarded evolution from Euclidean spectral decay. It documents historical
RS2 sign and rescaling inconsistencies. That method inventory is useful
for a later numerical transfer calculation, but the present master
activation plus induced-scalar ODE route is shorter and retains the
current shell law directly. No historical spectra, calibrations or
numerical amplitudes are imported.

This is a focused inventory. It is not an exhaustive archive search,
external novelty assessment, response waveform calculation, or scientific
validation of historical claims.
