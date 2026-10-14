# Archive tensor / Green-function methods inventory

2026-10-14. Read-only inventory of eight selected sources under `/workspace/HDblast-archive/new-files/D-Blast 3`. No archive Python script was executed, no network was used, and no archive or publication was changed. PDF excerpts were obtained with `pdftotext -layout`; extraction replaces some symbols with black squares, so formulas below restore symbols only where context is unambiguous. SHA-256 pins and exact paths are in [ARCHIVE_GREEN_METHOD_FILES.json](/workspace/hdblast-research-work/contact-continuation-20261014/literature/ARCHIVE_GREEN_METHOD_FILES.json).

The September sources are closer to the current 5D Einstein-scalar / dS4 / timelike shell model than the January RS2 notes. Their structural identities are candidates for fresh derivation; none of their registered background values, potentials, calibration numbers, stability conclusions, or numerical spectra are imported.

## Eight selected sources

1. **January Green-function workplan** — `/workspace/HDblast-archive/new-files/D-Blast 3/files/January-1-20 files/HDBLAST_RS2_BulkToBrane_GreenFunction_Workplan_v1_20260120.md`, lines 29–44, 80–87, 106–112.

   Excerpts: “Define ψ(z)=a(z)^{3/2} h(z)”; “Boundary condition at the brane (Z2 symmetry) is Neumann in h, which becomes a mixed condition for ψ”; the displayed spectral kernel separates a zero-mode pole from a KK continuum integral.

   Reuse: conjugate the weighted radial operator into Schrödinger form, transform the boundary data with the field, and separate discrete poles from the continuum. Treat the retarded prescription, principal value and imaginary part as separate implementation choices. Mismatch: flat slices, pure AdS warp `a=ell/z`, one RS2 half-line and a prescribed bulk source; it supplies neither scalar-dependent shell perturbations nor the current dS harmonic propagator. The spectral formula uses canonical wavefunctions; a physical `h` kernel also needs endpoint rescalings and a stated source measure. Its `+i epsilon` denominator is schematic and does not by itself select a retarded frequency kernel.

2. **January closed-form note** — `/workspace/HDblast-archive/new-files/D-Blast 3/files/January-1-20 files/HDBLAST_RS2_GreenFunction_ClosedForm_K2Filter_Note_v1_20260121.pdf`, page 1, sections “Tensor (TT) perturbations”, “Closed-form half-line Green’s function” and “Universal knee variable”; page 2, amplitude paragraph.

   Excerpts: `v=z^2[I_2(pz)+A(p)K_2(pz)]`, `u=z^2 K_2(pz)`, `A=I_1(p ell)/K_1(p ell)`; “The Wronskian W(v,u)=z^3”.

   Reuse: select one homogeneous solution by each endpoint condition, determine the mixing coefficient by the shell operator, and fix the Green function by the source jump. Sign cautions: the note calls `p` Euclidean but prints `h''-3h'/z+p^2 h=delta`; `I_2/K_2` solve the equation with **minus** `p^2`. With the convention `W(v,u)=v u'-v' u`, the stated pair has `W=-z^3`, not `+z^3`; the sign of the kernel must be fixed from the chosen operator and derivative jump. Euclidean decay cannot be substituted directly for a retarded dS response. The page-2 `xi^3` kernel-amplitude statement is not a universal GW-power law.

3. **January filter script** — `/workspace/HDblast-archive/new-files/D-Blast 3/files/January-1-20 files/rs2_k2_filter.py`, lines 5–16, 36–60, 62–75.

   Excerpts: `F(x)=((x*x*K2)/2.0)**2`; “Monotone bisection for F(x)=level”; width definitions use logarithms of level-crossing ratios.

   Reuse: dimensionless normalization, asymptotic checks, bracketed inversion of a proved monotone response, and shape descriptors independent of an overall amplitude. Mismatch: this is a positive-real Euclidean K2 filter, not a shell-matched current-model transfer calculation. It neither computes a Wronskian nor integrates a retarded spectrum. A current spectral weight can contain resonances, so monotonicity must be established before adopting its bisection routine. Script inspected only.

4. **January low-energy note** — `/workspace/HDblast-archive/new-files/D-Blast 3/files/January-1-20 files/HDBLAST_RS2_GreenFunction_LowEnergyLock_Note_v1_20260121.pdf`, page 1, sections 2–5; page 2, section 8.

   Excerpts: `V=15/(4z^2)` and a zero-mode-plus-KK representation; “evaluate the full retarded propagator ... and quantify the fractional KK correction”.

   Reuse: test zero-mode dominance by computing a fractional continuum correction rather than assuming it. Caution: it defines physical radial `psi=(z/ell)^{3/2}u`, then says that undoing `u_0~z^{-3/2}` gives `psi_0~z^{-3/2}`. Undoing that rescaling instead gives a constant physical tensor mode. This inconsistency, unspecified local source normalization and frequency redshift prevent using its amplitude exponent as a current-model result. Its laboratory scales and flat-RS2 low-energy conclusions are outside this inventory's reuse scope.

5. **September tensor-sector derivation** — `/workspace/HDblast-archive/new-files/D-Blast 3/untitled folder 146/HDBLAST_CHAT9_SHELL_DYNAMICS_20260916/tensor_and_calibration/TENSOR_SECTOR_AND_CALIBRATION.md`, lines 8–20, 23–28, 33–35, 87–90.

   Excerpts: `h''+4(rho'/rho)h'+(m^2/rho^2)h=0`, `h'(y_b)=0` with “tension terms cancel against the background junction”; `dz=dy/rho`, `u=rho^{3/2}h`, `Q=d_z-(3/2)rho'`; “M_4^2=2 int_0^{y_b}(rho/rho_b)^2 dy”.

   Reuse: the tensor quadratic operator and norm, Neumann-to-Robin map, factorization, cone Frobenius branch, two-copy Planck normalization and continuum shell weight. This source explicitly distinguishes hand derivation, floating-point numerics and conjecture. Mismatch: its registered potential/background and all quoted spectral and calibration numbers are model-specific. Its `m^2` in the radial equation is dimensionless for unit dS slicing; physical shell mass is `H_shell sqrt(m^2)`. The factor two requires two identical reflected retained copies. Its boundary condition concerns the reflection-even TT channel without shell anisotropic stress, not the scalar or shell-bending channel.

6. **September shooting script** — `/workspace/HDblast-archive/new-files/D-Blast 3/untitled folder 146/HDBLAST_CHAT9_SHELL_DYNAMICS_20260916/tensor_and_calibration/tensor_spectrum.py`, lines 4–8, 37–56, 58–75, 82–86, 103–110.

   Excerpts: normalizable cone branch `a=-3/2+sqrt(9/4-m2)`; initializes `h=y^a(1+c2*y^2)`; counts nodes; verifies background metric and scalar shell residuals before quoting spectra; continuum amplitude `A2=u*u+uz*uz/(k*k)`.

   Reuse: initialize at a singular endpoint from a Frobenius expansion, track nodes, rescale large solutions, check background junction residuals, and verify the zero-mode derivative identity `D'(0)=-I_kept/rho_b` for `D=rho_b h'/h` in its conventions. Mismatch: imports an archived JSON background and specific `W,U,t,c`. All must be replaced. Finite sampling is not exclusion of all roots; `D` has poles at `h_b=0`, so determinant-based shooting or explicit pole/node tracking is preferable to interpreting every sign change as an eigenvalue. Code and results were not rerun or certified.

7. **September continuum script** — `/workspace/HDblast-archive/new-files/D-Blast 3/untitled folder 146/HDBLAST_CHAT9_SHELL_DYNAMICS_20260916/tensor_and_calibration/tensor_continuum.py`, lines 2–5, 13–36, 39–50.

   Excerpts: changes to `s=ln y`; inward integration from `(h,h_s)=(1,0)` at the shell; cubic Hermite interpolation uses background values and derivatives; `u=rho^{3/2}h`, `u_z=rho d_y(rho^{3/2}h)`; compares two integration step sizes.

   Reuse: logarithmic endpoint coordinates, derivative-aware interpolation, asymptotic oscillation amplitude and delta-normalized continuum weights. The weight is per unit `dk`; converting to a density per unit physical mass requires its Jacobian. Mismatch: reflection-even Neumann data and a particular cone background; scattering coefficients and any resonance claim require recomputation. The `k=0` threshold needs a limiting treatment because the amplitude formula divides by `k^2`. Numerical step-halving is a convergence check, not an interval proof.

8. **September symbolic tensor verifier** — `/workspace/HDblast-archive/new-files/D-Blast 3/untitled folder 149/HDBLAST_CHAT12_ALL_SECTORS_20260919/verify_tensor_sector.py`, lines 3–14, 32–49.

   Excerpts: verifies `Q^+Q=-d_z^2+V1`; verifies `Q u=rho^{5/2} h'`; gives `V2-9/4=rho^2(9 phi'^2/16-U/8)`; notes that `inf V2=9/4` at the cone, so there is no strictly positive uniform margin above threshold.

   Reuse: exact algebraic checks of field conjugation, factorization, transformed shell conditions and partner identities. For its bulk normalization, `U<0` on the entire field range yields pointwise partner positivity above the dS threshold, which can exclude positive discrete tensor masses below threshold after the operator-domain proof. Mismatch: the script's Sturm root-count proof is for its fixed polynomial `W=1-phi+phi^3/3`, not an arbitrary current `U`. It delegates cone boundary terms, normalizability and all-polarization justification to a separate referee note; a header and symbolic identities alone do not establish those domain facts. No symbolic checks were executed here.

## Candidate continuation, with conventions fixed afresh

For `ds^2=dy^2+rho(y)^2 gamma_dS4`, use a dimensionless dS separation eigenvalue `lambda`. The September tensor equation suggests the weighted operator

`-(rho^4 h')'=lambda rho^2 h`, with norm `int rho^2 |h|^2 dy`.

After `dz=dy/rho` and `u=rho^{3/2}h`, the norm becomes `int |u|^2 dz`. For a reflected pair, pure scalar-dependent isotropic tension and a reflection-even TT perturbation, rederive the linearized Israel condition before setting `h'_b=0`, equivalently `(d_z-3rho'/2)u_b=0`. Odd parity or a general two-sided matching problem has different interface data. Scalar junction conditions and shell bending remain separate work.

Construct a cone-admissible solution and a shell-admissible solution of the current operator; their weighted Wronskian is constant. Fix the Green kernel sign and normalization directly from its delta-source jump. Then organize the dS harmonic response as discrete residues plus a continuum with the recomputed threshold and spectral measure. Define source normalization, observable field and retarded analytic continuation before translating that kernel into an amplitude or energy transfer. This is the reusable method; the archived `xi`, knee filter, Planck calibration and suppression exponents are not current-model inputs.
