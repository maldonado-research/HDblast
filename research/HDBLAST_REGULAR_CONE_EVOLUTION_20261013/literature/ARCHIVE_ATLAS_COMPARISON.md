# Archived negative-branch cone atlas versus the current positive branch

Read 13 October 2026. This is a comparison of private project manuscripts, not external primary-literature validation or a novelty assessment. No archived numerical verifier, release reproduction, PDE estimate, or upstream dependency was rerun. Exact input-body hashes and selective reading scopes are recorded in `PUBLIC_PRIMARY_READING_REGISTER.json`.

## Inputs and claim boundary

Main inputs: PHYS-M471 *Exact conformal first-cone normal form for the sixteen-field harmonic Einstein–scalar system* (complete reading, §§1–11), and PHYS-M465R1 *Cone transparency, finite focusing and scalar transfer* (delegated complete reading, §§1–11). Supporting inputs: M465 master certificate/release audit; M464 main endpoint/ambient theorem, selected endpoint/ambient/FRW passages; and M471 *Principal-scalar two-sided terminal scattering theorem*, complete reading. Current scope checked against the previous cone capsule's `RESULTS.md` and manuscript action, generic master, geometry and conclusion passages.

Archived `PASS_*` headings record claims of those lanes. This reading does not independently reproduce their algebra, authenticated releases, interval atlases, or cited existence theorems. M465's certificate explicitly distinguishes checked principal-symbol algebra from its non-machine-certified PDE energy estimate and comparison theorem.

## Same bulk action, different branch and shell realization

The archive superpotential `W=1−phi+phi³/3` is exactly the current `W=1/3+eta²+eta³/3` under `eta=phi−1`; the bulk potential `U=W_phi²/2−2W²/3` therefore is the same polynomial. The current bulk action is `R/2−(nabla phi)²/2−U`, giving the reusable Ricci-form equations `R_AB=phi_A phi_B+(2/3)U g_AB` and `Box_g phi−U_phi=0` quoted by M465 §2. The currently registered shell is `2W+delta(1+c+c eta)`. Frozen archived shell laws or endpoint certificates do not certify this shell's perturbative boundary conditions.

The endpoint tube and its sign/numerical gates differ:

| Stationary reference | phi=−1 archive | phi=+1 current reference |
|---|---:|---:|
| W | 5/3 | 1/3 |
| U | −50/27 | −2/27 |
| U_phi_phi | 76/9 | 28/9 |
| AdS scale kappa=sqrt(−U/6) | 5/9 | 1/9 |
| Reference homogeneous scalar U_phi_phi/kappa² | 684/25 | 252=14(14+4) |

These exact stationary evaluations are not numerical enclosures of a finite shell solution. The archive's near-−1 `U_phi>0` gate and endpoint/scale intervals must be recomputed. The current cone manuscript instead records `p=phi'<0` on its nonconstant positive branch; its regular endpoint has `p=(U_phi,h/5)y+O(y³)`, so the archived scalar-gradient sign cannot be copied. At the constant reference `p=0`, divided master variables require separate treatment.

## Reusable local geometry and PDE structure

M464's selected endpoint passages use `rho=y a(y²)`, `phi=Phi(y²)` and ambient `s=X·X` to write `g=C(s)eta+D(s)X_flat X_flat`, `C+sD=1`. Analytic coefficients with `a(0)=1` make this an analytic nondegenerate metric across the coordinate cone; tensor equations extend through it if their residuals are analytic and vanish on the exterior open set. This is a background extension, not arbitrary perturbed-data existence. The current manuscript already gives the corresponding exact conformal-radius architecture, Eq. (17).

M471 §§1–2 provides a useful explicit optical reparametrization:

`(log L)_s=(a⁻¹−1)/(2s)`, `Z=L(s)X`, `sigma=s L²`, `Omega=a/L`, hence `g=Omega(sigma)² eta` and `sigma_s=L²/a`.

The apparent singular ODE is regular because `a(0)=1`. This construction transfers structurally after proving the current coefficient regularity and invertibility. It does not transfer the archived common interval tube or its quantitative bounds. With fresh endpoint values `u=U_h,u1=U_phi,h,u2=U_phi_phi,h`, the archived connection/curvature jets are candidate algebra to rederive, not current checked receipts.

M465 §3 and M471 §§5–6 reduce **15 symmetric metric components plus one scalar** in harmonic gauge. M471 (5.1)–(5.2) retains the first-order scalar–metric couplings and algebraic curvature/Hessian terms; the diagonal wave principal part is normally hyperbolic on the smooth ambient metric. Its harmonic constraint and residual-gauge vector operator are (6.1)–(6.2). Compatibility and propagation of subsidiary data are essential: these formulas do not make sixteen arbitrary profiles on one cone admissible characteristic data.

M465's curvature diagnostic `delta R−F'(Phi) delta phi` is not a closed reconstructed physical scalar master. M471's scalar-probe discussion likewise is not the current generic Einstein–scalar master equation. A proposed regular-master transformation must be checked against the current constraints, reconstruction, exceptional harmonics, actual endpoint coefficients, and the `p=0` limit; normally hyperbolic full-system equations alone do not establish equivalence of completed physical master spaces.

## Linear local crossing versus stronger claims

M465 §3 applies a standard finite local wave estimate on an assumed smooth compact ambient patch. It is a linearized harmonic problem on the fixed background. No quantitative common-patch constant is validated by that certificate. M471 §§8–9 distinguishes its local analytic normal form from a complete weighted scattering theorem: the physical indicial spectrum, constraint quotient, uniform high-harmonic control and external/precone source map remain held. Complete postcone H4 slices inside the orbit tube are downstream surfaces; their existence does not specify incoming data that populate the first cone.

Those global holds do not themselves obstruct an ordinary local Cauchy construction on a smooth ambient patch. A new statement must specify its spacelike data, constraint and gauge compatibility, support/domain of dependence, and whether it concerns linearized perturbations or nonlinear solutions. A nonlinear upgrade requires quasilinear harmonic evolution, constraint-satisfying Einstein–scalar data and constraint propagation; the archived linear system does not supply that upgrade. Under perturbation, the fixed background surface `X·X=0` need not remain null. Crossing that coordinate surface and crossing the perturbed metric's own null hypersurface are different assertions.

Neither main archive supplies a constraint-preserving doubled-shell timelike IBVP. M465 §10 explicitly excludes shell/brane perturbations, radion stability, canonical-energy positivity and complete scalar-master reduction; M471 §11 excludes shell/junction perturbations. The current prior capsule likewise leaves registered radiative existence with shell matching and global gauge coverage unestablished. Interior local existence and later shell evolution should have separately stated results and hypotheses.

## Material warning about stationary transfer and terminal scattering

M465 §6's negative-branch stationary comparison has noninteger Gegenbauer degree `18/5` and its selected nonzero second-cone coefficient. Reusing the same reference ODE at phi=+1 gives instead

`chi(theta)=2F1(−14,18;5/2;sin²(theta/2))`.

This truncates to a polynomial regular at both cones, with `chi(0)=chi(pi)=1` and zero limiting `b⁴ chi_dot`. Thus the archived selected-mode obstruction and its transfer signs do **not** migrate to the positive stationary reference. This exact algebraic observation does not predict the finite branch's nonlinear fate or stability. The corresponding reference focusing scale is `9pi`, versus `9pi/5`; finite-branch comparison bounds require its actual `U_h`. M465 §7's floating-point crunch/scaling diagnostics remain unvalidated nonlinear claims.

The other M471 file is a **terminal** scattering lane on the archived M466 kinetic tail, not a first-cone crossing theorem. Its Eq. (1.10) explicitly begins on a regular postcone slice and relies on frozen M466–M470 quotient/tail inputs; §8 continues to hold the ambient/apex-to-slice map and nonlinear/shell claims. Neither its terminal background nor its frozen derivative/Sobolev bounds can be imported into the current positive shell branch without new analysis. Its dependencies were not independently audited here.

The reusable content is the action normalization, regular ambient-coordinate architecture and harmonic-wave/constraint structure. The archives provide no external novelty certification, new cosmological-origin evidence, or replacement for a freshly proved current-branch master/reconstruction/Cauchy theorem.
