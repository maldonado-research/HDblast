# A matter-production extension with a signed energy budget

22 September 2026. This is a proposed extension of the registered HDBLAST Einstein–scalar shell, with exact classical identities and approximation tests. It is not a reheating simulation or evidence of an observed higher-dimensional origin. The pure-tension model and its completed numerical calculations are not changed by this document.

The useful advance is a fully specified way to add matter without losing the gravitational and scalar boundary conditions. A separate exact check exposes a small but real defect in the often-used late-time comparison: an exactly constant bulk scalar at φ=1 cannot satisfy the registered pure-tension scalar junction when δc≠0. A nearby static branch has an analytic starting estimate, derived below, that can be tested by a full boundary-value solve.

## 1. Action, normal and dimensions

Use signature (−++++), two identical copies of the interior z<0, and the normal from each interior toward the shell. On the chosen copy n=e^(−B)∂z. The full doubled action is

\[
S=\sum_{a=1}^{2}\left\{\frac1{2\kappa_5^2}\int_{M_a}\sqrt{-g}
[R-(\nabla\phi)^2-2U(\phi)]\,d^5x
+\frac1{\kappa_5^2}\int_\Sigma\sqrt{-h}\,K_a\,d^4x\right\}
+\int_\Sigma\sqrt{-h}\left[-\frac{\sigma(\phi)}{\kappa_5^2}+\mathcal L_m\right]d^4x.
\]

There is one shell action, not two. The extrinsic curvature is Kμν=hμ^A hν^B∇A nB. The physical shell tension is λ=σ/κ₅². In the mass-dimension convention ℏ=c=1,

| Quantity | Mass dimension |
|---|---:|
| κ₅, κ₅² | −3/2, −3 |
| φ, Φ=φ/κ₅ | 0, 3/2 |
| U, σ, λ | 2, 1, 4 |
| canonical shell scalar χ | 1 |
| gb, ḡ=gb/κ₅ | −1/2, 1 |
| shell density ρ, scalar source j | 4, 4 |

Add one real shell scalar and, optionally, a light Dirac fermion:

\[
\mathcal L_m=-\tfrac12(D\chi)^2-\tfrac12m_\chi^2(\phi)\chi^2
+i\bar\psi\gamma^\mu D_\mu\psi-m_\psi\bar\psi\psi-y\chi\bar\psi\psi,
\]

\[
m_\chi^2=m_0^2+g_b^2(\Phi_b-\Phi_*)^2
=m_0^2+\bar g^2(\phi_b-\phi_*)^2.
\]

The fermion and Yukawa coupling are optional new physics; y is dimensionless. Interactions sufficient to thermalize its decay products are additional specifications. No value of a mass, gravitational scale, coupling, crossing point or cutoff is selected here.

For a classical field, define \(T_{\mu\nu}=-2(\sqrt{-h})^{-1}\delta S_m/\delta h^{\mu\nu}\) and j=−∂Lm/∂φ. The mass coupling gives

\[
j_\chi=\bar g^2(\phi_b-\phi_*)\chi^2.
\]

For produced quantum matter use the **renormalized** expectation values of Tμν and j, with one covariant renormalization prescription for both. The bare coincident quantity ⟨χ²⟩ is divergent. Vacuum polarization also generates local shell counterterms, including a φ-dependent tension and, generally, curvature operators. Retaining the registered σ as a renormalized input therefore requires explicit matching conditions; it does not follow automatically from adding the interaction. The equations below are the minimal action's classical or consistently renormalized truncation, with any neglected higher-derivative terms required to remain small.

## 2. Junctions obtained from the same action

Let Sμν=−λhμν+Tμν. The two Gibbons–Hawking–York variations and the single shell variation yield

\[
K_{\mu\nu}-Kh_{\mu\nu}=\frac{\kappa_5^2}{2}S_{\mu\nu},\qquad
n\phi=-\frac12[\sigma_{,\phi}+\kappa_5^2j].
\]

For the scalar, the boundary coefficient is −2nφ/κ₅²−σ′/κ₅²−j. Setting it to zero fixes the sign and the factor of two. Reversing the normal reverses the extrinsic curvature and normal scalar derivative; importing a formula with another normal convention without this conversion gives an incorrect sign. Standard scalar-brane derivations provide a comparison, but the expressions here are derived for this action and orientation. [Maeda and Wands, equations 2.15–2.19](https://arxiv.org/pdf/hep-th/0008188).

For a homogeneous isotropic shell with proper time τ, matter energy density ρ and pressure p, put H=dAb/dτ and v=dφb/dτ. In the conformal chart of the existing solver:

\[
\boxed{nA=\frac{\sigma+\kappa_5^2\rho}{6},\qquad
nB=\frac{\sigma-\kappa_5^2(2\rho+3p)}{6},\qquad
n\phi=-\frac{\sigma'+\kappa_5^2j}{2}.}
\]

The first two conditions are different once matter is present. Setting ρ=p=j=0 recovers all three registered boundary conditions. Multiplying these equations by eᴮ gives the coordinate Neumann data required by the existing conformal evolution equations. The bulk equations remain unchanged away from the shell.

## 3. Energy transfer has a fixed sign

The matter field equations imply the Ward identity

\[
D_\mu T^\mu{}_{\nu}=-jD_\nu\phi,\qquad
\boxed{\dot\rho+3H(\rho+p)=jv.}
\]

A growing particle mass can increase matter energy: if φ moves away from the mass minimum, jv>0. For a homogeneous classical χ this follows directly from χ̈+3Hχ̇+mχ²χ=0, ρχ=(χ̇²+mχ²χ²)/2 and pχ=(χ̇²−mχ²χ²)/2. The verifier checks that identity independently of the junction algebra. It also holds for the full covariantly renormalized quantum stress tensor and source, under the stated prescription; an arbitrary vacuum subtraction would not guarantee it.

Including tension, the signed shell budget is

\[
\frac{d(\lambda+\rho)}{d\tau}+3H(\rho+p)
=\left(\frac{\sigma'}{\kappa_5^2}+j\right)v
=-\frac{2}{\kappa_5^2}(n\phi)v=-2T^{\rm bulk}_{n\tau}.
\]

This last equality is also the Codazzi identity with the stated outward normals. The factor two represents the two sides. It is not permission to count bulk energy twice in a particle-production estimate.

After a particle description becomes valid, an illustrative decay ledger is

\[
\dot\rho_\chi+3H(\rho_\chi+p_\chi)=j_\chi v-Q,\qquad
\dot\rho_r+4H\rho_r=Q.
\]

The internal transfer Q cancels in the total budget. The action, quantum state and collision integrals must determine it. A decay approximation Q≈Γχρχ needs a nonrelativistic population and controlled decays; it is not an exact closure during the production event. Nonzero leakage would require its own flux term and consistent bulk source. Particle number, radiation energy, equilibrium and visible-sector identity are different claims.

## 4. Exact shell geometry and the unclosed bulk contribution

Define kₛ=nA, k₀=nB, w=nφ and the Weyl scalar 𝒲=−Eττ/3, where Eμν is the electric projection of the five-dimensional Weyl tensor. Spatial curvature is zero, as in the registered model. Direct projection of the bulk stress and Gauss equation gives

\[
\boxed{H^2=\frac{(\sigma+\kappa_5^2\rho)^2}{36}
+\frac{v^2}{12}-\frac{(\sigma'+\kappa_5^2j)^2}{48}
+\frac U6+\mathcal W,}
\]

\[
\boxed{\dot H=-\frac{\kappa_5^2}{12}(\sigma+\kappa_5^2\rho)(\rho+p)
-\frac{v^2}{3}-2\mathcal W,}
\]

\[
\frac{\ddot a}{a}=k_s k_0-\frac{v^2}{4}-\frac{w^2}{12}+\frac U6-\mathcal W.
\]

An intermediate check is that the projected scalar stress contributes Fττ=v²/4−w²/4+U/2 and Fij=(5v²/12+w²/4−U/2)hij. The extrinsic contributions are 3kₛ² and −kₛ²−2kₛk₀ respectively. The trace independently gives 6(Ḣ+2H²)=−v²−w²+2U+6kₛ(k₀+kₛ).

Differentiating the Friedmann identity, using energy exchange and the junctions, yields

\[
\dot{\mathcal W}+4H\mathcal W=
\frac16[4k_swv-4Hv^2-v\dot v+w\dot w-U'v].
\]

This is a consistency identity, **not** an autonomous four-dimensional model: v̇ depends on the normal bulk scalar evolution, and ẇ contains the evolving matter source. Solving only a radiation equation beside the old vacuum shell trajectory cannot supply those data. There is no justification for simply setting 𝒲=0 during a nonlinear event. Its sign is also not fixed by the local equations. Calling 𝒲 “radiation” would not establish a population of matter particles.

## 5. The exact constant-scalar endpoint fails one junction

In the registered dimensionless variables,

\[
W=1-\phi+\phi^3/3,\quad U=\tfrac12W'^2-\tfrac23W^2,
\quad\sigma=2W+\delta(1+c\phi),
\quad c=2/1.0357712571566784-4/3.
\]

At φ=1, U′=0, U=−2/27, U″=28/9, but σ′=δc. A bulk scalar identically equal to one has nφ=0. Therefore an exactly constant φ=1 **does not solve the pure-tension scalar boundary condition for the registered δ=.001**. It would require j=−δc/κ₅², or a changed tension, or a nonconstant bulk scalar. This rules out that precise ansatz, not a nearby de Sitter configuration or the measured tendency toward φ≈1.

A regular nearby static branch can be sought with φ=1+η. Around the AdS+ geometry ρ(y)=sinh(ky)/k, k=1/9, the linear scalar equation is

\[
\eta''+4k\coth(ky)\eta'-\frac{28}{9}\eta=0.
\]

Its regular solution is exactly

\[
\eta(y)=\eta_h f(y),\qquad
f(y)=\frac{C_{14}^{(2)}(\cosh ky)}{C_{14}^{(2)}(1)},
\]

because C₁₄⁽²⁾ satisfies (x²−1)fxx+5xfx−252f=0, and 252=14×18. Its logarithmic derivative tends to 14/9 at large radius. Matching to η′b=−2ηb−δc/2 at leading small δ gives

\[
\boxed{\eta_b=-\frac{9\delta c}{64}+O(\delta^2).}
\]

For a static warped maximally symmetric solution the exact radial Einstein constraint and junction give

\[
H_b^2=\frac{\sigma_b^2}{36}-\frac{\sigma_b'^2}{48}+\frac{U_b}{6}.
\]

Combining the matching estimate with this identity yields

\[
\boxed{H_b^2=\frac{\delta(1+c)}{27}
+\delta^2\left[\frac{(1+c)^2}{36}-\frac{c^2}{384}\right]+O(\delta^3).}
\]

The often-used constant-φ Randall–Sundrum comparison omits the negative δ²c²/384 correction and is not an exact solution of all the registered junctions. The series above is a candidate expansion, not an existence, uniqueness or stability proof. It also does not show that a time-dependent roll-off reaches the candidate.

For δ=.001 the finite-curvature *linear* matching, evaluated on the RS comparison geometry, gives the following starting values for an independent nonlinear boundary-value solve:

| Quantity | Seed estimate in registered units |
|---|---:|
| comparison ρb | 129.9237404812555 |
| comparison yb | 30.2766097392882 |
| regular radial amplification f(yb) | 6.28778836751138×10¹⁸ |
| ηb using the finite-curvature logarithmic derivative | −8.40485357060633×10⁻⁵ |
| ηh | −1.33669472942723×10⁻²³ |
| H² through second order in δ | 5.92401502678141×10⁻⁵ |
| constant-φ RS comparison H² | 5.92410802670494×10⁻⁵ |

Numerical work must evolve η=φ−1 directly: storing the cone value as 1+ηh in ordinary double precision rounds it to one and destroys this seed. The exact polynomial, its normalization and the H² series are checked by the companion verifier. No physical length unit is inferred from these numbers.

The separate full static boundary-value calculation in `../static_branch/PLUS_BRANCH_RESULTS.json` has now found a floating-point candidate at δ=.001:

\[
\phi_b=0.9999159473169134,\quad
\rho_b=129.92476284968,\quad H_b^2=5.9240147943293\times10^{-5}.
\]

The refined solve satisfies both junction residuals at approximately 10⁻¹⁷ and differs from the displayed second-order H² expansion by −2.3245×10⁻¹². These are numerical residuals, not rigorous error bounds. See the independent static-branch check in that directory for a second radial formulation. The candidate's perturbative stability and dynamical attraction remain unproved; it should not be substituted for a completed late-time evolution.

## 6. Production is testable without assuming oscillations

For one real χ initially in the appropriate adiabatic vacuum, a short approximately linear crossing of the mass minimum gives the following **after-crossing asymptotic occupation**, when the particle interpretation has again become valid:

\[
q=|g_b\dot\Phi_*|=|\bar g\dot\phi_*|,\qquad
n_k=\exp[-\pi(k_{\rm phys}^2+m_0^2)/q],\qquad
n_\chi=\frac{q^{3/2}}{8\pi^3}e^{-\pi m_0^2/q}.
\]

This is the established instant-preheating approximation, not new production physics. [Felder, Kofman and Linde, equations 3–4](https://arxiv.org/pdf/hep-ph/9812289). The verifier independently integrates this Gaussian over three-dimensional momentum. For m₀=0 it also checks the formal integral ∫d³k k nₖ/(2π)³=q²/(4π⁴). **This massless weighting of an asymptotic occupation is an energy proxy, not the instantaneous renormalized stress at the nonadiabatic crossing.** Establishing energy during production requires the mode functions and a consistent stress-tensor subtraction. After production, the particle energy must instead be evaluated with its actual evolving frequency.

Necessary approximation checks include |H|/√q≪1, |v̇|/(|v|√q)≪1, |Ḣ|/q≪1, q>0, and √q,mχ below the specified effective-theory cutoff. The initial quantum state, renormalization and backreaction must also be controlled. A trajectory which never crosses φ* does not trigger this approximation. Choosing φ* to fit a successful numerical outcome is a new parameter selection, not a prediction of the registered theory.

If production and decay occupy separated valid regimes and particles become nonrelativistic before decay, their density is approximately mχ,d nχ,* (a*/ad)³ until decay. The mass-growth work is paid for by the jv term. It cannot be added as free energy. Decay acting during the nonadiabatic event can modify production itself, so the undamped occupation formula cannot automatically be combined with instantaneous conversion. Thermal radiation after decay further requires sufficiently rapid scattering; a formal temperature (30ρr/π²g*)^(1/4) is physically a temperature only once equilibration has been demonstrated.

## 7. A positive residual vacuum imposes a quantitative radiation test

As a **conditional benchmark**, freeze the scalar and the matched renormalized tension at σf>0, neglect Weyl corrections and continuing exchange, and let a consistent vacuum solution have Hvac²>0. Then

\[
H^2=H_{\rm vac}^2+\frac{\kappa_5^2\sigma_f}{18}\rho_r
+\frac{\kappa_5^4}{36}\rho_r^2.
\]

The radiation contribution exceeds the residual vacuum term only if

\[
\boxed{\rho_r>\rho_{\rm crit}
=\frac{\sqrt{\sigma_f^2+36H_{\rm vac}^2}-\sigma_f}{\kappa_5^2}.}
\]

The exact radiation deceleration condition within this same benchmark is slightly different:

\[
\rho_r>\frac{\sqrt{\sigma_f^2+108H_{\rm vac}^2}-\sigma_f}{3\kappa_5^2}.
\]

At low energy κ₅²ρr≪σf, κ₄²=κ₅²σf/6 and ρcrit≈3Hvac²/κ₄². Radiation-dominated expansion lasting N e-folds after its last injection requires

\[
\boxed{\rho_{r,d}>e^{4N}\rho_{\rm crit},\qquad
N_{\rm max}=\tfrac14\log(\rho_{r,d}/\rho_{\rm crit})}
\]

when that ratio is above one. A positive fixed vacuum eventually wins over freely redshifting radiation. This is a restriction on this specified post-injection benchmark, not a no-go theorem for continued production, a changing vacuum, or the full dynamical bulk. Merely creating some particles does not establish a long hot epoch.

No physical scale is needed to formulate a parameter scan. If the physical coordinate unit is L, define α=ḡL, β=κ₅²/L³, μ₀=m₀L, q̂=qL²=|α dφ/dτ̂|, σ̂=Lσ and R=κ₅²Lρ. The formal massless-weighted occupation proxy and the geometric threshold are

\[
R_{\rm proxy}=\frac{\beta\hat q^2}{4\pi^4},\qquad
R_{\rm crit}=\sqrt{\hat\sigma_f^2+36\hat H_{\rm vac}^2}-\hat\sigma_f.
\]

If one *postulates* a toy deposited budget εRproxy with efficiency ε≤1, its radiation-duration criterion is εβq̂²/(4π⁴)>e^(4N)Rcrit. This is algebra conditional on that postulate, **not a validated production or prompt-transfer bound**. A physical transfer estimate must establish a regime separating valid particle production from decay, use the actual particle energies there, and include any mass-growth work. If decay overlaps production, the coupled process must be solved instead. A scan over (α,β,μ₀,y,φ*) would state these added assumptions explicitly; the existing pure-tension simulation fixes none of them.

### A deliberately restricted four-dimensional budget illustration

Assume, additionally, a closed low-energy four-dimensional system with a **single constant** reduced Planck mass M₄, initial total density 3M₄²H₀², the same constant residual vacuum afterward, and no subsequent external injection. Even converting 100% of all remaining initial energy to radiation gives

\[
\rho_{r,d}/\rho_{\rm vac}\le(1-f_\Lambda)/f_\Lambda,
\qquad f_\Lambda=H_{\rm vac}^2/H_0^2.
\]

Inserting the independently calculated candidate H² and the archived initial ρb=78.8281771422 gives fΛ=0.36811126, radiation/vacuum≤1.71657 and N≤0.135082 e-folds. These numbers are a **toy-budget consequence, not a bound on the actual registered five-dimensional system**. In fact, σ changes from 2.00095445 to 0.66826423 between the configurations: the local low-energy identification κ₄²=κ₅²σ/6 changes by a factor 0.333973. Thus even the fixed-M₄ premise cannot simply be borrowed for this transition. The bulk energy and its exchange are also unclosed. This illustration shows why a claim of a long hot era must identify its residual-vacuum depletion or additional controlled energy budget; it does not rule out one in a different, consistently derived effective description.

### An actual trajectory-based eligibility calculation

The companion `analyze_crossing_eligibility.py` reads the unchanged, completed external Chat14 archive. It performs no new PDE evolution or particle creation. For hypothetical crossings at five values of φ*, it estimates proper-time derivatives using quartic local fits with half-widths Δ(H₀τ)=.05 and .10. Three archived grids are compared: default, finer shell spacing, and finer far spacing. All listed crossings occur at coordinate t<6.65, before the earliest far-taper signal at t=8.5. The known initial junction mismatch and remaining constraint limitations of those trajectories still apply.

Define G=|gb/(κ₅H₀)|=|ḡ/H₀|, s=H₀τ, v₀=dφ/ds, h=H/H₀, q/H₀²=G|v₀|. Imposing three local approximation criteria at tolerance ε=.1 requires

\[
G\ge\max\left\{
\frac{h^2}{\epsilon^2|v_0|},\;
\frac{(dv_0/ds)^2}{\epsilon^2|v_0|^3},\;
\frac{|dh/ds|}{\epsilon|v_0|}\right\}.
\]

| Hypothetical φ* | Crossing H₀τ, finer far grid | Minimum G, finer far grid | Relative spread over three grids |
|---|---:|---:|---:|
| .10 | 4.446076 | 1414.30 | 0.0076% |
| .25 | 5.048908 | 397.61 | 0.0147% |
| .50 | 5.575005 | 136.22 | 0.0077% |
| .75 | 5.993831 | 110.69 | 0.0267% |
| .90 | 6.314273 | 1488.50 | 0.0729% |

The derivatives here come from the saved scalar trajectory, not by silently equating the archived dissipative velocity variable to its exact derivative. On the finer far grid, fitted and stored velocities agree to within 1.7×10⁻⁵ relative at these crossings; changing the derivative-fit window changes Gmin by at most 7.7×10⁻⁵ relative. Agreement among these estimators is useful but is not a certified continuum error bound.

Within this **fixed-trajectory local screen**, a mid-roll crossing is less demanding than a crossing near either slow part of the roll. This is an actionable parameter-screening result. It does not choose G, prove that the coupling lies below the effective-theory cutoff, set m₀, test backreaction, or show successful reheating. A large G relative to a small H₀ does not by itself imply a large dimensionless quantum interaction; that question requires the missing physical normalization and cutoff. The exact thresholds are recorded in `CROSSING_ELIGIBILITY_AND_TOY_BUDGET.json`, together with all three grid estimates and both derivative windows.

## 8. Current primary-source screen and the next calculation

The web screen on 22 September 2026 covered brane scalar reheating, bulk-inflaton energy loss, instant production and recent braneworld warm inflation. It cannot cover the entire web or establish absence of prior art. No result dated after the research date is relied upon.

- [Maeda and Wands, *Dilaton-gravity on the brane*](https://arxiv.org/abs/hep-th/0008188), 2000: action-based scalar/matter junctions and bulk contributions provide a convention check. Their normal and scalar normalization must be converted, not copied unchanged.
- [Himemoto and Tanaka, *Braneworld reheating in the bulk inflaton model*](https://arxiv.org/abs/gr-qc/0212114), December 2002, and [Tanaka and Himemoto, *Generation of dark radiation in the bulk inflaton model*](https://arxiv.org/abs/gr-qc/0301010), January 2003: primary precedents for separating brane energy deposition from bulk loss. Their dissipation model does not fix the present coupling or loss rate.
- [Yeasmin and Deshamukhya, *Warm inflation in a braneworld scenario*](https://link.springer.com/article/10.1140/epjc/s10052-026-16156-3), journal publication 27 August 2026, [preprint first posted 10 December 2025](https://arxiv.org/abs/2512.09389): a relevant recent consistency study. Its treatment assumes a nearly equilibrated, quasi-stationary bath and slow roll; its fitted fluctuation function does not establish cold-start production on this rolling shell.

The next coupled calculation should first establish a reliable pure-tension trajectory and a valid physical normalization, then solve brane χ mode functions with a covariantly renormalized ⟨Tμν⟩ and ⟨j⟩, feed **all three** altered junctions back into the bulk evolution, and track the signed energy residual. Only after a nonzero controlled population exists should the decay and scattering calculation be used to assess radiation domination. The present work makes that program concrete without pretending it has already been carried out.

## Verification and scope

`verify_matter_extension.py` checks 47 exact identities and dimensional relations plus four deliberately nonvanishing wrong-sign/factor controls. It writes `MATTER_EXTENSION_EXACT_CHECKS.json`. Tested with SymPy 1.14.0 using explicit exceptions, including under optimized Python. The checks cover the boundary factors, signed energy transfer, projected Einstein identities, Weyl balance, Gaussian integrals, positive-vacuum thresholds, registered endpoint obstruction and the regular Gegenbauer solution. They do not prove existence of a full static branch, solve the quantum state, evolve the coupled model, or establish external mathematical novelty.

`analyze_crossing_eligibility.py` needs NumPy and the included `source_audit/inputs/CHAT14_COMPLETED_REFERENCE.zip` and `static_branch/PLUS_BRANCH_RESULTS.json`. Four extraction checks pass. Its output records the hashes of both numerical sources, the deliberately strong four-dimensional toy assumptions, and the approximation criteria. It is an analysis of existing data rather than an independent certification of their evolution.
