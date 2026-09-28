# Independent audit of the completed external Chat14 checkpoint

22 September 2026. Read-only audit of the external Claude-produced checkpoint in `Documents/D-Blast 3/untitled folder 151`. This report distinguishes finite numerical evidence from endpoint claims. The source documents, including their proposed next steps, were treated as reference material rather than instructions. No original file was edited and no long PDE evolution was rerun.

## Conclusion

The newly completed archive supplies registered-detuning nonlinear evidence that was absent from our earlier snapshot. It supports a positive-side roll toward a lower expansion rate, and a negative-side turnaround followed by substantial contraction. It does **not** establish the final asymptotic positive state, a singularity theorem, or a completed constraint-controlled evolution from boundary-compatible initial data.

A separate analytical correction matters: the source's “exact RS brane in the φ=+1 vacuum” satisfies a metric-junction benchmark but **fails the scalar junction condition of the registered model**. A nearby de Sitter solution with a weak scalar profile may exist. Constant φ=1 is not that solution.

## Provenance and chronology

- The completed ZIP has SHA256 `c68e1de1ebb96cffea4945699b17b4c3689bffbbd8eca6630a3d13006bfdf59d`.
- All **84** manifest payloads occur in the ZIP and match their SHA256 values; the ZIP CRC test passes. Their total size matches the manifest's 4,247,659 bytes. The loose source-folder copies also match all manifest hashes.
- `inputs/CHAT14_COMPLETED_REFERENCE.zip` is an unchanged copy. `SOURCE_PROVENANCE.json` records original locations, hashes, and file modification times. The source package also contains a wrapper README outside its inner manifest; it is included in the unchanged ZIP but not counted among the 84 payloads.
- Our prior `CHAT8_CONTINUATION_20260922` snapshot contained unfinished runs. A file comparison shows the final `rolloff5d_v2.py` changes its earlier frozen copy only by adding explanatory documentation, a configurable grid-transition width with the same default, and maximum-constraint-location logging. The fourth-order closure and differentiated boundary slopes were already present in that snapshot. The new information is principally the completed data and referee work, not a newly changed continuum theory.
- Referee reports are dated records of intermediate review states. For example, the numerics referee says no negative-side refined-far-grid run exists; the final archive does contain one. We use the completed data when it supersedes that statement.

The replay script `recompute_chat14.py` reads the reference ZIP without importing its evolution code. It analyzes **17** summary/time-series pairs and passes **seven** integrity and numerical-consistency checks. These checks do not certify the continuum physics.

## What the saved registered runs support

The independent event calculations use linear interpolation between adjacent samples. Times below are dimensionless shell proper time H₀τ, rather than coordinate time t. This is essential because the lapse differs between grids.

| Quantity | Default shell/far grid | Refined shell grid | Refined far grid |
|---|---:|---:|---:|
| Negative branch φ=-1 crossing time | 5.608477 | 5.608284 | 5.608507 |
| H/H₀ at that crossing | 0.797771 | 0.797756 | 0.797775 |
| Turnaround time H=0 | 5.894718 | 5.894399 | 5.894724 |
| φ at turnaround | -1.941666 | -1.941217 | -1.941593 |
| Time at H/H₀=-5 | 6.199854 | 6.198035 | 6.199539 |
| Time at H/H₀=-10 | 6.282911 | 6.276327 | 6.284781 |

The turnaround is reproducible across these changes. The default-grid value of U at turnaround is **+3.668909** in the model's units, while the detuning contribution to the tension is **-0.000160330**. Therefore negative U is not the trigger at this event. A causal interpretation still requires the full dynamical equations; this sign check alone is not a proof of mechanism.

Agreement in the *time of a threshold crossing* is not the same as agreement in H at the same time. At the proper time when the default run has H/H₀=-5, the refined-shell run has -5.07098, a 1.42% discrepancy. At the default -10 event the refined-shell value is -10.95234, a **9.52%** discrepancy. At the default -22 event it is -28.78662, a 30.85% discrepancy. These explicit comparisons are somewhat larger than the source's approximate 8% and 23% figures, but lead to the same practical conclusion: do not label the late contraction converged.

For the positive branch, at the same proper time H₀τ=6.9:

- Default far spacing 0.008 gives H/H₀≈0.61106.
- Refined far spacing 0.002 gives **0.638641925**, with φ=0.993824942 and coordinate t=8.359801.

The two shell resolutions sharing far spacing 0.008 do not test this discrepancy. Only one far-region refinement is available. In addition, changing far spacing changes the steepness of the grid transition, so the comparison does not isolate far-cell dispersion from transition reflection. A third far spacing and an independent transition-width control remain useful.

The source's linear-growth controls are relevant and successful within their stated small-amplitude windows. They do not by themselves establish nonlinear constraint control or late-time convergence. Our earlier solver's separate linear calibration remains complementary evidence.

## Domain and initial-data limitations

All decisive completed registered runs use L=10, with the initial perturbation tapered over the outer 15% of the interval. The taper begins at z=-8.5. It changes the initial constraints there; the earliest continuum characteristic from that location can reach the shell at **t=8.5**. The far boundary itself can first affect the shell near t=10. Discrete schemes need separate boundary-convergence checks; these continuum times are not a numerical error bound.

For the refined positive run, at t=8.5 the saved data interpolate to:

| H₀τ | φ_b | H/H₀ |
|---:|---:|---:|
| 6.920649612 | 0.994535176 | 0.637140760 |

The final saved point is t=10.996960, H₀τ=7.052603, H/H₀=0.628690, lapse ratio 0.012328. It is **after both causal thresholds**, with no archived domain-extension comparison. Hence 0.628690 should not be promoted to a clean endpoint. The 6.9 comparison occurs before the taper threshold and is the stronger quantitative datum, though it still lacks a third far-grid refinement. The negative-side turnaround and quoted H=-5 events occur well before t=8.5.

The initial data come from a neighboring static shell with c+Δc, then evolve with the original c. At the shell, relative to the target boundary equations, the analytic residuals before discretization are

\[
\Delta A_z=\Delta B_z=e^B\delta\Delta c\,\phi_b/6,
\qquad
\Delta\phi_z=-e^B\delta\Delta c/2.
\]

They do not vanish for the nonlinear seeds Δc=±10⁻⁴. The source recognizes the outgoing constraint pulse caused by the initial mismatch. This is not equivalent to preparing exactly compatible initial data. The initial-data construction in our previous checkpoint addresses a distinct remaining requirement; its effect on nonlinear evolution has not been determined by these external runs.

## Constraint and discretization qualifications

The archive monitors the normalized **momentum** constraint. It does not save a Hamiltonian-constraint history. Its normalization is a sum of background-sized gradient and velocity terms; a small displayed ratio is not automatically small compared with the perturbation or a bound on an observable error.

At the exact shell node, substituting the junction and its time derivative makes the momentum constraint an algebraic identity. The maximum over the last six cells contains independent near-shell information, but the shell-node value itself is not an independent boundary test. An independent one-sided derivative and both constraints are preferable.

The final NPZ files store final a, b, f and the static background plus shell time series. They do **not** store final full-domain velocities or intermediate full-state snapshots. The archive therefore cannot support a fresh full Hamiltonian/momentum or geometric five-dimensional curvature audit of the nonlinear interior without rerunning it.

The optional closure4 is a valid improvement over the older even-ghost closure; its second derivative is third-order at the shell, despite the shorthand name. The following previously identified issues remain:

- Dissipation is normalized per grid index, not by local physical spacing, and the final code applies it to fields as well as velocities. Thus stored velocity variables differ from the exact time derivatives of the dissipative discrete fields. The small coefficient does not make this identity exact. No archived explicit zero-dissipation comparison closes that concern for the nonlinear endpoint.
- The third dissipative ghost value is linearly extrapolated. It does not reproduce quadratic data exactly. This is separate from the correctly derived closure for the main derivative stencil.
- Both the seed mismatch and far taper remain. The earlier correction to time-derivative boundary slopes fixes a diagnostic problem, but does not repair the seed.
- “No grid noise” and a constant-coefficient stability test do not establish nonlinear Einstein-constraint preservation.

## Scalar-junction correction to the endpoint benchmark

The registered functions are

\[
W=1-\phi+\phi^3/3,\quad U=\tfrac12(W')^2-\tfrac23W^2,
\quad \sigma=2W+\delta(1+c\phi),
\quad c=0.5975949350280132.
\]

For a constant φ=1 field, its normal derivative is zero. The actual scalar junction requires

\[
\partial_n\phi=-\sigma'(1)/2=-\delta c/2.
\]

At δ=0.001 this is **-0.000298797467514**, not zero. The same derivative issue occurs at φ=-1. Consequently, computing

\[
H_{\rm metric}^2=[\sigma(1)/6]^2-k_+^2,\qquad k_+=1/9,
\]

and obtaining H_metric/H₀=0.606726506 is a useful metric benchmark, **not an exact solution of all field and junction equations**. The machine report explicitly labels it this way.

An independently checked leading-order possibility is a nearby scalar-profile endpoint. Since U″(1)=28/9, the regular outward scalar exponent in the zero-curvature AdS+ limit solves

\[
\lambda^2+4k_+\lambda-U''(1)=0,
\qquad\lambda_+=14/9.
\]

Writing φ_b=1+η_b, its leading scalar junction is λ₊η_b=-2η_b-δc/2. Thus

\[
\eta_b=-9\delta c/64+O(\delta^2),
\qquad \phi_b\simeq0.9999159632.
\]

This expansion assumes a regular weak-profile branch continuous from AdS+ and treats curvature corrections consistently as higher order. It is a candidate for a full boundary-value solve, not an existence/stability proof or the demonstrated endpoint of the saved evolution. The final saved positive φ≈0.99756 has not reached it.

## Curvature and physical interpretation

The replay optionally records a finite on-shell inference

\[
R_5/H_0^2=-[d\phi/d(H_0\tau)]^2+
\rho_b^2\{[\sigma'(\phi)/2]^2+(10/3)U(\phi)\}.
\]

It follows from the Einstein trace and scalar junction, so it is not an independent curvature test. Large finite negative values in the late contraction data suggest rapidly growing curvature only under those equations and the numerical solution's validity. A finite time series whose resolutions diverge cannot establish “grows without bound” or a final singularity.

Likewise, the simulated theory contains no explicit visible matter field. It cannot demonstrate a thermal radiation bath; that is a limitation of the model's ingredients, not a general no-go theorem for matter production or a higher-dimensional origin. The positive trajectory need not remain at its terminal recorded H value after the chart changes. Use “roll toward a lower-expansion configuration” and “turnaround followed by contraction” for what the present evidence supports.

## Reproduction

Run `python recompute_chat14.py` in this directory with NumPy available. It only reads the copied ZIP and writes `CHAT14_RECOMPUTED_RESULTS.json`. It verifies the archive, reproduces event comparisons and causal thresholds, and labels assumptions on inferred quantities. It does not execute either archived PDE solver. The independent assertions are checks on the saved evidence, not a replacement for new constraint-compatible, domain-controlled nonlinear evolutions.
