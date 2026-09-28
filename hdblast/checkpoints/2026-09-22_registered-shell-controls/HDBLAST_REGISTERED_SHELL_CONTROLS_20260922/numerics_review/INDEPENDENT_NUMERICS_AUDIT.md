# Independent audit of the nonlinear HDBLAST continuation

22 September 2026. Scope: Chat13 archived solver, mathematical constraint subsystem, and the unfinished mapped-grid follow-up in folder151. No original files were edited and no long evolution was run for this audit.

## Finding

The three nonlinear wave equations and their exact-background subtraction are internally consistent with the stated Einstein–scalar model. The existing registered-detuning runs do **not** yet establish the nonlinear fate at detuning 0.001. Their outgoing cone-region constraint front has a specific continuum explanation that must be separated from the physical unstable shell mode: the constraint equations themselves transport and strongly amplify unweighted constraint errors away from the shell. This is an explanation of a possible error mechanism, not a retrospective proof of the cause of every failed run.

I independently differentiated the constraints, substituted the evolution equations, and verified the characteristic identities exactly with SymPy. Seven positive checks pass and six deliberately incorrect variants are rejected. A second script performs ten archive/stencil controls. These are newly executed checks, not copied historic verification counts.

## 1. Exact constraint transport

Use the source chart and conventions:

\[
ds^2=e^{2B}(-dt^2+dz^2)+e^{2A}d\mathbf{x}_3^2,\qquad z\leq0.
\]

The solver's Hamiltonian and momentum residuals are

\[
\begin{aligned}
H={}&-2e^{2B}U+6A_t^2+6A_tB_t-12A_z^2+6A_zB_z
      -6A_{zz}-\phi_t^2-\phi_z^2,\\
M={}&-3A_{tz}-3A_tA_z+3A_tB_z+3A_zB_t-\phi_t\phi_z.
\end{aligned}
\]

On the three evolution equations, the remaining Einstein residual components are `E_tt=E_zz=H/2`, `E_tz=M`, and `E_ij=0`. The scalar equation supplies stress conservation. Direct differentiation, equivalently the contracted Bianchi identity, gives

\[
H_t=2M_z+6A_zM-3A_tH,
\qquad
M_t=\frac12H_z+\frac32A_zH-3A_tM.
\]

Consequently, with `C_+=H+2M` and `C_-=H-2M`,

\[
\boxed{(\partial_t-\partial_z)[e^{3A}C_+]=0},\qquad
\boxed{(\partial_t+\partial_z)[e^{3A}C_-]=0}.
\]

`C_+` propagates toward decreasing `z` along `z=z_0-(t-t_0)`; `C_-` propagates toward the shell along `z=z_0+(t-t_0)`. The factor `e^{3A}` is the volume density of the homogeneous three-dimensional slices. It is a useful diagnostic weight; it does not make a nonzero Einstein residual physically acceptable.

On the static background `A=t+ln rho(z)`, an outgoing error launched at the shell at `t=0` satisfies

\[
\frac{C_+(t,-t)}{C_+(0,0)}=
e^{-3t}\left[\frac{\rho_b}{\rho(-t)}\right]^3.
\]

Evaluated on the archived registered background (`rho_b=78.82817709590013`):

| Coordinate time | Outgoing position | Unweighted error amplification |
|---:|---:|---:|
| 0.1 | -0.1 | 108.78 |
| 0.5 | -0.5 | 3,016.32 |
| 1.0 | -1.0 | 7,345.33 |
| 2.0 | -2.0 | 10,575.57 |
| 2.5 | -2.5 | 10,938.76 |
| 3.0 | -3.0 | 11,074.43 |

This table uses floating-point archived background values and interpolation; it is not an interval certificate. Near the regular cone, `rho~C exp(z)`, so the amplification tends to a constant rather than continuing to grow exponentially forever. A front reported near `z=-2.5` at `t=2.5` follows exactly the outgoing constraint characteristic from the shell. This spatial coincidence and the large gain motivate a targeted diagnostic. They do not by themselves identify where the discrete error was generated, exclude simultaneous physical perturbations, or establish that correcting the front will fix the measured shell growth rate.

This derivation is a project-specific application of standard Bianchi constraint propagation. No claim of new mathematics or an additional physical unstable mode is warranted.

## 2. Concrete issues in the existing mapped-grid follow-up

The inspected file is `untitled folder 151/HDBLAST_CHAT14_REGISTERED_DETUNING_ROLLOFF_20260922/rolloff5d_v2.py`; it imports its copy of the Chat13 model/background. Existing registered logs end at early times and do not document a successful registered nonlinear run.

1. **The shell momentum diagnostic is invalid at its last points.** It uses a zero Neumann slope for `pa`, although the actual shell velocity derivative is generally nonzero. The exact data are

   \[
   \partial_z A_t=\frac{e^B}{6}(\sigma B_t+\sigma'\phi_t),
   \quad \partial_z B_t=\partial_z A_t,
   \quad \partial_z\phi_t=-\frac{e^B}{2}(\sigma'B_t+\sigma''\phi_t).
   \]

   Substituting these and the field junction data gives `M_b=0` identically. Imposing that identity inside the monitor would hide errors; a useful independent monitor differentiates stored velocities one-sidedly and compares them with these required slopes. Chat13's referee already used a one-sided monitor, but V2 regressed to the incorrect zero slope while reporting a shell residual.

2. **The automatic registered growth-fit window can be empty.** For the logged `dc=1e-4` seed, `|delta phi_b(0)|≈8.979e-5`. The code asks simultaneously for `|delta phi_b|>30*seed≈0.002694` and `<0.002`. No sample can pass. A null fit is therefore not evidence of missing growth. Reduce the seed or select explicit, reported linear windows after transients, and compare signed seeds and several windows.

3. **Mapped-grid dissipation is not normalized like Chat13.** The V2 grid-scale KO damping is approximately `epsilon` per coordinate time; Chat13's is `epsilon/dz`. To maintain the latter local physical-time convention one would use `epsilon/z'` for unit computational spacing. Either convention can be defined, but the same numerical `--ko` is not the same damping experiment. Stability and convergence must be tested after any change; adding damping is not a proof of correctness.

4. **V2's third KO ghost fails a simple polynomial control.** The first two ghosts are generated from the Neumann closure; the third is set to `2*ghost2-ghost1`. For `f(xi)=xi^2`, the correct values are `(1,4,9)`, while V2 uses `(1,4,7)`. The sixth difference at the last point becomes `-2` rather than zero. A consistent polynomial/Hermite third ghost removes this artificial boundary forcing.

5. **The original symmetric Neumann ghost closure is low order for the second derivative.** It imposes the first derivative exactly but gives `D2(z^3)|_0=-(4/3)h` instead of zero. Thus the “fourth-order finite differences” description applies to the centered interior, not the full boundary discretization. The previous measured global second-order behavior is consistent with this issue. A quintic Hermite closure fitted to five interior/shell values plus the shell first derivative reproduces degree-five polynomials and allows a fourth-order second derivative at the boundary, but stability still needs testing.

6. **Exact zero deviation is an inert numerical test.** A zero seed often makes every deviation term exactly zero and stays zero even with an unstable discretized perturbation operator. It tests the subtraction identity, not perturbative stability. Explicit small nonzero seeds or a discrete eigenvalue test are needed.

7. **The seed does not satisfy the target junction exactly.** It is a static bulk solution for `c+dc`, evolved with tension slope `c`. Before any taper, its interior continuum constraints hold; its shell derivatives differ from target by

   \[
   \Delta A_z=\Delta B_z=\frac{e^B\,t_{\rm det}\,dc\,\phi_b}{6},\qquad
   \Delta\phi_z=-\frac{e^B\,t_{\rm det}\,dc}{2}.
   \]

   It is therefore an initial-boundary corner perturbation, not fully compatible data for the unchanged shell theory. At registered detuning and `dc=1e-4`, the scalar derivative mismatch is about `-3.94e-6`. A spatial taper also destroys exact bulk constraints in its transition region. An outer taper is causally irrelevant to the near-shell continuum solution only for `t<(1-taper)L`; the finite-difference scheme still needs a domain-size check.

8. **The outgoing front cannot casually be identified with the clamped far-end seed.** A far-end error near `z=-L` reaches the shell no earlier than approximately `t=L` in this chart, or `t=(1-taper)L` from the nearest taper edge. The reported front at `z≈-t` is consistent with a shell-origin disturbance moving outward. Its causal direction needs to be measured, not inferred from the presence of a far-boundary clamp.

9. **The coordinate map does not change characteristic causality.** The chain-rule formulas in V2 are correct: `f_z=f_xi/z'` and `f_zz=f_xixi/(z')^2-f_xi*z''/(z')^3`. Characteristic speed in the computational coordinate is `±1/z'`; in physical conformal `z` it remains `±1`. The near-shell effective spacing is larger than nominal `dz_fine` because the tanh transition has not saturated at `z=0`. The output logs correctly report effective spacing. Refinement studies must use actual spacing.

10. **A fixed small fine-grid region is insufficient for a complete roll-off.** As the domain wall leaves the shell, it moves into coarser regions. A grid that passes the initial linear test need not resolve its later nonlinear motion. Further full-run claims require moving/extended refinement or a demonstrably adequate coarse region, and an independent time-slicing treatment when the brane lapse collapses.

## 3. Dissipation and monitoring qualifications

In Chat13/V2, KO is added to both field equations and velocity equations. Then the stored `pa` differs from the actual numerical derivative `a_t` by the field KO term. The same issue applies to `pb` and `pf`. Observables and constraints computed directly from stored velocities inherit this finite-resolution ambiguity. Applying KO only to the acceleration equations preserves the interpretation of stored velocities, but it still modifies the continuum constraint equations off-shell.

For a momentum-only modification

\[
A_{tt}=E_A+Q_A,\quad B_{tt}=E_B+Q_B,\quad\phi_{tt}=E_\phi+Q_\phi,
\]

the additional sources are

\[
\begin{aligned}
S_H&=(12A_t+6B_t)Q_A+6A_tQ_B-2\phi_tQ_\phi,\\
S_M&=-3\partial_zQ_A+3(B_z-A_z)Q_A+3A_zQ_B-\phi_zQ_\phi.
\end{aligned}
\]

Thus `(partial_t ∓ partial_z)[exp(3A) C_±]=exp(3A)(S_H±2S_M)`, in addition to truncation effects on the discrete grid. The exact homogeneous characteristic identity must not be reported as an exact property of a dissipated finite-difference run. Monitor both absolute and normalized `H,M`, independent shell derivative defects, and weighted outgoing/incoming residuals; retain the complete fields **and velocities** for retrospective checks.

## 4. Bounded validation recommended before any registered-fate claim

Use corrected boundary diagnostics and compatible high-order ghosts; choose a small explicit seed; run both seed signs over a short linear interval; use at least two true near-shell resolutions, two dissipation strengths, and an adequate causally separate outer boundary. Compare growth against the independently derived physical eigenvalue near `1.657` in these dimensionless time units. Check stability of fit windows, not merely one fitted number. Keep the expected physical tachyon: a solver that erases it through excessive damping has not passed.

Only extend to nonlinear outcomes after this gate passes with small, converging Hamiltonian and momentum residuals and a front analysis. A short linear-gate pass cannot establish a hot universe, heating, a late de Sitter limit, or the nature of a singularity. The matter-free homogeneous model does not contain a visible thermal sector, so its classical evolution alone cannot demonstrate visible reheating.

## 5. Reproduction and provenance

- `verify_constraint_transport.py`: independent SymPy derivation, seven positive and six negative controls; output `CONSTRAINT_TRANSPORT_CHECKS.json`.
- `check_archived_numerics.py`: read-only Chat13 archive diagnostic and exact stencil controls; ten checks; output `ARCHIVED_NUMERICS_CHECKS.json`. Takes an optional archive path.
- The original solver read from the ZIP matches the unpacked Chat13 solver byte for byte. Its SHA-256 and the archived numerical background's SHA-256 are recorded in the diagnostic output. The archive hash is recorded too.
- No archived evolution was rerun in this audit. Historic numerical claims above remain attributed to the archive; the exact identities and ten small diagnostic controls were newly verified.
