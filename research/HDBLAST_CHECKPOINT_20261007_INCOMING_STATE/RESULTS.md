# Conditional incoming-state stress and direct-work transfer

Ricardo Maldonado · 7 October 2026 UTC · mathematical methods checkpoint

The unchanged linear metric model now has an exact incoming-state counterpart
to the earlier source-error kernels. It separates conserved normalization drift
from a correlated homogeneous phase amplitude, bounds their direct continuous
finite-band density and pressure contributions, and gives exact state-difference
primitives for continuity work and time-integrated pressure.

**The actual incoming state remains unenclosed.** This is a conditional theorem
about independently supplied state-error envelopes. No source callback, saved
quantum-array decode, physical trajectory or twelve-case calculation occurs.
The complete numerical pressure/contact certificate remains **UNRESOLVED**,
metric calibration remains **FAIL**, and a higher-dimensional origin of the
Big Bang remains **NOT_ESTABLISHED**.

For identical forcing and contacts, define `c=Re deltaU_a-Im deltaW_a/(2k)`
and `A=deltaW_a/(2ik)`. The normalized first-order Wronskian coefficient is
`2c`, corresponding to physical drift `2epsilon c`. Exact
first-order normalization establishes `c=0` only if independently verified.
The phase evolves as `A exp[2ik(t-a)]`; retaining that correlation cancels
the density's leading phase term while preserving the pressure's leading term.

Under a genuinely proved uniform `|A(k)|<=sigma`, `c=0`, `L<=2/7` and the
same declared measure constant `Pi>=3`, the all-time finite-band bounds are

\[
|\delta R_K|\le\left(K^3/189+K/686\right)\sigma,
\qquad
|\delta P_K|\le\left(K^4/108+5K^2/1764\right)\sigma.
\]

The [full theorem](INCOMING_STATE_TRANSFER.md) includes arbitrary `c` bounds,
spectral/infrared hypotheses, the exact closed density norm integral, and
direct endpoint-transfer formulas. The source, residual, arithmetic, geometry,
contact, quadrature, serialization and ultraviolet terms remain separate.

For scale only, the following entries illustrate the theorem *if* an independent
calculation were to prove `sigma=10^-16` and `c=0`. That enclosure has not been
proved for this project's incoming state. These are outward rounded sensitivity
bounds, not achieved errors or new acceptance conditions.

| K | Density upper | Pressure upper | Continuity-work upper | Time-integrated pressure upper |
|---:|---:|---:|---:|---:|
| 64 | 1.38710 × 10^-13 | 1.55357 × 10^-11 | 2.46593 × 10^-13 | 3.23653 × 10^-13 |
| 128 | 1.10963 × 10^-12 | 2.48556 × 10^-10 | 1.97266 × 10^-12 | 2.58912 × 10^-12 |
| 256 | 8.87688 × 10^-12 | 3.97685 × 10^-9 | 1.57811 × 10^-11 | 2.07127 × 10^-11 |

Direct differentiation of the independently specified density and pressure
proves `deltaR'=L(deltaR-3deltaP)` and `(-deltaR/(3L))'=deltaP`. The smaller
integrated-pressure envelope uses that exact oscillatory primitive; it does
not bound the integral of the absolute pressure error. Conserved incorrect
states still obey both formulas. The verifiers also exhibit an erroneous
paired omission of canonical geometry terms that passes Ward conservation
while changing both direct stress operators.

Two separately implemented routes verify the algebra and rational bounds in
ordinary and optimized Python. The standard-library Laurent-polynomial route
passes 34 exact checks and rejects 13 mutations; the independent symbolic
route passes 54 exact checks and rejects 17 controls. All 24 rational
coefficient values agree across the two routes. The same source
bytes are used in both modes; no check relies on Python assertions. This is
internal independent AI-assisted verification, not external peer review, an assessment
of mathematical novelty or a proof-assistant formalization.

[Reproduction](REPRODUCE.md) authenticates the complete new payload, reruns
both routes in both modes, compares every scientific receipt field and checks
that the payload bytes remain unchanged. Pinned inherited documents identify
the model and old execution restrictions. This new result is separate from
the immutable 3 October publication packet and 5 October proof.

The next calculation should audit the actual incoming-state prescription,
distinguish exact stored discrete inputs from the continuum state, and
construct justified spectral enclosures for `c(k)` and `A(k)`. Freeze a
complete independently reviewed implementation before physical evaluation.
Keep all twelve cases and the inherited `2e-8` full-integral gate unchanged.
No continuum, ultraviolet, heating or cosmological claim follows from these
conditional finite-band results. External novelty is **NOT_ASSESSED**.
