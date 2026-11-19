This source generator replaces the comoving pair band 60<ω<120 by a positive
quantum oscillator bath, while retaining the exact previous 64-node Gauss rule
on q<ω<60. No trajectory, outgoing yield or gravity history determines its
nodes or weights. The primary N=4 rule was selected before examining its new
source errors; the positive/negative near pulses and far pulse test transfer
from the quiet history.

For each of the 24 native q modes, set x=ω⁻² and use the positive measure

```math
d\mu_q(x)=\frac{1-q^2x}{8\pi^2x^2}\,dx,
\qquad 120^{-2}<x<60^{-2}.
```

The Gaussian nodes/weights `(x_j,μ_j)` define physical oscillators by
`ω_j=x_j⁻¹ᐟ²` and `g_j²=μ_j/x_j>0`. The two relevant moments are

```math
\begin{aligned}
C_{60,120}(q)&=\sum_j\frac{g_j^2}{\omega_j^2}
=\frac{120^2-60^2-q^2\ln(120^2/60^2)}{8\pi^2},\\
D_{60,120}(q)&=\sum_j\frac{g_j^2}{\omega_j^4}
=\frac{\ln(120^2/60^2)+q^2(120^{-2}-60^{-2})}{8\pi^2}.
\end{aligned}
```

Every N≥1 Gaussian rule integrates both C and D exactly: their integrands in
this measure have degree zero and one. These have dimensions energy squared
and dimensionless; x has dimension energy⁻², μ energy² and g² energy⁴. The
frequencies are fixed comoving energies in H₀ units and redshift as ω/a in
proper time. The implementation scales x to [−1,1] and applies reorthogonalized
Lanczos to a positive 512-node reference measure, avoiding Hankel moment
matrices. A separate 1024-node construction changes the response by less than
1.4e−14. Scaled polynomial moments through degree 2N−1 agree to 2.5e−14 or
better. `RULES.json` and `RULES.npz` document all nodes and tests.

N=1,2,4,8,16 source histories use this adapted rule. N=32 and 64 source
histories use independent direct Gauss quadrature in ω on [60,120], providing
an independent reference and refinement check. Their source differences are
below 1.3e−13 relative to the high-band peak. `RULES.npz` also exports adapted
N32/N64 rules for the independent kernel bounds.

Both the force response and vacuum noise come from the same oscillators:
the bath starts in its conformal vacuum with variance `1/(2ω_j)` and momentum
variance `ω_j/2`. Its vacuum correlation is
`Σ g_j² exp(−iω_j Δη)/(2ω_j)`, with `dη=ds/a`. The covariance hierarchy evolves
the parent response, daughter energy and parent-bath cross moments. This
preserves the quantum noise/response relation rather than replacing the
high-band force alone.

The primary static convention is the original positive completion
`a F² C_total χ²/2` per y², with
`C_total=Σ_low g²/ω²+Σ_high g²/ω²` evaluated from the actual bath. All high rules
have the same C and D to roundoff; changing the high order does not retune a
static or kinetic matching parameter. No kinetic counterterm is added. Bare
histories are also exported as a check. The low-only control includes only
the original low-band completion. Subtracting it from a total history isolates
the high-band response at order y².

The generator reuses the preceding fixed-trajectory scalar-clock covariance
equations. One C2 inverse of the fine quiet scalar calibrates `F(φ)` for every
grid/arm. The density includes parent, daughter, interaction and completion
energies; the pressures include zero interaction pressure and
`p_completion=−ρ_completion/3`; the scalar current includes the F′ terms.
The exported proper-time derivatives include F″. Native q weights are used
once. All source quantities are hatted physical coefficients per y²; the
native interface multiplies them once by `y² b/r_b`.

At the end of the gate, interaction and completion density/pressure/current
vanish exactly. The daughter conformal energy `a⁴ρ_daughter` subsequently stays
constant, while the parent response continues its free evolution on the
prescribed wall. Across all 80 histories, the largest nodal Ward residual is
4.1e−12, the largest interpolated residual 3.7e−6, and the integrated comoving
work error 4.3e−7. The after-gate daughter conformal-energy drift is below
8e−16. The ten low-only tables are bitwise identical to the preceding Λ60
source histories.

Quiet accuracy does not determine held-out accuracy. At the proper-time
readout s=2.1, the N4 high-band daughter error is −0.00241% on fine quiet,
+0.954% on near-minus, −0.801% on near-plus and −0.768% on far. Its largest
high-band source error over the whole history is 0.0246% on quiet and 1.162%
on near-minus. N8 reduces the largest held-out daughter error to 0.072%.
N16 reduces the largest high-band source error across all five walls below
9.5e−8 and readout daughter error below 1.8e−8. Errors here normalize by the
high-band reference, not by the much larger low-band contribution. These are
finite-history source comparisons, not uniform error guarantees for arbitrary
trajectories; kernel bounds and native gravity tests are separate checks.

`CHECKS.json` contains balance/closure checks and all comparison norms;
`COMPARISON.csv` is a flat LF-terminated table; `SOURCES.json` records per-case
switch times, endpoint, component readouts and the SHA256 of each temporary
source table. `NATIVE_MODES.json` authenticates the common native q rule.
Heavy source histories are generated in a temporary directory, not archived.

Reproduce from the preceding native directory containing FORK.json/WALLS.zip:

```sh
OPENBLAS_NUM_THREADS=1 python -IB source/source.py --native NATIVE --output TEMP/tables --cases fine_quiet fine_minus fine_far coarse_quiet fine_plus --high-orders 0 1 2 4 8 16 32 64 --jobs 4
OPENBLAS_NUM_THREADS=1 python -IB source/verify.py --native NATIVE --tables TEMP/tables --output TEMP/source
OPENBLAS_NUM_THREADS=1 python -IB source/quadrature.py --modes source/NATIVE_MODES.json --output TEMP/rules
```

The generator writes `{grid}_{arm}_{low|Nn}_{bare|static}.npz`. The native
interface is unchanged: arrays `s,rho,p,j,drho,dp,dj,rad,drad`; `SourceTable`
uses cubic Hermite interpolation of the supplied analytic slopes. `rad` is
daughter-only physical density. Additional arrays give all components,
comoving energies, background and gate derivatives. Every native Hermite
knot, field-gate event and s=2.1 is an integration/table knot; sources continue
to the archived endpoint near s=2.168. The verifier consumes these generated
tables without evolving them again. Python, NumPy and SciPy suffice.
