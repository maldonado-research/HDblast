# Registration protocol: uniform source and Duhamel operator prerequisite

**AWAITING_REMOTE_FREEZE_NO_REAL_SOURCE_EVALUATIONS.** This is a narrower prerequisite
for the proposed full twelve-case continuity certificate. It does not decode
archived modes, momentum nodes, weights, histories or any other physical array.
It does not compute action pressure, finite-band contacts or the original
ledger attribution. The full design is preserved in `full-integral-proposal/`
with its `2e-8` full-integral target explicitly unresolved.

## Exact target and fixed universe

Let `a=-9/2`, `b=-7/2`, `z=eta+4`, `L=-1/eta`. The mathematical source is
`h=exp(-z^2/(1-z^2))` for `positive_B`, or
`h=z*exp(-z^2/(1-z^2))` for `signed_uB`. Define
`g=4*L^2*h-2*L*h'-h''` from analytic differentiation of those functions.
The interval is inside the support; no flat support endpoint or larger history
is included. This real analytic source differs from treating previously
rounded native source samples as exact. No real source callback is permitted
until the new remote freeze has been verified.

Partition `[a,b]` into exactly 64 equal panels. For panel `j=0,...,63`,
`l_j=-9/2+j/64`, `r_j=l_j+1/64`, `c_j=l_j+1/128`, and half-width `h_p=1/128`.
The analytic remainder proof covers both sources, all these panels and every
real momentum `k` in `[0,256]`. It must also cover `k=0` as a removable limit.

The numerical probe universe is fixed before evaluation:

```
sources = [positive_B, signed_uB]
panels = [0,1,...,63]
k = [0, 1/2^40, 1/2^12, 1/4, 1, 16, 64, 128, 256]
```

All probes use those exact rational numbers; `k` is not read from an archive.
Each source/panel/momentum row encloses

```
M0(l,r)      = integral_l^r g(s) ds
Mexp(l,r;k)  = integral_l^r exp(2*i*k*(r-s))*g(s) ds
Mu(l,r;k)    = integral_l^r Phi_k(r-s)*g(s) ds
Phi_k(t)    = (exp(2*i*k*t)-1)/(2*i*k), with Phi_0(t)=t.
```

Also enclose these three moments over the complete `[a,b]` for every source
and momentum probe. This gives **1152 panel rows and 18 whole-interval rows**
per registered arithmetic configuration: 1170 moment rows in total. `M0` is
retained in every row so the output contract is uniform, despite its momentum
independence. No source/panel/momentum subset is selected after results.

Register the separately proved `|L*g|<=32` analytic-disk bound and the real
source-work channel
`J_Lg(l,r)=integral_l^r L(s)*g(s) ds`, with 128 panel rows and two whole-interval
rows per configuration. This channel is included in both routes. That scalar
is a source primitive, not the
full action-derived finite-momentum baseline work or contact integral.

The retained represented constants are

```
epsilon_rep = 3777893186295716171 / 37778931862957161709568
pi_rep = 14488038916154245685 / 4611686018427387904.
```

These constants remain exact represented rationals if the response operator is
reported. The unscaled moments above do not require momentum weights or pi.
The response operator for arbitrary incoming modes is the linear map

```
w(r) = exp(2*i*k*(r-l))*w(l) - epsilon_rep*Mexp
u(r) = u(l) + Phi_k(r-l)*w(l) - epsilon_rep*Mu.
```

No actual incoming state is supplied here. A bound on this operator's forcing
term is not a bound on any archived state or finite-momentum stress integral.

## Separate source, phase and arithmetic guarantees

The primary has fixed source Taylor degree 24, 64 panels, analytic
disk radius `R=1/8`, and a proved uniform Cauchy bound `M_g=64`. For
`h_p/R=1/16`, the exact source-model remainder is

```
R_g = 64*(1/16)^25/(1-1/16) = 1/(15*2^90).
```

The source-work bound `M_Lg=32` gives
`R_Lg=1/(15*2^91)`. These are uniform analytic remainder bounds, contingent on
the published source/complex-domain proof. The polynomial coefficients still
need outward arithmetic enclosures. A small probe ball alone does not prove
the Cauchy bound or its derivatives.

For local centered coordinate `x` in `[-1,1]`, set `q=2*i*k*h_p`; `|q|<=4`
and `|q*(1-x)|<=8`. Use degree-96 entire series for the phase and regularized
kernel, including zero momentum without division. A conservative exact phase
tail bound is `4096*8^97/97!`; a centered regularized kernel bound is
`2*4096*8^97/98!`. The proof must define the exact series index convention and
show these tails belong to the actual implementation. A separately proved
real-phase Taylor tail may be tighter, but changing a registered tail or degree
after source evaluation is prohibited.

These are **pointwise** bounds. Integrating a centered kernel on `[-1,1]`
adds a factor of two: its kernel-moment tails are
`2*4096*8^97/97!` and `4*4096*8^97/98!`, respectively. Subsequent physical
panel-moment formulas also apply their `h_p` or `h_p^2` factors and the proved
absolute source-polynomial coefficient bounds. Those conversions must appear
in the output error inventory, rather than treating a pointwise bound as a
completed integral radius.

The complete integral certificate includes source-model tail, phase/kernel
tail, coefficient balls, integrated arithmetic, global transport arithmetic
and exact output conversion. Keep all those terms separately in the output.
Uniform analytic model error and the total radius of an evaluated Arb ball
are different quantities. A uniform proof applies to the declared analytic
domain; numerical ball-width claims apply only to the 1170 evaluated rows.

Whole-interval moments must use valid composition of locally enclosed moments,
not a degree-96 Taylor expansion at a whole-interval phase as large as 512.
For distance `d_j=b-r_j`, a valid additive identity is

```
Mexp(a,b;k) = sum_j exp(2*i*k*d_j)*Mexp(l_j,r_j;k)
Mu(a,b;k)   = sum_j [Phi_k(d_j)*M0(l_j,r_j)
                    +exp(2*i*k*d_j)*Mu(l_j,r_j;k)]
M0(a,b)     = sum_j M0(l_j,r_j).
```

This preserves the exact target without repeated rectangular interval wrapping.
Global rotations are evaluated by a pinned validated transcendental routine.
Evaluate `Phi_k(d)` with its entire series when the exact rational phase
magnitude is at most one, including `k=0`; at larger phase it may use validated
`(exp(2*i*k*d)-1)/(2*i*k)` with a proved nonzero denominator. Record the
additional entire-series tail and rounding in the complete radius. Never
recover small-momentum `Mu` by subtracting `Mexp-M0` and dividing by `2*i*k`.

Error propagation is absolute and nonnegative. On the real domain,
`|exp(2*i*k*t)|=1` and `|Phi_k(t)|<=t`. Thus a uniform source error `R_g`
contributes at most `(r-l)*R_g` to `M0/Mexp` and `(r-l)^2*R_g/2` to `Mu`.
This statement is independent of the sign of the source. Carry kernel and
arithmetic contributions by explicit proved triangle bounds as well; agreement
of signed residuals is not an absolute error bound.

## Exact output and acceptance

The new narrow gate is **whole-interval total absolute enclosure radius at
most `1e-26`** for `M0`, `Mexp`, `Mu`, and `J_Lg`. The gate includes
every remainder and arithmetic/output contribution. It is not the unresolved
full continuity gate `2e-8` and does not establish it. All panel enclosures are
also retained and checked; final whole-interval values cannot conceal a panel
failure or omitted row. The primary configuration is Arb256/source24/phase96,
with source Taylor jet cap27. The independent configuration is dyadic512/
source24/forced ODE96. There is no post-result precision, polynomial-degree or
probe selection. Complete fabricated implementation tests must exercise these
same configurations before freeze.

Emit each real component as an exact reduced rational interval
`{"lo":"n/d","hi":"n/d"}`, with positive denominator and `lo<=hi`.
Extract actual Arb lower/upper endpoints or the independent method's exact
rational enclosure endpoints; general reduced rational denominators are
allowed. A printed ball midpoint is not an enclosure. For a complex rectangular
output define its conservative
absolute radius as `radius(Re)+radius(Im)`, where each component radius is half
its interval width. This exact rational upper bound avoids an uncertified
square root and cannot understate the Euclidean radius. The `1e-26` gate applies
to that complete bound. Imaginary `M0` and `J_Lg` are exactly zero and may be
represented as `[0,0]`.

The registered real-target projection makes the imaginary parts of `M0` and
`J_Lg` exactly zero, and the imaginary parts of `Mexp/Mu` exactly zero at `k=0`.
The proof establishes these exact analytic identities; the real enclosure is
retained. This removes generic complex-ball imaginary padding for a known real
quantity and does not project or alter an incoming physical mode.

Return `PASS_UNIFORM_MODEL_AND_REGISTERED_PROBE_CERTIFICATE` only after the
uniform proofs, exact row membership, full target enclosures, every registered
control, outward serialization, no unauthorized input decode and resource
receipts pass. A mathematically valid but wide ball is
`UNRESOLVED_PROBE_RADIUS`. A failed remainder/identity/containment comparison is
`CERTIFICATE_CONSISTENCY_FAILURE`. Missing pins, marker/schema failures,
unsupported arithmetic, timeout/RSS overrun or failed serialization are
`FAIL_EXECUTION`, preserving partial outputs. No tolerance or probe removal is
allowed after evaluation. No status in this checkpoint may be named a
completed twelve-case continuity certificate or full ledger repair.

## Chronology, fabricated controls and remote freeze

The earlier empirical active-source diagnostic, its original failures and its
registration-marker repair are known. This is a disclosed follow-up, not blind
confirmation. A single versioned stdlib contract is shared by both entry
wrappers and the outer supervisor. Its canonical marker
`chronology.new_target_evaluations_before_public_verification` requires exact
integer zero counts for every registered real-source route. Booleans, strings,
renamed keys and generic "all project evaluations are zero" claims are rejected.

Before freeze, run the exact numerical implementation on fabricated polynomial
sources and rational grids with the same 64-panel/nine-momentum/two-source
output dimensions and every registered control. Required controls cover zero
and tiny momentum, high phase, analytic polynomial moment containment,
source bias, reversed phase sign, reversed forcing sign, omitted or understated
source/phase remainder, incorrect whole-panel transport, inward-rounded output,
duplicate/omitted rows, marker rename, and complete entry-to-validator replay.
The marker-rename control must fail before the fabricated decode/source stub;
both actual production wrappers must exercise the exact finalized registration
schema so the previous execution-contract defect cannot recur.

Measure two **complete fabricated whole-route** runs for each real-source route,
including source coefficient generation, all panel/global moment channels,
every configuration, JSON serialization and final validation. Use non-root,
one numerical thread and no child processes. The external supervisor's `wait4`
wall/RSS record is authoritative. Each route has a fixed budget of **120 seconds
and 131072 KiB (128 MiB)**, selected before any real source evaluation. Both
complete fabricated implementations, including source-coefficient construction,
must demonstrate feasibility before remote freeze. If that fabricated evidence
requires a different limit, it must be chosen, documented and reviewed before
freeze. Do not inherit the old 900-second/262144-KiB limits or enlarge a frozen
budget after a real-source timeout. A prototype-only gate test is not a
whole-route resource benchmark.

After independent review, publish the exact source, schema, proof, parameters,
fixture outputs, registration and complete immutable inventory. Independently
download all registered public blobs at the full frozen commit and verify
SHA256/length/membership before any real source callback. Authenticate those
pins with stdlib before numerical-core import and use fresh external outputs.
Check registered bytes before, between and after route execution. Preserve all
attempts and actual chronology if a repair becomes necessary.

The inherited opaque input manifest may be copied as context/provenance, but
**its capsules are not execution inputs** to this narrow checkpoint. Its old
zero-evaluation notes belong to the old preparation and do not establish a new
prospective zero. No archived physical array is needed or permitted to be
decoded by the new experiment.

## Unresolved larger claims

The full twelve-case saved-discrete momentum continuity integral, direct
action pressure, finite-cutoff work/contact integration, `2e-8` integral
half-width, original metric failure repair and all physical conclusions remain
unresolved. The original 29 endpoint and 30 refinement failures remain FAIL.
This prerequisite can validate a numerical source/forcing component. It does
not certify an incoming state, a whole-history solution, continuum momentum,
coupled gravity, heating, a hot Big Bang or a higher-dimensional cause. No
external mathematical novelty claim is made.
