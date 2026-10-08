# Exact finite-weight incoming-state transport

Candidate implementation and fabricated-only verification, 8 October 2026 UTC.
No retained array has been read or decoded by this work. No real-source callback
has been evaluated. This contribution is conditional on the independently
certified state-comparison enclosures and the later registration/GO procedure.
It does not report actual numerical incoming-state errors.

The scope is the difference between two exact linear evolutions initialized at
`a=-9/2`, one by the retained incoming state and one by the specified BD target.
Both evolutions must subsequently use the identical real forcing and complete
contacts. Their inhomogeneous difference cancels; the state difference remains.
This does not make the actual forcing, contact, saved-solver residual, time or
momentum quadrature, continuum state preparation, or ultraviolet errors zero.
In particular, these bounds do not certify later retained numerical trajectory
values without a separate propagation/residual analysis.

## Fixed exact inputs and twelve positional prefixes

Keep every retained Gaussian `dk` weight `w_j` strictly positive and every exact
represented node `k_j>0`. Define

```
Pi = 14488038916154245685 / 4611686018427387904
epsilon = 3777893186295716171 / 37778931862957161709568
mu_j = w_j k_j^2 / (2 Pi^2)
q_j = w_j / (4 Pi^2)
```

`U_saved=u_1/epsilon`, `W_saved=w_1/epsilon` use these represented constants,
without binary64 conversion. The target engine supplies exact rational
rectangles for `U_BD,W_BD`; subtraction gives rectangles for the incoming
`delta U,delta W`. Its analytic proof, separately reviewed, must establish
`c_BD=0`. Retain the exact saved invariant and its stable residual:

```
c_j = Re U_saved - Im W_saved/(2 k_j)
d_j = 2 k_j c_j = 2 k_j Re U_saved - Im W_saved
A_j = delta W_j/(2 i k_j)
k_j A_j = -i delta W_j/2
```

No Wronskian projection is performed. The first-order Wronskian defect of the
saved direction is `2 epsilon c_j`; exact nonlinear normalization includes an
additional quadratic term and is not asserted here. The constant imaginary
part of `delta U-A` is recorded in the incoming error norm even though it does
not enter the displayed linear stress operators.

Each of `positive_B` and `signed_uB` has its own coarse and fine capsule.
The ordered `K=64,128,256` prefix counts are respectively `2048,4096,8192`
on coarse and `4096,8192,16384` on fine. Thus there are four separate capsules,
49,152 distinct capsule-node occurrences and twelve prefix rows (86,016
node appearances across prefixes). Equality of source-grid arrays is never
inferred. At evaluation the API verifies strict ordering, positivity, complete
counts and exact equivalence of every positional prefix to `k<K`, excluding
`k=K`. It reports the actual weight sum and moments through order four.
It never assumes `sum w_j=K` and never regenerates or renormalizes weights.

## Signed anchor intervals preserve the invariant

Set `L=-1/eta`, so `La=2/9`, `Lb=2/7`, and let `v=Re delta W(a)`,
`z=Im delta W(a)`. The direct operators are

```
delta R = ((2 k^2+3 L^2) Re delta U-k Im delta W-L Re delta W)/(2 k)
delta P = ((2 k^2/3-L^2) Re delta U-k Im delta W-L Re delta W)/(2 k).
```

Because `Re delta U=(d+z)/(2k)`, their exact weighted anchor forms are

```
mu delta R(a) = q [(k^2+3 La^2/2)d + (3 La^2/2)z - La k v]
mu delta P(a) = q [(k^2/3-La^2/2)d - (2 k^2/3+La^2/2)z - La k v].
```

Each affine expression is evaluated on the incoming `v,z` rectangle with exact
rational endpoint ordering. The resulting interval is exactly its range on
that rectangle. Summing interval endpoints with positive retained weights
therefore encloses the signed anchor prefix difference. It retains the
canonical `c/W` correlation; independently bounding `U,W` would discard it.
The result can still overestimate the true set if the two W components or
nodes are correlated beyond the supplied rectangle product.

The API checks that the supplied U-real rectangle intersects `(d+z)/(2k)`.
This catches inconsistency but does not prove the target rectangles or its
canonical identity. Those are explicit upstream mathematical hypotheses.

## Uniform finite-sum bounds

For exact subsequent homogeneous evolution,

```
delta W(eta)=2 i k exp(2 i k(eta-a)) A,
Re delta U(eta)=c+Re(exp(2 i k(eta-a)) A).
```

Let `D=|d|`, and let `V` be the sum of the two absolute maxima of the incoming
W-component intervals. Then `|A| <= V/(2k)` and `|kA| <= V/2`. With
`p=Re(EA), r=Im(EA)`, the direct operators reduce to

```
delta R = (k+3 L^2/(2k))c + 3 L^2 p/(2k) + L r,
delta P = (k/3-L^2/(2k))c - (2k/3+L^2/(2k))p + L r.
```

Use the positive rational density triangle majorant `L+3L^2/(2k)` and pressure
majorant `2k/3+5L^2/(4k)`. Their squares exceed the exact phase coefficient norm
squares by respectively `3L^3/k` and `21L^4/(16k^2)`. Multiplication by the
actual measure before division yields the following per-node nonnegative
bounds, without reciprocal-k denominators in any weighted sum:

| Component | Uniform weighted bound, `L*=2/7` |
| --- | --- |
| density canonical | `q (k^2+3 L*^2/2) D` |
| density phase | `q (L* k+3 L*^2/2) V` |
| pressure canonical | `q (k^2/3+L*^2/2) D` |
| pressure phase | `q (2k^2/3+5 L*^2/4) V` |

These are finite sums, not numerical approximations to any proved continuum
integral. Small positive nodes require no continuum infrared assumption;
there is no represented zero node. The fused algebra also avoids forming an
LCM over the odd denominators of all `1/k_j` factors. Exact individual `c_j`
is computed only for suprema, without summing such denominators.

The implementation additionally retains signed canonical prefix sums
`C_R(L)=sum q(k^2+3L^2/2)d` and
`C_P(L)=sum q(k^2/3-L^2/2)d`. Each is affine in `L^2`, so its absolute maximum
on the time interval is exactly the larger endpoint absolute value. Therefore

```
sup_eta |delta R_prefix| <= max(|C_R(La)|,|C_R(Lb)|) + sum R_phase_bound,
sup_eta |delta P_prefix| <= max(|C_P(La)|,|C_P(Lb)|) + sum P_phase_bound.
```

The separate canonical triangle sums are also reported. The signed-canonical
improvement uses exact arithmetic and does not infer cancellation among the
unknown oscillatory terms.

## Work and time-integrated direct pressure

Direct differentiation, checked as identities over exact Laurent polynomials,
gives `delta R'=L(delta R-3delta P)` and `(-delta R/(3L))'=delta P`.
Finite weighted sums commute with differentiation/integration without any
continuum interchange argument. With `b-a=1`,

```
delta I = integral_a^b L(delta R_prefix-3delta P_prefix) d eta
        = delta R_prefix(b)-delta R_prefix(a),
delta J = integral_a^b delta P_prefix d eta
        = (b delta R_prefix(b)-a delta R_prefix(a))/3.
```

The weighted canonical work coefficient is exactly
`q [3(Lb^2-La^2)/2] d`: the leading `q k^2 d` cancels.
The pressure-time canonical coefficient is
`q [k^2/3-(Lb-La)/2] d`. Valid nonnegative components are

| Component | Weighted bound |
| --- | --- |
| work canonical | `q 3(Lb^2-La^2) D/2` |
| work phase | `q [(La+Lb) k+3(La^2+Lb^2)/2] V` |
| integrated pressure canonical | `q [k^2/3+(Lb-La)/2] D` |
| integrated pressure phase | `q [2k/3+(La+Lb)/2] V` |

For each prefix the API sums the exact signed canonical coefficients first,
takes their absolute value, and adds the corresponding positive phase sum.
The separate triangle canonical component remains in the report. These bound
the absolute value of the signed direct-pressure time integral. The integral
of the absolute pressure difference is bounded separately by the uniform
pressure bound times the exact duration one. Pressure is never defined or
repaired from a conservation residual.

## API, receipts and reproducibility

`discrete_transport.py` is standard-library-only and has no data reader,
physical source evaluator, network action, trajectory evaluator or CLI that
could accidentally open a retained capsule.

- `Interval(lo: Fraction, hi: Fraction)` stores exact closed rectangles.
- `NodeError(k,weight,d,delta_u_re,delta_u_im,delta_w_re,delta_w_im)` stores one
  incoming difference. All scalar inputs must be `Fraction`, not float.
- `from_saved_target(k,weight,u_saved,w_saved,target,exponent=-96)` accepts
  already normalized saved complex pairs and the integer endpoint mapping
  supplied by `entire_target_engine`. It computes exact saved `d` itself.
- `aggregate_capsule(capsule,nodes,expected_count,prefixes)` supports explicit
  fixture plans and verifies every actual prefix/node relation.
- `aggregate_registered(capsules)` requires the four exact labels
  `positive_B/coarse`, `positive_B/fine`, `signed_uB/coarse`, `signed_uB/fine`.
  It produces only twelve rows. `encode_exact` serializes rational numerators
  and denominators as strings, with no rounding or per-node output.

Incoming error norms include supremum and weighted L1 rectangle upper bounds
for U and W, `sup |c|`, the exact first-order Wronskian-defect supremum, and
stable `kA` bounds. U/W error norm fields refer to the incoming anchor;
all-time stress fields are explicitly named `uniform`.

Run fabricated-only validation with fresh output filenames:

```
python verify_transport.py --full-fabricated --output NORMAL.json
python -O verify_transport.py --full-fabricated --output OPTIMIZED.json
```

The retained full fabricated receipts each contain 47 checks and 17 mutation
controls; both have scientific SHA256
`bf6fdcb498ef2766c2133a60a04302a070f4e3b3180b42a1fb9dfd21ff097d52`.
Exact 49,152-node fabricated summaries agree byte-for-byte in both modes,
produce twelve rows and occupy 76,747 bytes. Concurrent benchmark runs took
11.411 and 11.760 seconds on this environment. There are zero retained-array
reads/decodes and zero real callbacks in these receipts.

The verifier uses canonical sparse Laurent polynomials over Fraction for the
universal direct, canonical, fused, continuity, primitive and square-gap
identities. Fabricated rational rectangle corners and rational unit-circle
phases independently exercise the numerical API. Controls reject nonpositive
weights/nodes, floats, reversed intervals, incompatible state projection,
missing/extra nodes or capsules, wrong cutoffs, cutoff equality, duplicate
nodes, order errors, and altered signs. The paired canonical geometry omission
is explicitly shown to pass a Ward test while differing from both direct
operators. This is an internal exact algebra check, not external peer review
or proof-assistant formalization.

The historical metric calibration remains **FAIL_UNCHANGED**. The original
full twelve-case `2e-8` pressure/contact certificate, continuous integrals,
source/contact evaluation, time quadrature and ultraviolet completion remain
**UNRESOLVED**. A higher-dimensional Big Bang cause remains
**NOT_ESTABLISHED**, and external mathematical novelty remains
**NOT_ASSESSED**. No bound here promotes the full gate to PASS.
