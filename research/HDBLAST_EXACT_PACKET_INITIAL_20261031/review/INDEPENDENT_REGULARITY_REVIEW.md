# Independent regularity audit

The exact-packet displacement is compatible with the required initial data
regularity when the quadratic continuation is applied after completing and
packing every component. No fatal analytic regularity defect was found.
This audit does not certify the new numerical interval error gates or a
finite response.

The independently read registered algebra is
`HDBLAST_FULL_ANGULAR_RESPONSE_20261026/initial/full23_initial.py`, lines
151–208. Its scalar displacement uses `fy=dy(fs)`; `fs` includes `vy^2`,
and `vy=dy(v)` with `v` using the first raw packet derivative. Thus the
scalar displacement includes raw packet derivatives through order three.
The five means are overwritten only after packing, and the same dependency
orders hold for those overrides. The label mapping was independently read
from `HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/projection/inputs/linear_operator.py`.

| Complete component slots | Highest q derivative in displacement | Highest q derivative in analytic velocity |
|---|---:|---:|
| K:0,2,4 | identically zero | 1 |
| L:0,2,4 | identically zero | 3 |
| N:0,2,4 | 0 | 1 |
| J:0,2,4 | 2 | 1 |
| w:0,2,4 | 3 | 2 |
| v:2,4 | identically zero | identically zero |
| u:2,4 | identically zero | 2 |
| H:2,4, T:+, T:- | identically zero | identically zero |

Displacement through rr therefore requires the raw packet through order
five. An independently reconstructed analytic velocity through r would
need the packet through order four. The retained velocity is instead the
same complete smoothed P1 vector; the new exact q does not redefine it.
`HDBLAST_HIGHER_INITIAL_JETS_20261030/checks/check_initial_second_error.py`,
lines 79–104, identifies the candidate's final-secant affine continuation
and distinguishes it from an auxiliary analytic tangent continuation used
only in an error proof. The smoothing kernel has exact triple zeros at
both support endpoints. For a finite continuous piecewise affine function,
the second distributional derivative is a finite sum of slope-jump deltas,
so its smoothed second derivative is a finite sum of translated C2 kernels.
Consequently the retained smoothed P1 velocity is at least C4 and meets the
required first-jet regularity, including its exterior.

The exact Fraction Horner interval calculation over the complete
`z=ell*r` footprint `[0,7/4000000]` separates the nominal denominators:
`rho > 129.9214` and `p < -0.0001306529`. The exact endpoint fractions are
in both canonical JSON reports. Exponential flatness dominates every
finite-order rational pole in `d(r)` at r=1. With those separated nominal
denominators, all required reconstructed jets, overridden means, and zero
slots vanish at r=1. The former cutoff locations are no mathematical joins
for the exact candidate; any different evaluation branches must enclose the
same analytic expression.

Each completed component's exterior is defined by its exact three traces,
so all 23 value/r/rr trace differences are algebraically zero. This gives
C2, generally not C3, at R=7/4. A later estimate requiring global C3 would
need additional join treatment and cannot be justified by this definition.

A raw-packet quadratic continuation is effectively rejected. At R,
`q=exp(-4/3)` and `q_rrr/q=925696/729`. The coefficient of `q_zzz` in the
scalar mean displacement is `A*vy/(6p)`, where
`A=3alpha/(rho^3 p)`. Matching only q, q_r and q_rr while setting the
exterior q_rrr to zero yields an exact nonzero scalar mean value jump
`C*q^2`, where C is a positive exact rational, approximately `9.49649e31`.
The full rational witness is recorded in the reports. Thus the control
fails continuity even before considering reconstructed second derivatives.

`check_regularity.py` imports no producer module. It authenticates every
read source before and after, authenticates its own bytes and its source
map before and after, uses explicit failures without assertions, and
requires absent output files. Normal and optimized invocations generated
byte-identical `REGULARITY_NORMAL.json` and `REGULARITY_OPTIMIZED.json`.

Source closure beyond this narrow read set, equality with the producer's
actual new expression table and velocity byte ledger, full interval
coverage and benchmark gates belong to the parent independent audit.
