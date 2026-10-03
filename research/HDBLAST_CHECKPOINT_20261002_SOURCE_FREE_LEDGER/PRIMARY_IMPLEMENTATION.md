# Prospective primary saved-ledger implementation

This source decodes only the separately registered immutable saved-data capsules.
It never evaluates a compact source, solves forced modes, changes a stress or
alters an earlier failed experiment. The previous metric calibration remains
FAIL. The twelve cases use both saved sources and both resolutions with
K=64,128,256, on the single inclusive interval [-2.5,-1.5]. Every retained
density/pressure sample is compared:129 coarse or257 fine.

The registered entry point is `code/diagnostic_primary.py`, with arguments
`--checkpoint-root CPP --registration-sha256 SHA --freeze-commit COMMIT
--output-dir NEW`. The shared `ledger_integrity.verify_frozen` verifies the
registration bytes, public freeze receipt and every registered input/source
before any physical array is decoded, and repeats verification afterward.
Before importing that verifier, the primary entry itself checks the explicit
FULL_REGISTRATION byte hash and pins both its own registered source bytes and
`ledger_integrity.py` against that authenticated registration. The helper is
loaded from its exact pinned path; symlink ancestry is rejected. Registered
execution requires Python3.12.14,NumPy2.2.6 andmpmath1.3.0 as well as the
binary80 ABI and represented constants below.
Output is `diagnostic.json`, with `progress.jsonl` after each arithmetic pass.
A separate `--inputs DIR --output NEW --preflight-only` interface checks hashes,
CRCs and NPY headers without decoding array payloads. Unauthenticated physical
array execution is unavailable through the CLI.

## Fixed arithmetic and exact represented inputs

Each complete twelve-case primary execution includes BOTH decimal precision
levels80 and100 under one shared900-second/262144-KiB budget. A separately
executed independent complete replay has its own clearly labeled budget; no
combined runtime claim is made. Resource receipts include the route's process
wall/RSS and the outer replay's authoritative wall/RSS.

The MP context uses exactly80/100 DECIMAL significant digits, corresponding to
269/336 binary MP bits, with zero guard digits. Transcendental seed phases are
computed at that same precision. The integer recurrence has absolute fractional
scales Q=2^266/2^333 respectively, the unique exponents p with
2^(p-1)<10^dps<=2^p. Every component/product quantization uses nearest rounding
with ties to even, including negative values. Seed phase at the anchor is1.
One step rotation is computed for each actual momentum node per precision;
there is NO phase refresh across the129/257 samples.

Native finite NumPy long-double inputs are required before conversion. Scalars
already reduced to binary64 are rejected. `as_integer_ratio` retains each exact
represented binary80 value, including the producer's numerical epsilon and pi.
An odd-mantissa/`ldexp` roundtrip handles both extreme and subnormal exponents
without overflowing an integer denominator. Zero signs are checked during
conversion and negative-zero counts are retained as input metadata. Complex
inputs are split into their native long-double real/imaginary components.
Long-double16-byte storage and63 fraction bits are required. The analytic
background is L=-1/eta from the exact stored time ratio in each MP context.

The original NumPy long-double Simpson expression is applied to the full(N,3)
stored integrand. S_b and S_a are computed first and subtracted in long double;
only then is S_ab converted exactly. Direct-reset Simpson is separately
reported. Fine doubled-step Simpson uses GLOBAL F[::2] and the matching global
prefix indices before long-double subtraction. DeltaR instead subtracts the
two exact-ratio converted stored endpoint densities in MP arithmetic.

## Free mode, profile and canonical antiderivative

For the unrestricted retained anchor mode define Omega=2k, d=w_a/(i Omega),
c=u_a-d. Both parts of c survive. No normalization or Wronskian projection is
performed. The exact free flow is u=c+dE,w=w_a E with
E=exp(i Omega(eta-eta_a)). The independent primitive follows the expanded Ward
integrand and the cancellation of its J1 moment, as derived in the frozen
proposal. The CANONICAL I_ab is its direct per-mode phase formula, evaluated
with MP exp(i Omega(eta_b-eta_a)); it uses no measured final density.

Let mu=weight*k^2/(2*pi^2), alpha=mu*d/(2*k*epsilon),
Cq=sum mu*Re(c)/(2*k*epsilon), CE=sum mu*k*Re(c)/epsilon. Define

```
A=sum Re(alpha E), B=sum Re(i Omega alpha E),
U=sum (4k^2/3)Re(alpha E).
R=CE+3L^2(Cq+A)-L B
P=CE/3-L^2(Cq+A)-U-L B
F=6L^3(Cq+A)+3L U+2L^2 B.
```

The primary kernel accumulates these three moments with exact integers in
three disjoint momentum bins, then forms the original nested K prefixes. At
each sample it projects B and U using rounded fixed-point products and rotates
the complex alpha state using the old real/imaginary components simultaneously.
The pressure U moment and unrestricted Cq/CE are retained explicitly.

The recurrence-derived primitive is reported separately. Its absolute gap
against CANONICAL direct I_ab must be<=1e-12 in every case at EACH precision.
This is an arithmetic consistency gate, not an80/100 difference of a possibly
shared bias. The signed density flow projection uses actual endpoint defects
du=u_b-u_a-expm1(i Omega Deltaeta)*d,dw=w_b-E_b*w_a. Its triangle bound applies
only to this measured projection, not to inherited-state or physical error.

## Conditional integer-rounding envelope

Treat the work-context seed alpha and rotation r as exact dyadics for this
auxiliary analysis. The bound does NOT enclose transcendental/context rounding,
the momentum rule, mode integration, inherited state or total physical error.
For Q=2^p, complex component quantization errors are<gamma=2/Q. Let M be one
plus the integer square root of the squared quantized seed components, and T
the corresponding rotation integer norm. Then |alpha|<=(M+2)/Q=A and both
work-context and quantized rotation norms are<=rho=max(1,(T+2)/Q).

If e_j is the complex recurrence error, rounded products give
e_(j+1)<=rho*e_j+gamma*A*rho^j+gamma, with e_0<=gamma. For n<=256,

```
e_j <= gamma*[1+n*(A+1)]*rho^n
rho^n <= 1/(1-n*(rho-1)), provided n*(rho-1)<1.
```

The geometric majorant follows from binomial coefficients bounded by n^m.
The implemented defining expression uses exact integer/Fraction arithmetic.
With z_bound=A*growth+e_bound, projection errors satisfy
eb<=|Omega|*e_bound+gamma*(z_bound+1) and
eu<=|4k^2/3|*e_bound+gamma*(z_bound+1). Across the interval,
|deltaR|+|deltaP|<=4L_b^2*e_bound+2|L_b|*eb+eu per mode. Positive bounds sum
over actual prefixes. Converting the rational envelope to the MP report and
each positive accumulation is inflated by16 context epsilon to retain an
upper direction through finite context rounding. This reported conditional
envelope remains separate from empirical80/100 and direct-phase diagnostics.

## Fixed outputs, gates and interpretation

Each of four records has source,setting and levels80/100 with three K rows.
Every row reports full times,R,P,F,stored_R,stored_P,stored_F; signed
S_ab,DeltaR,I_ab,D_S,D_cont,E_Q,E_flow; triangle_bound; direct-reset and
same-fine-history doubled-step Simpson; profile maxima; decomposition and
flow/operator closure; direct/recurrence primitive diagnostics and the
conditional integer envelope. Coarse rows omit the fine-only doubled fields.

EVERY non-K scalar/profile value under `records[].levels` participates in the
80/100 precision gap and decimal serialization checks, including stored
profiles, auxiliary primitive gaps and bounds. MP dyadics and serialized
decimals are compared as exact Fractions. All scientific threshold comparisons
use these exact rational/decimal values; Python binary64 is used only for
process wall-clock reporting. No observed result changes a threshold.

The100-digit values classify the fixed2e-7 R/P profile,D_cont and signed E_flow
checks. A failed scientific consistency criterion is a complete negative
observation, with classification CONSISTENCY_FAILURE and exit0. Arithmetic
closure, per-precision primitive gap, serialization, precision gap, input,
schema, authentication or resource failures are fatal execution failures.
All consistency gates must pass before attributing a reset |D_S|>2e-6 with
|D_cont|<=0.1|D_S| and |E_Q|>=0.9|D_S|. The other classifications are
LEDGER_ERROR_DEMONSTRATED and NO_GATE_SCALE_ATTRIBUTION. Every outcome retains
`old_metric_status="FAIL"`; none grants a physical calibration PASS.

The synthetic preflight fabricates complex modes with nonzero unrestricted
Re(c), source-free histories and the exact comparable shapes. It opens no
physical input files. Pure guards cover ties-even signs, binary80 extremes,
subnormals/zero signs, a binary64-loss witness, direct-MP high-frequency phase
checks, distinct global-prefix/reset arithmetic and malformed array schemas.
