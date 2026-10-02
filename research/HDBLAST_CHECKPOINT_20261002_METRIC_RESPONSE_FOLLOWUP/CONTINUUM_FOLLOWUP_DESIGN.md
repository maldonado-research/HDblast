# Separately registered continuum integration followup

The original public freeze is
`57668b8fadd75df8738565e0bbd1eb852c1ebae8`, with full registration SHA256
`f8d6bbd17b540b46c5e4382ab21fd74e15f0b283956991a250f49994553a476d`.
Its full replay stopped at the first inside-support continuum integral,
positive_B at eta=-4.5, because the registered SciPy weighted-log routine
reported roundoff and the registered warning policy correctly treated that as
fatal. The original source, partial result and failure remain immutable.

This followup changes the continuum integration algorithm only. It retains
the same physical source, metric forcing, incoming state, observations,
cutoffs, epsilon, scalar-current convention, 900-second producer budget and
scientific gates. All eleven inherited arithmetic, source, finite-K and direct
stress function bodies remain byte-identical, as do source_jet.py and the
metric contact coefficients. The independent producer and validator are
unchanged.

At an inside-support observation define ell=eta+5 and endpoint f=g_n(eta).
The exact logarithmic integral is

```
integral_-5^eta g_n(t) log(eta-t) dt
 = ell*f*(log(ell)-1)
   + 2*ell*integral_0^1 z*(g_n(eta-ell*z^2)-f)
                      *(log(ell)+2*log(z)) dz.
```

The endpoint-subtracted integrand behaves as z^3 log(z), has continuous value
zero at z=0, and has two continuous endpoint derivatives. The source is flat
at t=-5. For post-support observations, the equivalent t=-5+2z affine
integral is smooth, with the flat compact-source endpoints explicitly zero.
The eta=-3 endpoint is also defined by this exact flat limit.

Two fixed mpmath calculations run independently at 50 and 70 decimal digits,
using tanh-sinh, maxdegree=10 and the exact panels [0, 1/2, 1]. Every bump
derivative, canonical g jet, endpoint value, Euler constant, mass/log scale,
prefactor and q derivative contact is evaluated within its own mp context.
Wrapping the original binary64 source routine is prohibited. Integer source
polynomials are reused exactly; Horner, exp and rational divisions are mp
operations. The reported qjet is the 70-digit result converted once to float.

The empirical allowance for each qjet component is the maximum of the actual
50/70 precision difference, both normalized quadrature estimates and the
final float conversion error. True raw mp.quad estimates and all complete
precision qjets/endpoint jets/histories remain recorded as decimal strings.
All raw mp values are checked finite before formatting. To preserve the
unchanged validator's float receipt arithmetic, the effective raw-memory
equivalent allowance is serialized upward and normalized with the exact
binary64 expression `raw_error/(8*math.pi**2)`. This is explicitly labeled
as an empirical equivalent allowance rather than a literal mp.quad estimate.
The precision diagnostics are source-pinned implementation evidence; the
unchanged scientific validator does not separately audit these extra fields.

The original finite_q_jet function body remains untouched. Only its caller's
output-error serialization uses the validator's division expression before
the unchanged direct-response error propagation. This corrects a possible
one-ulp difference between multiplying by PREF and dividing by 8*pi^2. It
changes no source, physical value, finite integral, panel, tolerance, warning
rule or scientific gate. Fabricated error vectors test this consistency.

The exact prospective `primary.continuum` contract is exported in
continuum_metric_mp.py and copied verbatim into the new experiment. The
producer requires its hash-pinned helper and contract before any physical
call. Provenance adds continuum_algorithm_sha256 and continuum_algorithm;
each continuum block adds continuum_precision_diagnostics and explicit
log_history_error_scope. Main route, nine quantities and legacy error-map
fields are retained. The replay remains 32 commands; the existing 45-check
primary preflight also includes a separately reported 22-check algorithm
audit, without launching additional physical calculations.

Preparation evidence consists of exact symbolic identities, AST arithmetic
checks, fabricated receipt vectors and blocked authorization calls. No new
source sample, continuum or finite-K response quadrature, or mode calculation
was evaluated before the followup's separate public freeze. Precision gaps and
quadrature estimates remain empirical diagnostics, not certified total errors.
