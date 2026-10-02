# Exact derivative norms for the prospective stress tail

Analytic staging, 2 October 2026. No pulse/source values, response integrals,
mode solutions, or floating interval primitive values have been evaluated.
Only exact polynomial arithmetic, rational critical-root isolation, and the
optional rational analytic global-envelope certificate were executed.
The callable interval library is for use after the encompassing stress
registration has been frozen.

Use the inherited profiles `B(u)=exp(1-1/(1-u^2))` and `F(u)=uB(u)` for
`|u|<1`, extended smoothly by zero outside. Here `u=eta+4`, unit half-width,
unit amplitude, and all endpoint derivatives vanish. For j=0,1,2 define

    N_j(u)=|f^(j+2)(u)|+integral_-1^u |f^(j+3)(v)|dv.

The integral is zero for u<=-1 and saturates at u=1. These are the exact
normalized derivative norms required for the first three sine-transform
integration-by-parts remainder bounds. They are not response quadratures.

## Derivatives through fifth order

Write `f^(n)=P_n(u)B(u)/(1-u^2)^(2n)`. The exact recurrence is

    P_(n+1)=(1-u^2)^2 P_n'+[4n u(1-u^2)-2u]P_n,
    P_0=1 for B, P_0=u for F.

The needed higher derivatives are

    B'''=-4u(6u^6+3u^4-10u^2+3)B/(1-u^2)^6,
    B''''=4(30u^10+45u^8-132u^6+58u^4+6u^2-3)B/(1-u^2)^8,
    B'''''=-8u(90u^12+270u^10-765u^8+300u^6
                  +264u^4-170u^2+15)B/(1-u^2)^10,

    F'''=-2(3u^8+24u^6-26u^4+3)B/(1-u^2)^6,
    F''''=4u(6u^10+81u^8-92u^6-46u^4+70u^2-15)B/(1-u^2)^8,
    F'''''=-4(30u^14+615u^12-570u^10-1235u^8
                   +1738u^6-555u^4-30u^2+15)B/(1-u^2)^10.

These are derivatives with respect to u. The module also supplies exact
polynomials for orders zero through two. For a general half-width sigma
and amplitude epsilon, the physical norm is
`|epsilon| sigma^(-(j+2)) N_j`; the endpoint source derivative of order n
has factor `epsilon sigma^(-n)`. Amplitude is not silently inserted twice.

## Critical points and exact total variation

For m=j+2, the integral in N_j is the total variation of `g=f^m`.
The critical points of g are the zeros of `f^(m+1)`. Remove its nonzero
constant and optional factor u, and set z=u^2. The resulting polynomials
and root brackets are below. Each displayed integer n denotes the exact
interval `(n/256,(n+1)/256)` in z.

| Source | m | Critical polynomial in z | Bracket numerators | Additional root u=0 |
| --- | --- | --- | --- | --- |
| B | 2 | `6z^3+3z^2-10z+3` | 95, 205 | yes |
| B | 3 | `30z^5+45z^4-132z^3+58z^2+6z-3` | 59, 184, 224 | no |
| B | 4 | `90z^6+270z^5-765z^4+300z^3+264z^2-170z+15` | 27, 166, 213, 233 | yes |
| F | 2 | `3z^4+24z^3-26z^2+3` | 116, 207 | no |
| F | 3 | `6z^5+81z^4-92z^3-46z^2+70z-15` | 76, 188, 224 | yes |
| F | 4 | `30z^7+615z^6-570z^5-1235z^4+1738z^3-555z^2-30z+15` | 47, 170, 214, 233 | no |

Exact Sturm counts certify respectively 2,3,4 roots in (0,1) for each
profile's three rows. Every bracket has exactly one sign-changing simple
root; polynomial gcd checks exclude repeated roots. Parity supplies the
negative-u roots, and the explicitly listed zero roots complete the list.
There are no other critical points inside the pulse. No floating root or
source evaluation was needed for these statements.

List the critical points below the observation u in increasing order as
`c_1,...,c_l`, and set `c_0=-1`, `c_(l+1)=min(u,1)` with `g(c_0)=0`.
The exact variation and norm are

    integral_-1^u |g'| = sum_(i=0)^l |g(c_(i+1))-g(c_i)|,
    N_j=|g(u)|+sum_(i=0)^l |g(c_(i+1))-g(c_i)|.

At observations after the pulse, g(u)=0 and the final value is also zero.
At observations before the pulse, N_j=0 exactly. A critical point coincident
with the observation may be included; it contributes a zero final increment.
This construction uses actual extrema rather than interval-length times a
large pointwise derivative envelope.

## Deferred interval implementation

`pulse_derivative_budgets.py` has no source or interval evaluation at import.
Its primary API is

    pulse_interval_data(iv, source, eta)
      -> {'jets': [f,f',...,f^(5)],
          'budgets': [N_0,N_1,N_2],
          'details': serializable metadata}.

The caller supplies an `MPIntervalContext` with exactly 60 decimal digits.
Source names are `positive_B` and `signed_uB`; the returned quantities have
unit amplitude and unit half-width. The caller can retain the same interval
context through the complete stress-tail arithmetic, including geometry,
reference mass, pi, and normalization `a^4 delta rho/epsilon` or pressure.

Each algebraic root is enclosed by 256 exact rational bisections of its
declared bracket. Its rational endpoints are outward rounded before taking
the interval square root. Interval evaluation of the derivative polynomial,
exponential, and positive denominator encloses every extremal primitive
value. Summing absolute interval increments encloses the exact variation.
Root ordering relative to the observation uses rational comparisons of u^2;
if the fixed root precision cannot settle an ordering and the observation
is not an exact root, the module raises an explicit error.

The six inherited observations are all exactly ordered by the original
coarse brackets, before any bisection or source evaluation. The convenience
wrapper `normalized_derivative_budgets` converts bounds to binary64 and
advances each nonzero upper endpoint toward positive infinity. It preserves
exact zero before the pulse. Parent tail arithmetic should use the supplied
interval API until its final upper-endpoint conversion.

## Exact proof and scope

`verify_pulse_derivative_budgets.py --output NEW_JSON` verifies the derivative
recurrence, parity, simple roots, exhaustive Sturm counts, each rational
bracket, and exact ordering at inherited observations. It imports only the
module's declarations and never calls the interval/source evaluator or its
root-bisection routine. Explicit exceptions keep all checks active under
Python optimization. The output pins the actual module and proof source
hashes. Normal and optimized actual logs and JSON reports accompany it.

## Optional global integer envelopes from rational inequalities

For `g=f^(j+2)`, `N_j'=sign(g)g'+|g'|>=0` wherever differentiated, with
the same monotonic conclusion through zeros by absolute continuity.
Consequently the global maximum of N_j is exactly the full total variation
of g, because g vanishes after the pulse. No additional supremum of |g| is
needed. The following conservative integer envelopes hold at every time:

| Profile | N0 upper | N1 upper | N2 upper |
| --- | --- | --- | --- |
| B | 97 | 2926 | 161866 |
| uB | 96 | 2760 | 154078 |

`prove_global_rational_envelopes.py` certifies these bounds using only exact
rational operations. It refines each already certified root bracket by
40 rational bisections, encloses its square root by 40 rational bisections,
and bounds the derivative polynomial and denominator by rational interval
arithmetic. For `t=z/(1-z)>0`, it bounds `B=exp(-t)` using the reciprocal
of positive-exponential Taylor inequalities. With N=80 and t<N+2,

    S_N(t) <= exp(t)
      <= S_N(t)+[t^(N+1)/(N+1)!]/[1-t/(N+2)].

The lower and upper rational t endpoints give outward reciprocal bounds for
B. The ratio bound on all subsequent positive Taylor terms proves the upper
remainder. Rational extrema intervals give an upper bound for the exact sum
of absolute increments, and an exact integer inequality certifies each final
ceiling. No floating exponential, observed source sample, response integral,
or interval-module source evaluation is used. These are analytic global
envelopes for prospective hard-budget feasibility; the registered runtime
module retains its separately frozen 256 root bisections and 60-digit
interval arithmetic for tighter time-specific norms.

`GLOBAL_RATIONAL_ENVELOPES.json` and its optimized counterpart preserve the
actual analytic certificate outputs and hashes. This optional certificate
does not change the pulse library or any source parameter.

The symbolic and rational proofs do not execute or validate mpmath interval runtime.
The deferred implementation requires independent static review and later
registered runtime checks. These norms bound analytic ultraviolet remainder
terms; they supply no mode, finite-momentum quadrature, logarithmic time
quadrature, finite precision, or cancellation error allowance.
