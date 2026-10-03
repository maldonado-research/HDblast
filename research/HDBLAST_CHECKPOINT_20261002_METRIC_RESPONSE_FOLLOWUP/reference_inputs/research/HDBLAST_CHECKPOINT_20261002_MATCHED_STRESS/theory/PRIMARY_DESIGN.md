# Prospective primary matched stress implementation

This document describes an implementation prepared before the new public
registration. The inherited source is commit
`e784b18128012825f126c5bd96a4bf2033ab63d9`. No physical source samples,
response quadratures, mode integrations, stress values or trial cutoffs were
evaluated during implementation. Only exact symbolic algebra and AST checks
were executed. The shared `EXPERIMENT.json` selects all numerical inputs and
gates; the producer reads it directly. `PRIMARY_SETTINGS.json` is an unused
early design draft and need not enter the final checkpoint.

The model is the inherited minimal scalar, H=1, fixed r=2, a=L=-1/eta,
unchanged incoming BD state. The source is s=epsilon*f, epsilon=1e-4,
with f=B(eta+4) or (eta+4)B(eta+4), support (-5,-3). Exactly six observations
and cutoffs K=64,128,256 are reused by the new explicit registration. This
reuse is not authorization for a run before the new public freeze.

## Analytic source derivatives

`source_jet.py` implements the normalized conformal-time jet f0 through f5.
For d=1-u²,

    B^(n)=B P_n(u)/d^(2n),
    P_0=1,
    P_(n+1)=d² P_n' + 4n*u*d*P_n - 2u*P_n.

The integer polynomial coefficients are frozen in the module. Signed jets
use (uB)^(n)=uB^(n)+nB^(n-1). Both sources and every jet are exactly zero
outside support. The bump is computed as exp(-u²/((1-u)(1+u))) to reduce
subtraction at u=0 and at the support endpoints. Horner evaluation and an
underflow guard avoid zero times large reciprocals. Floating-point source
jets remain numerical evaluations; the separately supplied interval derivative
norms do not convert those point evaluations into certified values.

All six derivatives of both sources, and all five recurrence transitions,
are independently checked by symbolic differentiation in
`verify_primary_algebra.py`. This verifier imports polynomial constants only;
it never calls the physical source function.

## Continuum variance and its derivatives

Use q=a²deltaQ and z_n=q^(n)/epsilon. The registered logarithmic functional
is evaluated with normalized sources:

    F[f_n]=integral f_n(t)[ln(M(eta)(eta-t))+gamma_E+1]dt,
    z0=-F[f1]/(8pi²),
    z1=-(F[f2]+L*f0)/(8pi²),
    z2=-(F[f3]+2L*f1+L²*f0)/(8pi²).

The moving reference-scale terms L*f0, 2L*f1 and L²*f0 are retained. For
observations inside support, QUADPACK `alg-logb` integrates the logarithmic
endpoint; after support, ordinary quadrature integrates a smooth logarithm.
The constant ln(M)+gamma_E+1 multiplies f_(n-1)(eta) exactly, since every
initial derivative vanishes. No raw Leibniz differentiation of a singular
endpoint is used.

The physical variance derivatives divided by epsilon are

    deltaQ/epsilon       = z0/a²,
    deltaQ'/epsilon      = (z1-2L*z0)/a²,
    deltaQ''/epsilon     = (z2-4L*z1+2L²*z0)/a².

## Independent finite-band expressions

Let v=K/sqrt(K²+M²), A=asinh(K/M)-v, and

    H_K[g]=integral g(eta-tau) 2sin²(K*tau)/tau dtau.

The exact finite-band expressions are

    z0K=[f0*A-H_K[f0]]/(8pi²),
    z1K=[f1*A+f0*A'-H_K[f1]]/(8pi²),
    z2K=[f2*A+2f1*A'+f0*A''-H_K[f2]]/(8pi²),
    A'=-L*v³,
    A''=L²*v³*(2-3v²).

The differentiated source convolution follows by integration by parts on the
smooth, initially zero source, at fixed finite K. It is algebraically equal
to differentiating the finite sine/cosine history expressions; it does not
use a removed-cutoff response as its answer. Stable 2sin² replaces 1-cos.
Each deterministic time panel is at most pi/(2K) wide. Summation uses
math.fsum and the global absolute tolerance is divided by the panel count.

For all logarithmic and finite-band integrals, the root registration selects
epsabs=epsrel=1e-12 and QUADPACK limit=300. IntegrationWarnings and nonfinite
outputs fail the run. Reported estimates propagate linearly with absolute
coefficients into stress/current estimates. Exact local formulas are given
zero *quadrature* uncertainty; that does not assert zero floating-point error.
The distinct root roundoff allowance and cross-route checks remain necessary.

## Exact baseline and direct stress contacts

The physical finite-band baseline field square has the stable closed form

    Q0K=[v²/(2(1+v))-v³/48-v⁵/16]/(2pi²).

Its derivative uses v'=-L*v*(1-v²) and

    d/dv[...]=v(2+v)/(2(1+v)²)-v²/16-5v⁴/16.

This form follows by analytically integrating the complete order-two
subtraction at fixed r. The verifier differentiates it with respect to K and
matches the original combined integrand, then checks its lower endpoint and
the limit Q0=1/(12pi²). These are symbolic checks, not numerical baseline
quadrature. Q0K is kept physical and is not divided by epsilon.

Define physical-momentum moments at P=K/a:

    J5=v³/6,
    J7=(v³/3-v⁵/5)/4,
    J9=(v³/3-2v⁵/5+v⁷/7)/8.

For R=a⁴delta_rho/epsilon and Pstress=a⁴delta_p/epsilon, direct reduction of
the matched physical stress and inherited subtractions gives

    R=(3L²*z0-L*z1)/2 + a²*Q0K*f0/2 + C_rho,
    C_rho=(3L²*f0-L*f1)*v³/(96pi²),

    Pstress=(z2-3L*z1-3L²*z0)/6 - a²*Q0K*f0/6 + C_p,
    C_p=[70L²*f0*J9-(30L²*f0+10L*f1)*J7
         +(f2-L*f1-9L²*f0)*J5]/(48pi²).

Pressure uses the directly integrated fourth-order subtraction mismatch.
It is **not** obtained from a trace reconstruction. Density is **not**
obtained by solving the Ward identity. The coefficients multiplying Q0K are
the reduced formulas' mass contacts; the independent direct-mode stress
retains its original explicit +/-s|v|² operator variations before reduction.

In the continuum, set v=1, J5=1/6, J7=1/30, J9=1/105. Then

    C_p,infinity=(f2-3L*f1-11L²*f0)/(288pi²).

The independent finite-band local trace remainder is recorded as

    a⁴deltaA_K/epsilon=
      [(f2-12L²*f0)*J5-(30L²*f0+10L*f1)*J7+70L²*f0*J9]/(16pi²).

The equality -C_rho+3C_p to this remainder is a separate exact check.
Every component is retained in result JSON, so omission or sign mutations
can be diagnosed without a new primary physical run.

## Current and independent conservation diagnostics

For the separately frozen translation b=1, delta_phi=delta_x/r,

    a²delta_j/epsilon = z0 + Q0K*f0/4.

The second term is the quadratic mass-law contact. Q0 is used in the
continuum and Q0K at finite K. No scalar-current term replaces a stress test.

The producer also records a direct analytic derivative of its density:

    a⁴delta_rho'/epsilon=
      -3L³*z0+3L²*z1-L*z2/2
      +a²[(Q0K'/2-L*Q0K)*f0+Q0K*f1/2]+C_rho'-4L*C_rho,

with C_rho' computed by differentiating its explicit local formula. The
symbolic verifier obtains the same derivative directly from R, including the
derivative of a^-4. Ward cancellation is then checked as an additional
identity. A separately implemented mode-history ledger supplies the numerical
conservation check; the density is never defined by that ledger.

## Tail bounds, provenance, and stopping

`theory/stress_tail_bounds.py` supplies deferred directed-interval bounds for
q, q_prime, q_second, rho, p, current, anomaly, and Q0-Q0K. It uses analytic
jets and total-variation budgets through fifth derivatives. Its values are
called only after the public freeze and recorded separately from quadrature
estimates. The actual combined derivative/stress bounds, not the earlier
variance bound alone, govern removed-cutoff comparisons. The derivative-norm
and UV goal design may be checked analytically before freeze; no response
results may inform a later change.

The producer accepts `--output`, `--registration-sha256`, and
`--public-freeze-commit`. Every file in the prospective full manifest is
verified before source evaluation. Fresh outputs must be outside the frozen
checkpoint and must not exist. The producer preserves active points,
completed finite-K values, completed rows, command/versions/source hashes,
elapsed time, and any exception. The 900-second wall budget is checked before
and after every bounded quadrature. Failed gates require a bounded diagnosis
or a new prospective amendment, never silent tuning.

`PRIMARY_SYMBOLIC.json` and `PRIMARY_SYMBOLIC_OPTIMIZED.json` record thirty
pure symbolic implementation identities under normal and optimized Python.
No symbolic development failures occurred in this primary implementation
attempt. Earlier inherited symbolic failures remain separate historical
evidence in the inherited package and are not replaced by these new checks.

This calibration remains a fixed-geometry homogeneous scalar-to-stress/current
response at the exact reference. It does not derive metric kernels, evolve a
coupled shell, evaluate an actual shifted-root Hessian, or establish quantum
stability, particle yield, radiation transfer, thermalization, or heating.
