# Candidate verified source and Duhamel operator

Status: **PREPARATION_ONLY_NO_REAL_SOURCE_EVALUATIONS**. This document describes
an executable candidate and proof, not a physical result or a prospective
registration already frozen. All preparation entry points use fabricated
sources/nodes. No checkpoint array was decoded and no real bump callback ran.

## Target and limits

The completed experiment specifies the exact interior analytic functions

`eta in [-9/2,-7/2]`, `z=eta+4`, `L=-1/eta`,
`B=exp(1-1/(1-z*z))`, `h=B` or `h=z*B`,
`g=4*L*L*h-2*L*h_prime-h_second`.

The inherited represented pi and epsilon ratios remain the exact rational
constants recorded in EXPERIMENT.json. The new moments below are unscaled
source integrals, so pi does not enter them. Multiplication by epsilon, when
used to construct a forced response, must use its registered represented ratio
3777893186295716171/37778931862957161709568, not an idealized decimal replacement.

The original two sources, two saved momentum rules, three cutoffs and twelve
ledger cases remain the later finite target. This narrower prerequisite checks
the source/operator only and does not decode or classify those cases. It
supplies no incoming-state error, saved-mode error, discrete-stress/contact
certificate, direct-pressure integral, momentum-continuum limit or Big Bang
mechanism. The old FAIL findings are unchanged.

## Settled candidate geometry and probe universe

Use exactly 64 equal panels. Panel j has center
`c_j=-9/2+(2*j+1)/128`, half-width `H=1/128`, normalized coordinate
`x=(eta-c_j)/H` in `[-1,1]`, and Taylor degree24. The analytic radius is `R=1/8`.
All geometry is rational. The moment phase degree is96 and Arb midpoint
precision is256 bits; these values must be part of the future public freeze.

Root's proposed narrow physical universe consists of both sources, all64
panels, and exact rational momenta
`{0,1/2^40,1/2^12,1/4,1,16,64,128,256}`. Each panel has M0/Mexp/Mdrift output;
each source/momentum also has the whole-interval triple. This gives1152 panel
rows and18 whole rows per registered precision. Source-work Lg integrals add128
panel outputs and two whole sums. The complete output schema, final resource
budget and nonroot external receipts are owned by the future driver.

The proposed whole-interval gate is certified complete radius at most1e-26
for each unscaled moment, computed as the sum of real and imaginary outward
endpoint half-widths for a complex moment. This is a new engineering
gate for this prerequisite, not a replacement of the inherited2e-8 full-ledger
proposal or its unchanged attribution gates. A finite interval exceeding this
gate is an unresolved certificate; no precision, mesh or degree may be selected
after seeing a favorable physical result unless a new separately disclosed
registration is made.

## Analytic proof and actual coefficient enclosures

On each complex disk `|eta-c|<=1/8`, with real c anywhere in the whole interval,
`|z|<=5/8`, `|1-z*z|>=39/64`, and `|L|<=8/27`. For `w=z*z`, the exact inequality
`Re[1/(1-w)] >=1/(1+|w|)` gives `|B|<=exp(25/89)<=89/64`. This uses a real-part
bound, not the incorrect absolute-reciprocal upper bound. Exact differentiated
formulas give `|g_positive|<64` and `|g_signed|<64`. Therefore the common Cauchy
majorant is64 and the degree24 normalized Taylor source remainder is

`R_g=64*(1/16)^25/(1-1/16)=1/(15*2^90)`.

The same disks satisfy `|L*g|<=512/27<32`, so the source-work Taylor model uses
`R_Lg=1/(30*2^90)`. The bounds are proved uniformly, not inferred from sampled
values. The full proof and exact-rational arithmetic checks are in
`cauchy-proof/CAUCHY_REMAINDER_PROOF.md` and `check_bound_constants.py`.

The real callback uses Arb formal series of h through degree26, differentiates
with respect to eta by dividing normalized derivatives by H, and extracts
g and Lg through degree24. Every returned coefficient is an Arb ball, not an
unvalidated native number. Ball uncertainty remains in all operations.

Important backend requirement: `arb_series(...,prec=27)` alone does not override
the global `ctx.cap`, which defaults to10. The callback explicitly sets cap27,
requires result precision at least25, and restores the prior cap in a finally
block. An independent fabricated audit caught this silent-truncation hazard
before physical evaluation. The repaired implementation also requires at
least256-bit arithmetic before evaluating source coefficients or kernels.

Uniform source error can be presented separately for an exact rational
midpoint polynomial: if coefficient ball j has exact midpoint m_j and radius
r_j, the exact source differs from sum(m_j*x^j) by at most
`sum(r_j)+R_g` everywhere on its panel. A tighter L1 source error is
`2H*sum(r_j/(j+1))+2H*R_g/26`. The source model remainder alone is never the
complete evaluated moment radius.

## Entire, zero-safe local phase kernels

For lambda=2ik and z=lambda*H, define

`E_j(z)=integral[-1,1] x^j*exp(z*(1-x)) dx`,
`Q_j(z)=integral[-1,1] x^j*(exp(z*(1-x))-1)/z dx`.

Q is evaluated as its entire series, including at z=0. Neither local kernel
uses a division by k or a recurrence containing1/k. If the source polynomial
is sum(a_j*x^j), the three local moments are

`M0=H*sum(a_j*J_j)`, where `J_j=2/(j+1)` for even j and0 for odd j;
`Mexp=H*sum(a_j*E_j)`;
`Mdrift=H*H*sum(a_j*Q_j)`.

Let `A_jm=integral[-1,1] x^j*(1-x)^m dx`. The implemented exact coefficient is
`A_jm=2^(m+1)*sum[ell=0..j] binomial(j,ell)*(-2)^ell/(m+ell+1)`.
E has coefficient `A_jm/m!` and Q has coefficient `A_j,m+1/(m+1)!`.
The degree96 coefficient arrays are rational balls and acb_poly evaluates their
finite polynomials in compiled FLINT. The family is cached per rational k and
uniform half-width, then reused across all panels and both sources.

For k in[0,256], `|z|<=4` and `|z*(1-x)|<=8`. The exact rational proof
`e<=49/18<11/4` gives `exp(8)<(11/4)^8<4096`. Thus after degreeT=96,

`tau_E=2*4096*8^(T+1)/(T+1)!`,
`tau_Q=4*4096*8^(T+1)/(T+2)!`.

Both real and imaginary kernel components are inflated by these nonnegative
bounds; a complex disk is safely contained in this square. Kernel coefficient
rounding and polynomial-evaluation rounding remain enclosed by acb. No fit
comparison, high-precision agreement or sampled maximum substitutes for a
remainder proof.

The centered source coordinate is retained in E/Q; there is no translation to
the left endpoint. Therefore the source coefficient l1 bound is
`64*sum((H/R)^j)`, not the larger translated `(3H/R)^j` bound.

## Moment source tails and whole-interval composition

The true local analytic source differs from its ball polynomial by at mostR_g.
Real-axis phases have modulus1 and `|Phi(d)|<=d`, where
`Phi(d)=(exp(2ikd)-1)/(2ik)` is interpreted continuously at k=0. Consequently
the added local tails are `2H*R_g` for M0 and Mexp, and `2H^2*R_g` for Mdrift.
Source-work uses its separateR_Lg. All coefficient-ball, phase, kernel and
accumulation arithmetic is additionally retained in the final Arb balls.

For a panel ending at r and a final endpoint b, its global contributions are

`A += exp(2ik*(b-r))*Mexp_local`,
`B += Phi(b-r)*M0_local+exp(2ik*(b-r))*Mdrift_local`.

This is the additive form of the exact Duhamel semigroup. It avoids repeatedly
rotating already inflated rectangular state balls. The identity
`Phi(d+s)=Phi(d)+exp(2ikd)*Phi(s)` proves the formula. The whole-interval M0 is
the sum of local M0; A and B are the whole Mexp and Mdrift.

Arb exp encloses global rotations even when the phase reaches512. There is no
degree96 global power series at phase512. Phi uses its separate entire degree96
series only when `|2k*d|<=1`, with tail `3*d/98!`; otherwise Arb exp and division
by the exact nonzero represented frequency enclose Phi. This handles k=0 and
k=1/2^40 without subtracting Mexp−M0 and dividing by tiny k.

For the anchored forced equations `u'=w`, `w'=2ik*w-epsilon*g`,
`w_b=phase*w_a-epsilon*A` and `u_b=u_a+Phi(T)*w_a-epsilon*B`. The generic response
routine also returns positive source-only error bounds. On a unit interval and
uniform source tailR_g these are `abs(epsilon)*R_g` for w and
`abs(epsilon)*R_g/2` for u. These are error contributions, not the complete ball
radius or an error on any inherited incoming state.

## Exact API and serialization

`future_registered_source_bundle(c,source,authorize)` returns a mapping with
`g` and `Lg`, each a `PolynomialPanel(center,half_width,coefficients,remainder)`.
The caller-owned authorize callback must authenticate the reviewed public
freeze, all source/schema/registration hashes, and the full input contract
before any callback evaluation. Without it the function raises. This module
does not implement or bypass the future complete freeze guard.

`EntireKernelFamily().at(k)` returns the two25-element kernel tuples;
`panel_moments(panel,kernels)` returns M0/Mexp/Mdrift and exact source-tail terms;
`forced_response(panels,k,acb(0),acb(0),Fraction(1),family)` returns the whole
unscaled moments. Real coefficients and incoming modes must be passed as
proper exact-rational-to-ball enclosures, never binary64 conversion.

`rational_endpoints(arb)` calls Arb's outward lower/upper methods and exports
their exact fmpq values as rational strings. `complex_endpoints(acb)` does this
for both components. Decimal midpoint formatting is never treated as a bound.
The consumer must compute gates from the exact interval endpoints and retain
both analytic model budgets and complete final ball radii.

## Fabricated evidence and remaining requirements

`test_fabricated_moments.py` passes in normal and optimized Python. It tests
exact zero phase, tiny phase, phase-cap256, constant-source global moments
against separately adaptive validated Arb integrals, source-coefficient
uncertainty, biased forcing, phase/forcing sign mutations, omitted phase
remainder, high degree source jets, no-freeze rejection and exact endpoint
serialization. Independent `cauchy-proof/audit_fabricated.py` compares ten
manufactured polynomial response cases against1024-bit closed polynomial
integrals and checks25 exact g and25 exact Lg coefficients.

`benchmark_narrow_fabricated.py` covers the complete proposed1152+18 row universe
and128 work rows, with nonroot uid1000 and no real callback/array decoder. Its
receipt is `FABRICATED_NARROW_RECEIPT.json`. Full fabricated-grid moment-only
benchmarks at8192 and16384 nodes are recorded separately; they do not establish
pressure/contact or twelve-case driver feasibility. The executable source was
audited at SHA256 de995d2edeca3944400fd90b892521f05c941b6d8a1f75c2a820740b82d995c0.

Before a physical run, root must settle and validate a complete driver/schema,
exact public freeze/readback, nonroot full-entry normal/optimized tests and
outer resource receipts, versioned backend/platform guarantees, all finite
probe membership and failure semantics, and a read-only reproducible output
package. No preparation evidence here authorizes actual physical evaluation.
