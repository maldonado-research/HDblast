# Independent fixed-geometry stress review

Analytic review, 2 October 2026. No physical source, mode, response, or tail
sample was evaluated. The inherited proposal and its 21-check exact proof
were read from `HDBLAST_CHECKPOINT_20261002_CAUSAL_RESPONSE/theory/`.
Pure symbolic algebra in `derive_finite_cutoff_ward.py` independently checks
the finite-cutoff issue discussed below. No checkout file was edited.

## Physical operators and subtraction

The proposed physical density and pressure are those of the minimally coupled
field. Expanding `D v=v'-L v`, using
`delta|v'|²=-k² A`, and varying the explicit mass gives precisely

    2a⁴ delta e_k=3L² A-L A'+s/(2k),
    2a⁴ delta p_k=-(4k²/3+L²)A-L A'-s/(2k).

The explicit mass terms are essential. The zero first-order free canonical
energy does not set either physical stress component to zero. At small k,
`A=O(1/k)` and `A'=O(1/k)`; multiplication by the measure gives integrable
terms. Their apparently singular mode factors require no infrared cutoff.

An independent way to audit the subtraction is to expand the exact WKB
stress expression

    4a⁴ e_W = W+(k²+a²x)/W+L²/W+L W'/W²+W'^2/(4W³),
    4a⁴ p_W = W-(k²/3+a²x)/W+L²/W+L W'/W²+W'^2/(4W³)

through adiabatic order four, with `W=w+W2+W4`. Treat `s` as order two,
keep r fixed in w, and vary only W2/W4 and the explicit x. This reproduces
the stated u, v, j4, density, and pressure subtractions. In density the
`2Uu/w-sU/w²` terms cancel. In pressure the W4/source-second-derivative
terms survive. The proposed direct operators and finite contacts pass this
analytic audit.

## The finite-comoving-cutoff Ward identity survives exactly

Let `P=K/a`, `v=P/sqrt(P²+r)` and `r=2H²`. A closed expression for the
matched baseline is

    Q0,K = H²/(2pi²) [v²/(2(1+v))-v³/48-v⁵/16].

To derive it, at unit scale factor put `Omega=sqrt(P²+r)` and use

    U=-H²/Omega-3rH²/(4Omega³)+5r²H²/(8Omega⁵).

The baseline integrand, excluding `P²/(2pi²)`, is then

    1/(2P)-1/(2Omega)-H²/(2Omega³)
       -3rH²/(8Omega⁵)+5r²H²/(16Omega⁷).

Its logarithmic primitive terms cancel at r=2H². The displayed Q0,K
vanishes at K=0 and tends to H²/(12pi²). Its exact positive tail is

    Q0-Q0,K = H²(1-v)²(v+2)(3v³+3v²+10v+4)
                                  /[96pi²(v+1)],

which is O(K^-4). At fixed comoving K,

    v_dot=-Hv(1-v²),
    Q0,K_dot=-H³v²(1-v)²(5v⁴+15v³+21v²+23v+16)
                                  /[32pi²(v+1)].

The time dependence cannot be omitted when differentiating the density.
Write `C_K=v³(H²d-H d_dot)/(96pi²)` for its local contact and `A_K` for
the proposal's finite-cutoff trace remainder. Inserting the closed density
and trace into the cosmic-time Ward identity leaves the possible defect

    R_K=(Q0,K_dot/2+H Q0,K)d+C_K_dot+4H C_K+H A_K.

It is **identically zero**. There is no additional cutoff surface term for
these fixed-comoving-K definitions. The symbolic script also independently
checks the subtraction Ward identity for each individual comoving mode.
At unit a, a per-mode quantity is a^-3 times a function of physical P, so
the subtraction identity reduces to

    [D_t-H P partial_P]delta e_sub+3H delta p_sub
                                      =Q0,sub,k d_dot/2.

Both routes give zero exactly. This strengthens the inherited proof, which
checked continuum Ward and finite-cutoff contact primitives but did not
explicitly check this time-dependent finite-K cancellation.

For an independent numerical ledger, integrate direct density/pressure and
compare the boundary change of `a⁴ delta rho_K` with an independently
quadrature-evaluated integral of
`a⁴[Q0,K d'/2-L delta T_K]`. Do not define density from that ledger or
differentiate a formula that already encodes the same identity and present
the result as an independent evolution test. A mode-direct route and a
separately evaluated analytic-contact/variance route supply a meaningful
comparison. Freeze all ledger quadratures and derivative methods first.

## Derivative and stress-tail envelopes

For any smooth f with the required vanishing past derivatives define

    B_K[f]=[3|f|M²/2+(|f''|+integral_past^eta |f'''|)/4]
                                                       /(16pi² K²).

This bounds the canonical variance response tail `R[f]-R_K[f]`.
Differentiation of the combined cutoff expression, retaining the scale M,
gives exact identities

    q_K'=R_K[s']-Lv³s/(8pi²),
    q_K''=R_K[s'']-[2Lv³s'+L²v³(3v²-2)s]/(8pi²).

Consequently the continuum-minus-cutoff tails are bounded by

    b0=B_K[s],
    b1=B_K[s']+|Ls|(1-v³)/(8pi²),
    b2=B_K[s'']+[2|Ls'|(1-v³)+L²|s|(1+2v³-3v⁵)]/(8pi²).

These require source derivatives through fifth order. They are not obtained
by differentiating an absolute error inequality. They follow by first
differentiating the exact response and then bounding each resulting tail.

With `D0=Q0-Q0,K`, density is bounded by

    |delta rho-delta rho_K|
      <=(3L²b0+|L|b1)/(2a⁴)+|d|D0/2+|C_infinity-C_K|.

Pressure is bounded by

    |delta p-delta p_K|
      <=(b2+3|L|b1+3L²b0)/(6a⁴)+|d|D0/6
        +[|A_infinity-A_K|+|C_infinity-C_K|]/3.

The moment differences `Jn(infinity)-Jn(P)` are positive exact polynomials
in v from the proposal. Applying the triangle inequality to its three
trace-contact terms supplies a constructive bound for `|A_infinity-A_K|`.
The theory producer separately supplies source-derivative norm certificates.
These formulas cover removed momentum tails only; quadrature, time stepping,
roundoff, and cancellation need independent allowances.

## A signed post-pulse polarization diagnostic

For a nonnegative, nonzero compact s entirely before observation, set
`Delta=eta-u>0`. Every local source/contact term vanishes, and

    q=-(8pi²)^-1 integral s(u)/Delta du < 0,
    q'=(8pi²)^-1 integral s(u)/Delta² du > 0,
    q''=-(4pi²)^-1 integral s(u)/Delta³ du < 0.

Since L>0 in the expanding de Sitter patch,

    delta rho=(3L²q-Lq')/(2a⁴)
      =-1/(16pi²a⁴) integral s(u)[3L²/Delta+L/Delta²]du < 0.

This is a signed **first-order change relative to the same BD reference**,
due to coherent polarization. It neither asserts negative total physical
energy nor describes a heated fluid. Positive canonical excitation energy
begins at second order, so there is no conflict. Pressure has no universal
sign here. This is a useful nonlocal minimal-stress diagnostic that survives
every purely local source-contact convention after the pulse.

## Controls and remaining scope

Distinct negative controls should omit the explicit stress mass contacts,
omit the pressure W4 variation, freeze Q0,K when taking its time derivative,
replace finite contacts by continuum ones, vary r with x, reverse the box-Q
term, substitute improved conformal stress for minimal stress, and vary the
initial state. A compact initial Bogoliubov perturbation produces homogeneous
stress that a source-only calculation cannot reproduce; an altered occupation
state changes its retarded response. Passing Ward alone does not detect all
such wrong observables or states.

No mathematical blocker was found. The formulas remain a fixed-geometry
scalar-to-stress response, with no metric response or coupled stability claim.

## Symbolic development provenance

`EXACT_SYMBOLIC_REVIEW.log` records successful exact checks, with no physical
sampling. A first draft of an extra baseline-antiderivative check accidentally
typed its last coefficient as 5/8 instead of 5/16; the check failed and the
original source and actual rerun output are retained in `development/`.
The primitive, its Ward cancellation, and the inherited physical formulas
were unchanged. The corrected check differentiates the same baseline
integrand displayed above and passes. This was an algebra-audit coding typo,
not an evaluated physical run or an adjustment of the model.
