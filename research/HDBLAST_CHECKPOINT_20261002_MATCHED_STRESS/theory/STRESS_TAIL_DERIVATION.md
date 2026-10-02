# Constructive cutoff bounds for the fixed-geometry stress calibration

Analytic preparation only. This document selects no numerical acceptance
tolerance and evaluates no source response, mode, quadrature, or stress.
The implementation and all numerical inputs must be publicly frozen before
its deferred interval functions are called. H=1, r=2, a=L=-1/eta,
L'=L^2; s=epsilon f is one of the two inherited smooth compact pulses.

## 1. Normalizations and exact finite-cutoff derivatives

Define q=a^2 delta Q/epsilon, R=a^4 delta rho/epsilon,
P=a^4 delta p/epsilon. The current normalization for the declared b=1,
r=2 quadratic mass law is J=a^2 delta j/epsilon=q+Q0 f/4.
All primes here are conformal-time derivatives; in particular q' is the
derivative of a^2 delta Q/epsilon, not a^2 delta Q'/epsilon.

Let M=a sqrt(r), v=K/sqrt(K^2+M^2), C=1/(8 pi^2), and

    R_K[g]=C integral_0^K dk [
       g(eta) k^2/(k^2+M^2)^(3/2)
       -2 integral_eta_i^eta g(t)sin(2k(eta-t))dt ].

R[g] denotes the matched removed-cutoff limit, with the same time-dependent
M. Since the source and all initial derivatives vanish,
`(I_g)'=I_(g')`, where I_g is the sine convolution. With
`A_K=asinh(K/M)-K/sqrt(K^2+M^2)`, exact differentiation gives

    A_K'=-L v^3,
    A_K''=L^2 v^3(2-3v^2),
    q_K=R_K[f],
    q_K'=R_K[f']-C L v^3 f,
    q_K''=R_K[f'']-C[2L v^3 f'+L^2 v^3(3v^2-2)f].

The continuum formulas replace v by 1. These scale-derivative contacts must
be retained when differentiating the matched variance.

## 2. Bounds for q, q' and q''

For j=0,1,2 define the normalized derivative budgets

    N_j(eta)=|f^(j+2)(eta)|+integral_eta_i^eta |f^(j+3)(t)|dt,
    B_j(K)=[3 M^2 |f^(j)(eta)|/2+N_j/4]/(16 pi^2 K^2).

Three integrations by parts, using the zero initial jets, give

    |I_g-g(eta)/(2k)| <= [|g''(eta)|+integral|g'''|]/(8k^3).

The subtraction inequality
`|k^2/(k^2+M^2)^(3/2)-1/k|<=3M^2/(2k^3)` then gives
`|R[g]-R_K[g]|<=B[g]`. The elementary momentum integral is
`integral_K^infinity k^-3 dk=1/(2K^2)`. Differentiation of a big-O remainder
is not used. Apply the bound separately to f, f' and f'', and retain the
exact local differences to obtain

    T0=B_0,
    T1=B_1+C L |f|(1-v^3),
    T2=B_2+C[2L |f'|(1-v^3)+L^2|f|(1+2v^3-3v^5)],

with `|q-q_K|<=T0`, `|q'-q_K'|<=T1`,
`|q''-q_K''|<=T2`. The last contact polynomial is nonnegative:

    1+2v^3-3v^5=(1-v)(3v^4+3v^3+v^2+v+1), 0<=v<=1.

Only derivatives through f^(5) are required. These are distinct derivative
budgets; the original variance tail cannot simply be reused for T1 or T2.

## 3. Exact reference and local-contact cutoff tails

The actual positive-reference background variance at the same cutoff is

    Q0,K=H^2/(4pi^2)[v^2/(1+v)-v^3/24-v^5/8].

Although the continuum Q0=H^2/(12pi^2) is constant, Q0,K is time dependent.
Its exact, positive, stable tail is

    DQ0=Q0-Q0,K
       =H^2(1-v)^2(v+2)(3v^3+3v^2+10v+4)/[96pi^2(v+1)].

It is O(K^-4). Here H=1 in the implementation. Define

    D_rho=L(3L f-f')/(96pi^2),
    c_rho=D_rho(1-v^3).

The exact finite-cutoff density contact is v^3 D_rho, so c_rho is its
continuum-minus-cutoff difference in the R normalization.

For the trace remainder, put

    d5=(1-v^3)/3,
    d7=(1-v)^2(3v^3+6v^2+4v+2)/15,
    d9=(1-v)^3(15v^4+45v^3+48v^2+24v+8)/105.

These are exactly `r^((n-3)/2)[J_n(infinity)-J_n(K/a)]`, for n=5,7,9,
where `J_n(P)=integral_0^P p^2/(p^2+r)^(n/2)dp`.
In the normalized variables the exact anomaly-contact tail is

    c_A=a^4(delta A_infinity-delta A_K)/epsilon
       =1/(2pi^2)[(f''-12L^2f)d5/16
            -(30L^2f+10L f')d7/32+70L^2f d9/64].

The pressure-contact tail is `c_p=(c_A+c_rho)/3`. This sum may cancel;
the interval implementation combines it before taking an absolute bound.
It does not bound the two large unrenormalized stress terms separately.

## 4. Stress and current tail budgets

The continuum stress closures in these normalizations are

    R=(3L^2 q-L q')/2+a^2Q0 f/2+D_rho,
    P=(q''-3L q'-3L^2 q)/6-a^2Q0 f/6
        +(f''-3L f'-11L^2f)/(288pi^2).

The finite-K closures replace q and its derivatives by q_K, Q0 by Q0,K,
and the local terms by their actual finite-K versions. Consequently

    B_rho=(3L^2 T0+L T1)/2+a^2|f|DQ0/2+|c_rho|,
    B_p=(T2+3L T1+3L^2 T0)/6+a^2|f|DQ0/6+|c_p|,
    B_j=T0+|f|DQ0/4,

bound `|R-R_K|`, `|P-P_K|`, and `|J-J_K|`, respectively.
These are explicit finite-point absolute error bounds on the omitted
momentum band. Before the source they give exact zero source-response tails,
even though the unrelated background reference tail DQ0 need not vanish.

The same linear coefficients propagate independently supplied errors in q,
q',q'', Q0,K and contacts. For example, finite-K density integration error
is bounded by `(3L^2 e0+L e1)/2+a^2|f|eQ0/2+e_local_rho` when those inputs
have actual enclosures. Empirical quadrature estimates remain estimates.
Additional source-jet evaluation errors require their own propagated terms.

A fixed acceptance rule may compare continuum/finite-K disagreement with
the computed B_rho or B_p plus separately frozen numerical allowances. It
must not silently replace the computed bound by the observed disagreement.
The bounds are O(K^-2) with finite explicitly computed pulse constants;
their actual K=256 values belong to the registered run. Any additional
hard maximum for those values must be frozen prospectively, with failure
preserved and no output-driven relaxation or cutoff extension.

## 5. Exact finite-K Ward identity

The direct minimally coupled mode observables and the full inherited
second-/fourth-order subtractions satisfy the per-mode source identity.
Integrating at a fixed comoving K therefore gives exactly

    delta rho_K'+3L(delta rho_K+delta p_K)=Q0,K delta x'/2.

Q0,K is the same matched finite-K reference above. The time dependence of
Q0,K and v^3 does not require an additional freely selected local contact:
its derivatives cancel against the fixed finite-K trace remainder. Using
continuum Q0 or A in that exact finite-band comparison would instead leave
a cutoff artifact. The separately prepared exact Ward audit verifies this
cancellation using the explicit Q0,K and J5,J7,J9 formulas.

The displayed stress tail bounds do not bound a time derivative of stress.
A continuum Ward residual needs further derivative control. A finite-K Ward
test or integral energy ledger must use independent numerical derivative or
integration accuracy controls; defining pressure or the energy derivative
through Ward would make the check true by construction.

## 6. Computable pulse constants and implementation boundary

`pulse-norms/pulse_derivative_budgets.py` writes f^(n) as a rational
polynomial times the normalized bump, for n=0,...,5. For each N_j it isolates
all critical points of f^(j+2) using an exact polynomial in u^2 and complete
root counts on (0,1). Total variation between those critical points equals
the required integral of the absolute next derivative. Exact rational root
brackets are refined by 256 bisections; primitive values are enclosed with
60-digit directed interval arithmetic. Every supplied partial interval is
ordered relative to the rational observation time without approximate sign
decisions. This gives a constructive norm bound rather than labeling a
quadrature estimate rigorous.

A separate exact rational certificate bounds the full total variations.
Because the primitive f^(j+2) vanishes at both pulse endpoints,
`|f^(j+2)(eta)|+TV(past)<=TV(full pulse)`. It proves

| Source | Global N0 upper | Global N1 upper | Global N2 upper |
| --- | ---: | ---: | ---: |
| B | 97 | 2926 | 161866 |
| uB | 96 | 2760 | 154078 |

These certificates use rational root enclosures and Taylor inequalities for
the exponential; they evaluate no source sample or floating exponential.
Using their common upper envelope, `|f|<=1`, `|f'|<4`, `|f''|<49`,
`L<=2/3`, `M^2<=8/9`, and pi>3 gives the additional rational certificates

    B_rho(256)<3e-5,  B_p(256)<1e-3,  B_j(256)<3e-6.

For example `|B'|<=8/e<3` and `|(uB)'|<=|B|+|B'|<4`;
the endpoint-zero second derivative satisfies `sup|f''|<=TV(f'')/2`.
`prove_global_stress_ceilings.py` checks the remaining rational inequalities
and records them in `GLOBAL_STRESS_CEILINGS.json`. These are conservative
analytic feasibility statements, not a changed numerical precision goal or
new acceptance tolerances. The deferred runtime functions return the tighter
actual interval budgets at each registered observation.

`stress_tail_bounds.py` carries all source-jet, derivative-budget, geometric,
pi, reference and local-contact operations through the same interval context,
then advances binary64 upper endpoints toward positive infinity. Its public
function is `stress_tail_bounds(source,eta,K)`. It returns q, q_prime,
q_second, rho, p, current, Q0_difference and anomaly bounds with diagnostic
intervals and explicit normalizations. Imports do not evaluate any source.

This preparation supplies error controls for a fixed-geometry linear
stress/current calibration. It makes no coupled shell evolution, shifted-root
kernel, quantum-stability, physical particle-yield, heating, or novelty claim.
