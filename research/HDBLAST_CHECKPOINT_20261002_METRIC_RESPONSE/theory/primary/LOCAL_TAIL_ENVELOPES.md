# Analytic omitted-band envelopes for the new metric contacts

These are algebraic inequalities to support a prospective registration. No physical source point, kernel integral, mode evolution, cutoff result, or trial gate was evaluated here. A numerical producer must use a separately frozen directed-arithmetic implementation and the registered source derivative budgets; this document by itself supplies no certified numbers.

Write `v=K/sqrt(K²+2L²)`, `0<v<1`, `t(v)=v²/(1+v)`, and use the exact coefficient inventories in `FORMULAS.json`. All expressions below omit the common `1/(2pi²)` until it is explicitly restored.

For every odd `n>=5`,

    ΔI_n=∫_K^∞ k²dk/(k²+2L²)^(n/2) > 0,
    ΔI_n <= K^(3-n)/(n-3).

This follows from `omega>=k`. It is valid for every positive K and L. It bounds only the omitted UV band. Its conservative value need not be close to the actual tail.

Moreover,

    0<1-v<=L²/K²,
    0<t(1)-t(v)<=3(1-v)/4,

because `t′(v)=v(2+v)/(1+v)²` increases on `[0,1]` and has maximum `3/4`. If the contact coefficients are `cq_n`, `dr_n`, `dp_n`, then

    2pi² |Cq∞-CqK| <= Σ_(n>=5)|cq_n| K^(3-n)/(n-3),
    2pi² |Er∞-ErK| <= (3/4)|Lh′+h″/4| L²(1-v)
                       + Σ_(n>=5)|dr_n| K^(3-n)/(n-3),
    2pi² |Ep∞-EpK| <= (3/16)|h″| L²(1-v)
                       + Σ_(n>=5)|dp_n| K^(3-n)/(n-3).

The finite coefficient tables stop at n9, n13, and n15 respectively. Coefficients at omega^-1 and omega^-3 are already included in the rational `t` block and must not be double-counted.

The baseline lower-order primitives are `Rlead=L⁴v²(2v+3)/[4(1+v)²]` and `Plead=-L⁴v²(4v+3)/[12(1+v)²]`. The exact primitive derivative can be bounded on `[v,1]` by directed rational interval arithmetic. Alternatively absolute coefficient bounds on their numerator/denominator derivatives give a simpler conservative envelope. Add the absolute high-order baseline coefficients times `ΔI_n`, restore `1/(2pi²)`, and retain the metric normalization factors. The complete component bounds then have the form

    bound(q_metric tail)=bound(q_mass[g] tail)
                         +2|h|L² bound(Q0 tail)+bound(Cq tail),
    bound(R_metric tail)=bound(R_mass[g] tail)
                         +4|h|bound(R0 tail)+bound(Er tail),
    bound(P_metric tail)=bound(P_mass[g] tail)
                         +4|h|bound(P0 tail)+bound(Ep tail).

The existing matched-stress omitted-band analysis can be applied to canonical forcing `g` only after deriving the appropriate `g` derivative budgets. The forcing is not the old bump source. For a compact `h`,

    g^(m)=4 Σ_(j=0)^m binom(m,j)(2)_j L^(2+j) h^(m-j)
           -2 Σ_(j=0)^m binom(m,j)(1)_j L^(1+j) h^(m-j+1)
           -h^(m+2),

where `(p)_j=p(p+1)...(p+j-1)` and `(p)_0=1`. This identity uses `L′=L²`. The lower forcing and local structure through g’s third derivative requires h derivative budgets through fifth order. That is insufficient for the full q_second omitted-band bound. Its inherited norm is `N2(g)=|g^(4)(eta)|+integral|g^(5)|`; therefore the complete bound requires g through fifth derivatives and h through seventh derivatives, with the corresponding certified supremum and L1/total-variation budgets. The new `metric-tail-theory/DERIVATIVE_CERTIFICATE.json` and its separately frozen directed tail module are the authoritative complete bound. On support `eta∈(-5,-3)`, `L<=1/3`, so sup/L1 budgets follow from the triangle inequality and the frozen h-jet budgets. For example,

    g‴=96L⁵h+60L⁴h′+12L³h″-2L²h‴-2Lh⁗-h⁗′.

Finite numerical quadrature error, refinement, ledger discretization, arithmetic error, and the pointwise accuracy of sampled source jets remain distinct. A directed interval UV tail enclosure does not certify them. If no frozen certified h/g derivative budget or directed local-tail implementation is available, the continuum comparison must stay descriptive rather than claim an interval-enclosed omitted-band acceptance result.
