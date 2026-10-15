# Matching the variation argument to the finite mixed-boundary evolution proof

15 October 2026 UTC. Addendum to
`INITIAL_GAUGE_SUBSIDIARY_AND_PARAMETER.md` (SHA-256
`74bce09080d4015b2eccf8c183b33f47b80516216a9547837e3e0e34251d7b03`).
This addendum records the finite-regularity passage explicitly. Its application
to the full contact theorem remains subject to the independent review of the
two evolution notes cited below. No earlier pinned proof is changed.

## 1. The actual dependency supplied by the evolution notes

The relevant arguments are Sections 2--6 of
`astra/QUASILINEAR_VARIABLE_MIXED_EVOLUTION.md` and Sections 2--8 of
`sol/INDEPENDENT_QUASILINEAR_MIXED_PROOF.md`. They construct smooth frozen
mixed-boundary solutions, control the variable scalar principal symbol, and
derive time/angular commutation plus recovery of normal derivatives. Their
tame estimates, specialized to uniformly bounded actual solutions in order
14, also control linearized solutions at integer orders 10 through 13.
Although the Astra note singles out orders 12 and 14 in its main application,
the same finite Leibniz count works at all intervening/lower orders: the
maximum total derivative order of the two commutator factors is at most 14
at order 13, and a factor of order at most seven is controlled in `L^infinity`
by the declared low index 10. The independent Sol proof obtains these orders
by its spatial Moser estimates and exact time-jet recursion. Neither argument
requires a new principal boundary estimate for the variations.

Let `X_r(V;t)=sum_{j=0}^r ||partial_t^j V(t)||_{H^{r-j}}`. Distinguish
`mathcal X_r=C_t X_r` (strong continuity of every declared jet) from
`mathcal B_r=L^infinity_t X_r`. The nonlinear solution has the strong
`mathcal X_14` regularity supplied by the smoothing comparison. The variation
argument needs only `W` bounded in `mathcal B_13` and `Z` bounded in
`mathcal B_12`, with strong continuity one order below. It does not need
top-order strong continuity of either variation. Accordingly the formal
high-order variation values in the original note should be read as these
bounded high-order estimates; its final asserted parameter topology remains
the strictly lower `mathcal X_10`.

## 2. Why the rough linearized coefficients fit those estimates

Use the notation `A(V) partial^2 V=F(V,partial V)` and `L_V` from the original
note. Its extra zero-order term is

    A_V[W] partial^2 V.

At order `r=13`, the time-commuted source norm differentiates this product
at most 12 times. The highest derivative of its `V` factor is therefore
of total order 14, exactly the nonlinear bound already available. Its `W`
factor is differentiated at most 12 times and is controlled by `X_13(W)`.
The remaining products are covered by the same finite Leibniz/Moser bounds.
For normal recovery, the corresponding pointwise source is differentiated
only 11 times, giving a lower variation norm `X_12(W)` and the known
nonlinear bound. The identity

    X_12(W;t) <= X_12(W;0) + integral_0^t X_13(W;s) ds

removes that instantaneous lower term before Gronwall. Thus the zero-order
coefficient does not require a fifteenth derivative of `V`. The first-order
terms `F_P partial W` have exactly the same order budget. The same argument
at `r=10,11,12` uses fewer derivatives.

The second-variation right side contains `W partial^2 W`. At order 12 its
highest time/source jet has 11 additional derivatives and uses at most 13
derivatives of `W`, already controlled. Its other worst term
`W W partial^2 V` uses at most 13 derivatives of `V`. Consequently the
second-variation estimate gives a uniform `mathcal B_12` bound. Constants
depend on the common nonlinear `mathcal X_14` bound, the fixed interval,
hyperbolic margins and finitely many initial curve derivatives. They are
independent of smoothing and angular truncation.

## 3. Constructing the actual variations, without an unproved rough solver

Smooth the **completed reduced initial data curve** by a fixed linear
operator and fixed support cutoff, retaining a smaller common zero shell
collar as in the evolution proof. Denote the resulting smooth initial curve
by `d_N(epsilon)` and its nonlinear solution by `V_N(epsilon)`. Smoothing
commutes with parameter differentiation of this initial curve. The exact
canonical constraints need not hold for these analytic approximations.
Their limiting initial curve is the original exact completed datum.

For each fixed smooth `V_N`, solve the two smooth linearized equations of
the original note with initial data `d_N'` and `d_N''`. The smooth frozen
solver applies to the scalar principal operator and its smooth bounded
first/zero-order linearized terms. The initial collar is zero for all
parameters. The nonlinear solution and the variation sources are therefore
zero on a common initial corner wedge, so their recursively determined
initial time jets obey every needed mixed boundary identity. No arbitrary
finite-regularity coefficient is smoothed independently of its boundary
geometry; using the actual `V_N` preserves `g^{za}|_b=0` and the corner
wedge exactly.

The preceding bounds give `W_N` uniformly in `mathcal B_13` and `Z_N`
uniformly in `mathcal B_12`. The nonlinear smoothing comparison gives
`V_N -> V` in `mathcal X_14`. Difference estimates for the linearized
equations then show

    W_N -> W in mathcal X_12,
    Z_N -> Z in mathcal X_11,

while retaining the respective bounded high norms by weak compactness.
For the first line the worst difference sources contain
`(V_N-V_M) partial^2 W_M` or `W_M partial^2(V_N-V_M)`.
The source order for an order-12 estimate is 11, controlled by the
`X_13` nonlinear difference and bounded `X_13(W_M)`. For the second line
the new worst source is `W_M partial^2(W_N-W_M)`, controlled in order 10
by the first line, giving order-11 convergence of `Z_N`. The initial
variation data converge in higher orders than these comparisons require.
Initial source jets are controlled by the same local recursion. Thus
`W,Z` are actual linearized solutions with the stated lower strong
regularity, rather than merely formal differentiated expressions.

## 4. The difference-quotient table does not assume top variation convergence

At two parameter values `epsilon` and `epsilon+delta`, the nonlinear
one-derivative-lower difference estimate gives

    ||V(epsilon+delta)-V(epsilon)||_{mathcal X_13} <= C |delta|.

Subtracting the first-variation equation from the nonlinear quotient
equation gives a residual whose worst term is

    [A(V(epsilon+delta))-A(V(epsilon))]/delta
         * partial^2[V(epsilon+delta)-V(epsilon)].

It is `O(|delta|)` at source order 11. Integral Taylor remainders of
`A,F` have the same bound. The order-12 linear estimate therefore gives
the quotient remainder `O(|delta|)` in `mathcal X_12`. This uses the
nonlinear `mathcal X_13` difference, not convergence of variations in
their top order. Comparing the first-variation equations gives
`W(epsilon+delta)-W(epsilon)=O(|delta|)` in `mathcal X_12` as well.

In the quotient of those first variations, the remaining worst terms
contain either

    partial^2[W(epsilon+delta)-W(epsilon)]

or

    partial^2[(V(epsilon+delta)-V(epsilon))/delta-W(epsilon)].

Their source order is at least 10, because both bracketed differences
are controlled in `mathcal X_12`. The declared order-10 target only
requires source order nine; this supplies a full additional derivative
of margin. Taylor coefficient remainders tend to zero in those same
product spaces. The quotient of first variations therefore converges
to `Z` in `mathcal X_10`. Comparing the second-variation equations
gives continuity of `Z` there, with the same margin. All comparison
initial jets converge by smoothness of the fixed initial curve.
Every boundary remainder is exactly zero, because the transformed
Dirichlet/Neumann rows are homogeneous and independent of the parameter.

This establishes the original conservative conclusion

    epsilon -> V(epsilon) is C^2 into mathcal X_10,

on the common controlled interval, once the reviewed order-14 evolution
theorem is integrated with the exact initial family. Composition with the
smooth inverse state flow preserves it. The same Sobolev embedding argument
as in the original note gives the piecewise spacetime `C^4` parameter
topology. No estimate assumes top-order differentiability of the nonlinear
solution map, top-order strong variation continuity, or a common spatial
`C^infinity` family.
