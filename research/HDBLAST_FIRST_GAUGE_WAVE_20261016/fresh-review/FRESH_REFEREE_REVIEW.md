# Fresh scientific review of the first relative-harmonic gauge wave

16 October 2026. **Accept the conditional first-order gauge-wave checkpoint
within its stated scope. No blocking mathematical defect remains in the
reviewed files.** The exact scientific corpus is authenticated by
`FRESH_REVIEW_PINS.json`. This is a separate internal AI mathematical review,
not external specialist peer review or a formal verification of every proof.

The accepted conclusion is a unique raw-amplitude solution
`U=(f,b,d)` on `0<=x<=3, r>=0, x+r<=8`, with `ell=10^-6`, zero Cauchy data,
Dirichlet f and Neumann b,d, a C1 contraction below 0.011412, continuous
one-sided mixed derivatives through order three, and the stated continuous
error and five signed point enclosures. The smooth retained background and
actual integral-defined master/contact solution remain hypotheses. Interval
cells alone do not construct that master or permit arbitrary differentiable
selections. The earlier complete existence arguments are inherited, not
re-proved by this review or by the new replay programs.

## Operator, actual forcing, and domain

I checked the raw-amplitude convention against the optical conversion
`D_y=-rho^-1 D_z`, including the different radial drift coefficients 3 and 5,
the off-diagonal couplings, and the scaled factors ell and ell squared.
The reviewed operator acts on the contravariant amplitudes in
`xi^y=fY`, `xi^t=bY`, `xi^i=d D_kappa^iY`; no hidden radius rescaling is used.
The construction uses the actual source, with

`u=(alpha/rho^2) cosh(t)^(-3/2) w`.

An earlier operator draft and symbolic source identity omitted the evolving
`cosh(t)^(-3/2)` factor. I reported this defect. The final operator note and
checker restore it in both source terms and include a targeted wrong-factor
control. The coefficient and solver source constructions already retained
it. The final proof and the numerical inputs now agree.

The whole-collar argument is necessary: the older prose estimate on
`0<=z<=8e-10` does not establish the required `0<=z<=8e-6` region. My separate
Fraction checker recomputes the potential derivative by interval Horner
evaluation, the four optical ODE images, and the first-exit tube from the
original shell box. For each variable v it checks, exactly,

`[v_b^- ,v_b^+] + [0,8ell] F_v(B) subset interior(B_v)`.

It records both strictly positive face margins and matches the coefficient
certificate. The ODEs used are
`rho_z=-rho^2 a`, `eta_z=-rho p`,
`p_z=-rho(U'(eta)-4ap)`, and `a_z=1/rho+rho p^2/3`.
Thus the full collar follows by first exit for the retained regular branch,
not by extending the earlier prose interval without proof. The actual master
source derivatives use its authenticated C3 bounds and time conjugation.
Finite propagation supplies the required source extension; every predecessor
of the stated causal triangle remains in the certified region. No additional
inner boundary condition has been introduced.

## Existence and boundary regularity

For `U_xx-U_rr=A U_x+B U_r+D U+s`, the physical mixed retarded kernels give
`||K F||<=X^2||F||/2` and each first derivative at most `X||F||`.
Consequently, in
`N_i=max(2||U_i||/X,||U_ix||,||U_ir||)`, the nonnegative majorant is
`M=X(A*+B*+(X/2)D*)`. I independently checked this weighting, all entries
of the reported matrix, and `(I-M) R=I` with `R>=0`. Its maximal row sum is
approximately 0.011411030502. A second, coarser derivation directly from the
outer background box gives 0.011752418218, still below 0.012. The analytic
contraction therefore constructs the actual solution on the full triangle;
the finite arithmetic checks only instantiate its inequalities.

I checked the reflection derivatives, including the Dirichlet second
spatial trace `-1_(x>r) F(x-r,0)` and the Neumann third spatial trace
`-1_(x>r) F_r(x-r,0)`, with their differentiated companions. These terms
cannot be removed by asserting smooth odd/even physical coefficients.
The noncompatible constant Dirichlet and linear Neumann forcings are useful
counterexamples. Here the actual source and every iterate vanish in the
initial wedge `x+r<1`, which supplies the requisite corner compatibility.
The conclusion is physical one-sided C3 regularity at the shell, not a
smooth parity extension of every coefficient across it.

The C3 proof establishes convergence, not only bounded differentiated
iterates. With `V_n=partial_x U_n` and `W_n=partial_x^2 U_n`, the same C1
contraction acts at each level; the new forcing depends on already controlled
lower levels. Their difference estimates have summable polynomial-times-
`q^n` tails. The finite initial-data offset is correctly retained:
`W_1` has data `(s0,s_x0)`, whereas `W_n`, for `n>=2`, has
`(s0,A(0,r)s0+s_x0)`. The coefficient here is A, not A_x.
The limits U, U_x and U_xx in C1, together with
`U_rr=U_xx-F`, `U_xrr=U_xxx-F_x` and `U_rrr=U_xxr-F_r`, yield all continuous
mixed partials through order three.

The higher-norm certificate uses those same equations and resolvent. I
checked the free-data maximum, the factor two in all first-coefficient
commutators, the field-value factor 3/2, `B U_rr` and `s_r` in the third
radial recovery, and the conversion `ell^(-i-j)`. I also read and replayed
the separate independent sixteen-term Neumann-plus-tail audit. It verifies
all thirty mixed component bounds and twelve displayed total-order caps
without importing producer code or computing the producer's inverse.
The large unscaled upper bound `9.59747e14` at order three is a conservative
coordinate derivative bound on a scale ell; it is neither a measured
derivative nor evidence of a physical instability. These estimates are not
angular, orthonormal-frame, proper-y, covariant-vector or energy norms.

## Comparator, actual error, and interpretation

I derived the free-wave comparator from the commutator identity
`(partial_x^2-partial_r^2)(r W_disp)=-2 partial_r W_disp`.
The odd displacement is `P_o(s)=s Q(|s|)`, not `s Q_o(s)`. The stated formula
has exactly zero Cauchy data and the physical shell rows, including through
the reflected characteristic. Its source is an auxiliary free wave.
The error estimate compares it to the actual retained master using all
three missing contributions: initial velocity, the authenticated boundary
derivative, and the retarded potential term. The exact boundary derivative
bound is not replaced by the much smaller boundary value bound.

The residual-to-error equation and positive resolvent imply the displayed
global amplitude/error bounds. My checker independently verifies all five
actual signs with a coarser scalar contraction and directly bounded first-
column couplings; it does not use the solver's inverse or point intervals.
Their signs are `+,-,+,+,-` at `(x,r)=(1/2,1),(1/2,2),(3/4,3/4),(1,1/2),
(2,1/2)`. Thus the signed statement concerns the forced solution, not merely
the free comparator. The comparator's continuous residual, rather than a
collocation grid, supplies the all-point error argument.

The vector is tangent to the shell. Its tangential values are free; Neumann
conditions do not set them to zero. Geometric linear shell laws are preserved
with the corresponding pullback of fields and shell data. The gauge identity
then removes the first relative-harmonic defect. This does not provide the
complete transformed quadratic bulk/shell sources, a second-order evolution,
an exact nonlinear family, or a nonlinear remainder. Coordinate amplitude
signs are not force measurements or an invariant physical second response.

The operator note's optional second-variation covariance identity is an
identity for sufficiently smooth fields, or a formal identity. Its term
`G=2 L_xi h+L_xi^2 Phi0` contains second derivatives of xi, so separately
evaluating the classical second-order operator `E'[G]` can require fourth
derivatives of xi. The certified C3 caps alone do not supply that regularity
or a norm for this separate term. I concurred with the final clarification;
the identity is preparation for later work, not an additional quantitative
conclusion of this checkpoint.

The supplied bounded literature comparison correctly distinguishes this
direct construction from scalar quarter-plane formulas, displacement-error
estimators, and the broader generalized-harmonic Einstein boundary problem.
Standard gauge-reduction identities do not by themselves certify this mixed
shell problem or the later Einstein-scalar evolution. I reviewed the supplied
comparison and access limits, but performed no new external literature
search. No priority or exhaustive-literature conclusion is warranted.

## Review evidence and replay

`check_fresh_referee.py` imports only the standard library and no producer
module. It authenticates the exact pinned scientific corpus, independently
recomputes the full-collar first-exit inequalities, checks a coarse contraction,
verifies the tight matrix certificate by multiplication, and independently
encloses the five actual signs. Its normal and optimized outputs are identical.
I also replayed the coefficient certificate and independent higher-norm audit;
their complete reports match the registered normal reports.

From an assembled capsule, choose new external output files and run:

```sh
python -I -S -B fresh-review/check_fresh_referee.py --root . --output /tmp/new-fresh-review-normal.json
python -I -S -B -O fresh-review/check_fresh_referee.py --root . --output /tmp/new-fresh-review-optimized.json
```

The pin set covers the scientific files and inherited input evidence, not
every earlier proof by a fresh derivation. It excludes this review directory,
the final manifest, common replay runner and top-level roster to avoid a
recursive receipt. Those packaging/publication objects need their own final
authenticated replay. Finite identities, controls and hashes supplement the
analytical arguments; their counts are not theorem or peer-review counts.
