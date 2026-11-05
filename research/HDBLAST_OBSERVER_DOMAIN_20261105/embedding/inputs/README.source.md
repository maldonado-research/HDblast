# Conditional nonlinear shell-observable budgets

5 November 2026, America/Los_Angeles. Working notes and exact algebra checks.

This checkpoint makes the missing nonlinear estimates more explicit for the
specified proper-time, actual-area quadrupole moment. It applies established
determinant calculus, parameter-dependent ODE estimates and integral Taylor
remainders. It supplies neither a new Einstein–scalar existence theorem nor
evaluated nonlinear constants or a numerical amplitude window.

## Observable and retained assumptions

On the common, noncrossing intrinsic shell-observer cuts, retain the fixed
transported template `Y=n1^2-n2^2`, the full pulled-back scalar
`A_e=Phi_e^obs-phi_b`, and area density `dA_e=j_e dA_0`. The observable is

```
M_e(tau) = integral Y A_e j_e dOmega / integral Y^2 j_e dOmega.
```

The positive background area factor cancels. Parameter derivatives hold actual
proper time and transported labels fixed. The retained first tangent is
`A_0=0`, `A'_0=S(tau/rho_b)Y`, so `M_0=0`, `M'_0=S`. The full scalar and
cut metric include clock, embedding and angular-transport variations. Replacing
them by coordinate-time fields would require a separate justified join.

The earlier [nonlinear bridge][bridge] states the exact-family, preparation,
regularity, gauge and observer-domain premises. These remain assumptions here.
Its sufficient guard `|e| K(a)<106`, where `K(a)=sup|M''_e|`, uses the older
coarse central proper-time bracket
`I=[-rho_b*6145ell/4096,-rho_b*6143ell/4096]` and endpoint margins `63/53`.
It is not automatically a guard for the later, narrower linear bracket.

## What this checkpoint adds

1. [Metric and area budget](geometry/AREA_DENSITY_BOUND.md). For a positive
   three-dimensional cut metric, the quadratic determinant term obeys the
   sharp bound
   `|tr(A)^2/4-tr(A^2)/2| <= (5/4)||A||op^2`.
   The continuous proof covers the entire eigenvalue cube. Grouping correlated
   terms before taking absolute values gives a direct normalized-moment bound:

   ```
   K(a) <= (6/pi) chi [E2+6p E1+a E1(3q+9p^2)],
   p=P/lambda, q=Q/lambda,
   chi <= (Lambda/lambda)^(3/2).
   ```

   Here `lambda I<=G_e<=Lambda I`, `G_0=I`, `||G'_e||<=P`,
   `||G''_e||<=Q`, and `E_r=sup|partial_e^r A_e|` on the same admitted
   domain. The bound retains area and denominator terms without requiring
   `a sup|j'|<=1/2`. An optional angular-oscillation form cancels angularly
   uniform area rescalings exactly. Only the determinant coefficient `5/4`
   is claimed sharp; the final moment constants are sufficient bounds.

2. [Odd-family remainder](odd/ODD_REMAINDER.md). For the explicitly selected
   equivariant family, `M` is odd in the fixed affine scalar-amplitude
   parameter. Integral Taylor bounds can use a continuity modulus rather
   than treating `sup|M''|` as a fixed constant. A residual formulation
   subtracts `e S Y` before estimating the numerator, preserving its exact
   cancellation for arbitrary positive area density. The note distinguishes
   anchored envelopes, full two-point Hölder moduli, and stronger Lipschitz
   assumptions. `C^2` alone does not provide a cubic remainder or an effective
   numerical constant.

3. [Observer-flow interface](observer/OBSERVER_FLOW_BUDGET.md). First and
   second parameter variations of an augmented geodesic/Jacobi flow supply
   conditional bounds for the full scalar and metric pullbacks. This makes
   the clock, embedding and label-Jacobian costs explicit. The required
   coefficient, initial-data, chart, common-domain and flow constants still
   need evaluation.

These are complementary interfaces, not three independent proofs of one
complete nonlinear shell development. A source audit separates genuine retained
uniform variation estimates from numerical quantities known only at parameter
zero. The latter cannot simply be substituted as uniform nonlinear budgets.

## Reproduction and interpretation

The closed inventory and replay plan authenticate the included sources and
run three standard-library algebra-checker pairs in normal and optimized Python
modes. Exact finite checks accompany the written continuous arguments; they
are not proof-assistant certification or external peer review. Source/CLI
corruption controls test transport integrity separately from the mathematics.
See [REPRODUCE.md](REPRODUCE.md), [REFERENCES.md](REFERENCES.md), and
[NEXT_STEPS.md](NEXT_STEPS.md).

No nonlinear PDE simulation, physical amplitude calibration, incoming-source
selection, radiation production, or observational confirmation is supplied.
The earlier frozen papers and capsules are preserved. This incremental working
checkpoint does not include a new Zenodo edition or journal submission.

[bridge]: https://github.com/maldonado-research/HDblast/blob/ae272e3ec8f888f5ea53f3cd573272a870f68e92/research/HDBLAST_CHARACTERISTIC_ODD_TRACE_20261104/nonlinear/NONLINEAR_BRIDGE.md
