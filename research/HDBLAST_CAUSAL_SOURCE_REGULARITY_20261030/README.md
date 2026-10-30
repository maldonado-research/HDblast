# Fourth master derivatives and a causally filtered numerical candidate

This checkpoint closes a quantitative fourth-derivative estimate for the
retained actual first master field. Separately, it defines a numerical
candidate with bounded second derivatives and explicitly bounds its natural
equation residual. These are different results with different scopes.

## Actual first-master result

On the retained triangle `0<=x<=3, r>=0, x+r<=8`, all five scaled mixed
derivatives of total order four of the same actual first master obey

```
|partial_x^a partial_r^b U_master| <342395.408,  a+b=4.
```

Its shell trace satisfies `|y''''|<3.235927`. The proof differentiates the
actual wave equation twice in time, uses its separable temporal potential,
and retains the reflected crossing and incoming-integral endpoint terms.
It requires no new third radial derivative of the potential. It does not
differentiate interval selections. See the
[derivation](theory/MASTER_FOURTH_BRIDGE.md), [exact report](reference/MASTER4.json)
and [independent internal review](review/INDEPENDENT_MASTER_FOURTH_REVIEW.md).

The actual first gauge field still needs its fourth-derivative ledger. Thus
this does not yet bound the complete differentiated quadratic bulk source
`F_x`. The first master here is not the formal second-order target below.

## A distinct auxiliary candidate with bounded second derivatives

Keep the actual second-coefficient target, first field, coefficients, physical
bulk/shell sources `F,G`, Robin matrix `M`, and PR100 quadratic-continued
initial data fixed. Only the new auxiliary finite approximation uses

```
k_epsilon(s)=30 s^2(epsilon-s)^2/epsilon^5, 0<=s<=epsilon,
F_epsilon(x,r)=integral_0^epsilon k_epsilon(s) F(x-s,r) ds,
epsilon=1/1024,
Z_epsilon=P[Bs_quadratic+K(La W0_quadratic+F_epsilon)+LG].
```

Negative-time `F` is extended by zero. All nonnegative past samples stay
inside the certified domain. This is component smoothing in the registered
coordinates, not a covariant physical averaging rule. The complete source
is coupled before output reflection.

The exact kernel variation bounds the derivative of the averaged source
without assuming a bound on actual `F_x`. The complete raw forcing cap is
`FH=17387513156573`. The new candidate change satisfies

```
X(Z_epsilon-W1_quadratic) <5.217e13,
X(U_target-Z_epsilon) <3.541e25.
```

The second inequality explicitly adds the source-change charge to the exact
previous error majorant. It retains the same readable conditional C1 upper
bound, not equality or improvement of the unknown true error. Previous point
tables are not inherited as exactly unchanged.

The differentiated principal-wave construction includes the initial
acceleration `u_s,rr+f(0)`, the actual `G_x`, and all actual coefficient time
derivatives. On `D+={0<=x<=3/2,max(0,x-5/4)<=r<=7/4-x}` it gives

| Candidate quantity | Uniform positive H-norm strict upper bound |
| --- | ---: |
| Value | 9.37e32 |
| Each first scaled coordinate derivative | 1.250e33 |
| Each second scaled coordinate derivative | 1.315e35 |

The second derivatives are continuous on the physical side, with one-sided
shell limits under the retained regularity premises. No C2 regularity across
the artificial doubled shell is asserted. The new initial target-error
forcing is generally `F(0,r)`, not zero. The physical corner remains flat.
See the [filter proof](theory/CAUSAL_FILTER_NOTE.md),
[second-derivative proof](theory/CANDIDATE_SECOND_NOTE.md),
[candidate report](reference/CANDIDATE.json) and
[independent review](review/INDEPENDENT_REGULARITY_REVIEW.md).

## Explicit natural-equation coefficient residuals

The geometric conversion retains the radial scale product rule, frame
connection, harmonic constraint derivatives, complete physical sources, and
`F_epsilon-F`. All 20 metric source identities and all 23 forcing
normalizations are checked exactly. The natural equations are

```
E_ab=Ric_ab-dphi_a*dphi_b-(2/3)V(phi)g_ab,
S=Box_g phi-V'(phi).
```

For the formal second parameter coefficient built with `Z_epsilon`, the
positive background-frame, probability-normalized angular L2 bounds are

```
||E_second[Z_epsilon]|| <2.918e43,
||S_second[Z_epsilon]|| <2.118e30.
```

These are very large absolute upper bounds. They do not demonstrate small
residuals, relative accuracy, or that the candidate satisfies the natural
equations. `E` is the stated trace-reversed Einstein–scalar equation, not
the Einstein tensor `G_ab`. These coefficient residuals are not residuals
of a finite-amplitude nonlinear spacetime. See the
[conversion and limits](theory/NATURAL_RESIDUAL_NOTE.md) and
[exact ledger](reference/NATURAL.json).

## Reproduction and limits

From a complete checkout, using only Python's standard library:

```
python -I -S -B research/HDBLAST_CAUSAL_SOURCE_REGULARITY_20261030/checks/run_checkpoint.py --repo-root . --output-dir /tmp/hdblast-source-regularity-fresh
```

Use a fresh output directory. The runner reproduces the producer and
independent-audit reports in normal and optimized Python, compares each
byte for byte with its reference, and checks the deliberately corrupted
master formulas. Inputs are authenticated by SHA-256 before use. These are
internal mathematical/computational checks, not outside specialist approval.

The [focused primary-reference comparison](theory/LITERATURE_COMPARISON.md)
discusses a 2026 version of a nonlinear vacuum boundary theorem and why its
boundary conditions cannot be imported here. It is not an exhaustive
literature or novelty review.

All results retain the registered background, exact first-field equations,
harmonic-gauge neighborhood, source identities, smoothness and compatible
corner premises. Complete actual gauge fourth jets, actual bulk `F_x`, C2
accuracy to the unsmoothed target, tangential constraint propagation and the
full compatible second-order problem remain open. No exact nonlinear family,
cosmological prediction or evidence of a higher-dimensional Big Bang origin
is established. Frozen manuscripts and older candidates are preserved.
No new Zenodo version or journal submission accompanies this checkpoint.
