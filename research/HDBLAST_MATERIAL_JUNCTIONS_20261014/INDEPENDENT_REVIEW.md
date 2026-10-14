# Independent review of the prospective material-junction supplement

14 October 2026. Internal mathematical review by an independent OpenAI
Codex agent. This reviews the proposed equations and their scope; it is
not an existence theorem or external peer review.

The complete draft
`interface-source-20261014/MATERIAL_SECOND_ORDER_JUNCTIONS.md` was
subsequently read and compared with this independent derivation. Its
equations (5), (9)–(16), pullback conventions, tensor qualifications and
initial-boundary-value-problem limitations are accepted.

## Exact tensor equations

Use a separate outward-normal Gaussian collar on each side of the timelike
shell. The proper-distance coordinate increases toward the outward normal,
and the interior occupies the inward side of `r=0`. In family-adapted
coordinates,

```
g_±(ε)=dr_±²+h_±(ε,r_±,y),
K_±,ab=½∂_{r_±}h_±,ab,   p_±=∂_{r_±}φ_±.
```

The two induced metrics and scalars at `r_±=0` are exactly common.
For shell dimension four and a nondegenerate common induced metric, the
Israel equation

```
I_ab=Σ(K_ab−K h_ab)+σ(φ)h_ab=0
```

is exactly equivalent to `ΣK_ab=σ(φ)h_ab/3`: its trace gives
`ΣK=4σ/3`. This equivalence is nonlinear and tensorial; it does
not assume that the perturbation is scalar or that the metric is fixed.
The scalar equation is `P=Σp+σ′(φ)=0`.

With the derivative convention
`h(ε)=h₀+εh₁+ε²h₂/2+…` and
`φ(ε)=φ₀+εψ₁+ε²ψ₂/2+…`, the exact second derivatives give

```
ΣK₂,ab − σ₀h₂,ab/3 − σ₀′ψ₂h₀,ab/3
    = (2σ₀′/3)ψ₁h₁,ab + (σ₀″/3)ψ₁²h₀,ab,
Σp₂ + σ₀″ψ₂ = −σ₀‴ψ₁².
```

The factor two on `ψ₁h₁` and the third tension derivative in the
scalar source are required. No extra `K₁h₁` contraction belongs in
this covariant-tensor equation: the exact trace elimination has already
cancelled those terms. Such contractions would reappear if one expanded
the uneliminated mixed-index or trace Israel equation, with all inverse
metric variations then retained.

For an exact `Z₂`-symmetric family, each side has
`K_ab=σh_ab/6` and `p=−σ′/2`; dividing the displayed sum equations
by two gives the correct one-copy conditions. First-order opposite-copy
parity alone does not justify imposing this exact-family `Z₂` ansatz.

## Pullback and gauge qualification

The Gaussian formulas are legitimate after pulling back a hypothetical
moving-shell family and constructing its outward-normal Gaussian collars.
Their induced fields include embedding motion, normal normalization,
bending and evaluation effects. In a fixed ambient chart, the second
pullback derivative contains

```
g₂ + 2L_{χ₁}g₁ + L_{χ₂}g₀ + L_{χ₁}²g₀,
```

and the corresponding scalar terms; further Gaussian-coordinate changes
must also be included. One cannot replace the transformed `h₂,ψ₂`
by the bare fixed-chart perturbations and thereby discard these terms.

Normal Gaussian coordinates and background wave-map coordinates are
different gauge choices. The supplement's equations are geometric boundary
identities, not a proof that one can impose both gauge systems together
without further gauge-boundary analysis.

## Normal force and tangential compatibility

Let `E=G−T`, `ΔK=K_+−K_-`, `Δp=p_+−p_-`, using common
induced fields. The two Gauss equations give directly

```
2(E_+,nn−E_-,nn)
 = ΣK ΔK − ΣK_ab ΔK^{ab} − Σp Δp
 = σΔK + σ′Δp − I_abΔK^{ab} − PΔp.
```

Thus full bulk equations and the tensor/scalar junctions imply the normal
force equation, including through second order. The Codazzi equations
likewise give `ΣE_{na}=D_bI^b_a−P D_aφ`. These residual identities
must be differentiated with all traces, contractions and normal vectors
belonging to the common parameter-dependent geometry. They establish
compatibility identities; they do not solve an initial-boundary-value
problem.

If the first- and second-order initial data vanish near the cut, all
quadratic material boundary sources and their initial time jets vanish
there. This supplies the initial corner compatibility for the separated
packet. At later contact, those sources must be computed from the evolved
first-order fields. Their compatibility with gauge conditions and all
incoming hyperbolic modes is a separate unresolved task.

## Angular sectors and legitimate claims

The metric source is a full symmetric tensor. The product `ψ₁h₁`
can populate tensor sectors for generic first-order induced data. However,
when `h₁=f(y)h₀`, both displayed metric source terms are pure trace.
Therefore the supplement must not infer nonzero tensor boundary forcing
for a specific scalar packet without computing its material induced
metric. Tensor components generated in the bulk can still require a
second-order tensor sector satisfying homogeneous tensor junctions.

The scalar-source equation alone is insufficient to close the full
coupled problem. A valid next step is to combine the complete quadratic
bulk source, these material boundary sources, exact matching, localized
initial corrections, and gauge propagation in one forced hyperbolic
interface system, proving its well-posedness or identifying its remaining
compatibility obstruction. The proposed equations are useful progress
toward that target, not a completed shell-contact solution.

## Acceptance of equation (17): rescaled tensor boundary condition

The later addition in Section 5 is also accepted. The inherited
reflection-even scalar-master representative sets its first-order even
normal graph displacement to zero (the prior full-junction sufficiency
theorem explicitly has `ζ_e=0`). At the background shell `s=b²`,
every tangent satisfies `X·dX=0`. Thus the reconstruction
`δg=2F(sη−3X♭⊗X♭)` restricts to

```
h₁=(2b²F_b/Ω_b²)h₀.
```

A Gaussian collar extension that keeps the same intrinsic identification
cannot change this induced first-order tensor. Therefore both metric
source terms in (9) are pure trace in this representative.

For an exact two-copy reflection-symmetric continuation, the trace-free
part of the one-copy law (11) is consequently

```
(½∂_r−σ₀/6)h₂^TF=0   at the shell.
```

The background has `h₀(r)=ρ(r)²γ_dS` and
`a_b=ρ′_b/ρ_b=σ₀/6`. The trace-free projector using `γ_dS` is
independent of `r`; with `h₂^TF=ρ²t₂^TF`,

```
(½∂_r−a_b)(ρ²t₂^TF)|_shell
 = (ρ_b²/2)∂_r t₂^TF|_shell.
```

Hence every projected angular transverse trace-free tensor component
satisfies the necessary homogeneous condition
`∂_r t₂,TT|_shell=0`, as in equation (17). Here `r` is outward
proper normal distance and the notation `b` labels evaluation at the
shell; it must not be confused with using the ambient areal coordinate
as the normal coordinate. Imposing `∂_r h₂,TT=0` would miss the
warp-factor term.

This deduction supplies a necessary tensor junction law in the stated
material gauge and exact reflection-symmetric class. It does not imply
cancellation of the full bulk tensor forcing, solve that forcing, or
establish a tensor initial-boundary-value existence theorem.
