# Independent internal review of actual first-gauge fourth derivatives

30 October 2026. This review concerns the same retained actual three-amplitude
first gauge field, not a newly selected comparison field or an error bound.
No external specialist, originality or physical interpretation is asserted.

**Assessment: PASS for the stated conditional physical-side C4 theorem and
its quantitative ledger.** The separate [convergence audit](CONVERGENCE_AUDIT.md)
states the regularity assumptions and precise initial-iterate offsets.

## Ordinary derivative and initial-data algebra

The checker independently differentiates polynomials using
`dH/dt=1-H^2` and `dq/dt=-3Hq/2`. It recovers the ordinary third derivatives
of `H`, `J=1-H^2`, `HJ`, and `q=cosh(t)^(-3/2)`. In particular

```
q'''=[(39/4)H-(105/8)H^3]q.
```

This is an ordinary derivative, not a normalized third Taylor coefficient;
the latter differs by `3!=6`. The strip `0<=H<=3ell` with `ell=1e-6`
justifies the producer's rational temporal polynomial majorants. The static
whole-collar coefficient report supplies the positive `rho,a,C,C_z` caps.
Combining these with the independently established actual master fourth
jets gives the complete source `s_xxx` bound through the ordinary Leibniz
formula. It does not differentiate interval selections.

An independent noncommutative word algebra recursively differentiates the
original PDE at the initial face, preserving matrix multiplication order.
With zero original Cauchy data it proves

```
V_xx(0)=s0,
V_xxx(0)=P0*s0+s_x0,
V_xxxx(0)=s_rr0+P0*(P0*s0+s_x0)+2P_x0*s0
           +Q0*s_r0+R0*s0+s_xx0,
V_xxxr(0)=P_r0*s0+P0*s_r0+s_xr0.
```

The matrix product `P0*P0*s0`, radial-source terms and factor two on
`P_x0*s0` are retained. A separate exact word-algebra differentiation proves
the third-time commutator and its binomial coefficients. The compatible
flat initial shell wedge is an essential retained premise.

## Independent numerical reconstruction

Only authenticated upstream inputs and Python exact rational arithmetic are
used. No scientific producer module is imported. The audit independently
reconstructs the static source factors, pure third-time matrix/source caps,
same scalar contraction and inverse, all four time-differentiated C1
majorants `N0..N3`, and the stabilized free-wave Cauchy-data bounds.

It then bounds `T_xx,T_xr,T_rr`, where
`T=P V_x+Q V_r+R V+s`, and uses the physical PDE sign

```
V_xxrr=V_xxxx-T_xx,
V_xrrr=V_xxxr-T_xr,
V_rrrr=V_xxxx-T_xx-T_rr.
```

All five fourth mixed caps agree exactly with the producer, as do all
15 mixed derivative slots through order four for each of the three raw
coordinate amplitudes. The strict uniform scaled-coordinate bounds are

| Total derivative order | Uniform absolute component bound |
| --- | ---: |
| 0 | 3.37e-6 |
| 1 | 2.25e-6 |
| 2 | 3.82e-5 |
| 3 | 0.000971 |
| 4 | 0.0499 |

The audit runs 61 named checks, including eight deliberately wrong-formula
controls. Normal and optimized reports are byte-identical. These algebraic
and arithmetic checks support the accompanying analytic convergence proof;
finite controls alone do not prove an all-point derivative theorem.

## Scope

The theorem controls the actual gauge field on
`0<=x<=3,r>=0,x+r<=8`, with continuous physical-side derivatives and one-sided
boundary traces. It does not require or prove a smooth extension across an
artificial doubled shell. Physical `t,z` derivatives of total order `k`
require division by `ell^k`; subsequent conversion to other coordinates or
metric/scalar components requires the stated background and frame factors.

These actual-magnitude caps do not themselves bound the complete quadratic
bulk-source derivative, second-response C2 error, natural equation accuracy,
nonlinear realization or a cosmological origin. The complete source-chain
rule and any later error gate require separate checks. The original field,
source, coefficient system and physical boundary law remain fixed.
