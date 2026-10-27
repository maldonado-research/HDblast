# Exact C² obstruction in the published 23-mode numerical candidate

The frozen all-mode candidate `V=B0+B1` is intentionally C¹. The existing
C¹ error certificate is unaffected by the observation below. This note proves
that one of its second coordinate derivatives has a jump at an interior point;
that jump is in the *approximation*, not in the unknown exact field, and is not
an Einstein-equation residual or evidence of a physical singularity.

The published representation in `DERIVATION.md` uses the exact smooth-in-the-
relevant-cell rational initial displacement and a continuous piecewise-linear
initial velocity `v` in the free wave `B0`. The feedback term `B1` uses the
constant matrices `A,B,D`, independent P1 profiles `u,du,v`, the reflected
Duhamel kernel, and the homogeneous Robin lift. Consider

```
(x*,r*)=(1/1000,1499/1000),   r*+x*=3/2,
q*=r*-x*=749/500,              r*>x*>0.
```

This is an interior point of the retained `D+`, outside the reflection
region. The P1 profile has the knot `3/2` with adjacent nodes `383/256` and
`385/256`; `q*` lies strictly between the preceding knot and 3/2. The published initial reconstruction sets the exact `K:0` displacement
identically to zero; hence it contributes no jump at either endpoint. In a
small neighborhood, therefore, the only possible jumps below come from the
`r+x` profile knot. The boundary lift vanishes because `x<r`.

Write `Δv'` for the right-minus-left jump in the radial slope of the P1
profile and similarly for `u,du`. The velocity part of the free solution is
`B0_v=(1/2) ∫_(r-x)^(r+x) v(s)ds`, so its xx-jump is `Δv'/2`. For B1, the
published function `angular_ridge.CoupledRidges.F` has

```
F(s) = ((A+B)(du(s)+v(s)) + D(u(s)+∫0^s v))/2.
```

The published non-reflected Duhamel formula contains `(x/2) ∫0^(r+x) F(s)ds`;
its other terms are C² across `r+x=3/2` in the selected neighborhood. Hence
its xx-jump is `x ΔF'/2` and the total jump is

```
Δ(V_xx) = Δv'/2 + (x/4)[(A+B)(Δdu'+Δv') + D Δu'].
```

For the `K:0` component, exact rational arithmetic on the published node and
matrix tables yields

```
Δ(V_xx) = 122949678051127779469013307270796835622755434515356118845205572291510230097770149650546701
          / 336999333339382997433337688587745383420464305281757156013795128115200
         ≈ 3.6483656164e20 != 0.
```

The sign is right-minus-left as x crosses x* at fixed r*. Both the free and
feedback contributions are included; cancellation has been checked rather
than assumed away. The amplitude and derivatives use the model's normalized
coefficient coordinates and do not represent detector strain or physical
energy. A pointwise classical second-order residual for this candidate cannot
be evaluated *on this characteristic* without a piecewise/weak convention or
a smoother candidate. This does not rule out residual computations away from
mesh characteristics, a valid weak residual, or classical smoothness of the
exact formal response.

## Reproduction

`check_candidate_kink.py --repo-root /path/to/HDblast` authenticates the frozen
profile, operator and initial-adapter bytes, checks the point, nodes and stencil, evaluates the
formula exactly, then independently obtains the same B1 jump from the frozen
`angular_ridge.py` by exact one-sided four-point differentiation. B1_x is
piecewise cubic on the chosen stencil, so this differentiation is exact.
Python's standard library is sufficient. Run normal and `-O` modes; explicit
checks remain active. It does not rerun the older physical source campaign.
