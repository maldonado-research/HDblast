# A source-bound test of completion dependence

31 October 2026. Analytic scope check. This does not contradict the
producer's selected weighted curve. It explains why that selection is essential.

The statement "the same first tangent suffices" is false under the retained
initial-data construction. There are symmetry-preserving exact constraint
families with that first tangent and different initial TT second derivatives.

Let `C=Tplus+Tminus` in the actual frozen angular normalization. The retained
`HDBLAST_SCREEN_MAGNETIC_CURVATURE_20261026/review/symmetry/SYMMETRY_AUDIT.md`
records the all-point identity

```
Pi_TT((DY) tensor (DY))=-(2/3) C,
<Tplus,Tplus>=<Tminus,Tminus>=1/6, <Tplus,Tminus>=0.
```

Thus `C` is nonzero, with squared probability-sphere norm `1/3`. It is also
invariant under the **whole** selected coordinate-reflection group: `Y` is
invariant under every coordinate sign change, `(DY) tensor (DY)` is an
invariant polar tensor, and the round-sphere orthogonal TT projector commutes
with each isometry, including orientation-reversing isometries. No claim
that either individual chiral template is reflection invariant is needed.

Choose any nonzero smooth radial `f` compactly supported strictly within the
selected weighted annulus, and let the physical metric direction be

```
B_yy=B_ya=0, B_ab=rho(y)^2 f(y) C_ab,
```

with all other canonical sectors of `B` zero. On the spatial warped metric
`dy^2+rho^2 kappa`, this direction is trace free and divergence free in the
**full four-dimensional** spatial sense. In particular, the radial divergence
is minus `a` times its angular trace, hence zero; the angular divergence is
the sphere divergence of `C`, hence zero. The background angular Ricci block
is a trace and the background scalar gradient is radial, so

```
DR[B]=div div B-Delta tr B-Ric0:B=0,
B(dphi0,dphi0)=0, P_g B=0.
```

The linearized momentum constraint has no metric dependence at
`pi0=varpi0=0`; therefore `D Chat|x0[B]=0` in both constraints. The direction
is smooth, compactly supported and reflection invariant, with the same open
background shell collar as the original seed.

Now use the **same** weighted correction operator and inverse in the retained
Banach implicit-function argument, but replace its trial family by

```
x_tilde(e)=x0+e v+e^2 B/2+A u_tilde(e).
```

This is a smooth change of the prescribed parameter polynomial; it preserves
the mapping/positivity/support assumptions. The derivative with respect to
the unknown remains the identical isomorphism `L`. Hence the same IFT gives
an exact canonical constraint family, conditional on exactly the retained
analytic hypotheses, with

```
u_tilde(0)=u_tilde'(0)=0,
u_tilde''(0)=-L^-1[w^-2(D^2 Chat[v,v]+D Chat[B])]
            =u''(0).
```

Its first tangent is the original `v`, while its second metric derivative is
the original second derivative plus `B`. Every weighted-range metric
correction has zero angular TT, so

```
Pi_TT gamma_tilde2=rho^2 f C !=0.
```

This differs from the zero initial TT displacement of the specially selected
curve `x0+e v+A u(e)`. The harmonic completion and prescribed-slab shell
theorem apply after the same qualitative parameter shrinkage. The resulting
TT evolutions differ by a homogeneous TT solution with this nonzero initial
displacement. This argument does not assert a particular numerical difference
at the frozen event for an unspecified `f`; such an event claim would require
propagating that homogeneous contribution. It already disproves the inference
from a common first tangent to a common TT initial pair or history.

The valid transfer is therefore for the **specified fixed weighted-range
initial curve**, or another family whose completed TT initial pair has been
matched separately. It is not a theorem for every nonlinear completion of
the common first quadrupole. The producer's present theorem has this proper
restriction.

Sources: exact initial theorem, `Fixed weighted spaces and linear inverse`,
`Exact fixed-regularity conclusion`, and `Concrete commuting quadrupole
consequence`; metric-only adjoint/scalar-gradient normalization in
`HDBLAST_SECOND_ORDER_LOCALIZED_20261014/manuscript/SECOND_ORDER_LOCALIZED.md`;
angular TT identity and norms in the cited symmetry audit. These are analytic
and all-point angular identities, not an inference from parity alone.
