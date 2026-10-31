# Independent reduced-source sign and premise review

31 October 2026. Private, read-only review of the retained conventions. No
scientific producer was imported or executed for this review. This note is an
analytic review, not an executed coefficient certificate or external peer
review. No sign or missing-premise defect was found in the proposed identities.

Let `g(e)=g0+e q`, let `P^d_ab` denote the first connection variation, and let
`C^d(e)=g(e)^{ab}(Gamma(g(e))^d_ab-Gamma(g0)^d_ab)`. All index movements below
use `g0`. The inverse expansion is `g(e)^-1=g0^-1-e q^sharp+O(e^2)`, and
the connection difference is

```
Gamma(g(e))^d_ab-Gamma(g0)^d_ab
  = e P^d_ab-e^2 q^{dc} P_cab+O(e^3).
```

Writing `C1_b=g0_bd C'(0)^d`, the upper gauge Hessian, lowered by the fixed
background metric, is therefore

```
QC_background_b = -2 q^{ac}P_bac-2 q_b{}^d C1_d.
```

The product rule for the full lowering `g(e)_bd C(e)^d`, with `C(0)=0`, gives

```
QC_actual_b = QC_background_b+2 q_b{}^d C1_d
            = -2 q^{ac}P_bac.
```

Thus the factor two and the Lorentzian contraction are compulsory. In the
stored frame `q_b{}^d C1_d=sum_d signature[d] q_bd C1_d`; a positive Euclidean
contraction changes the temporal term. The same point applies to
`PC1_ab=2 sum_d signature[d] P_dab C1_d`.

For the general reduced tensor, the second variation of the gauge term is

```
[sym nabla(g(e)) C_lower(g(e))]''(0)_ab
  = sym nabla0 QC_actual_ab-2 P^d_ab C1_d.
```

The connection acts on the lower covector and enters with a minus sign.
Consequently `QEred_general=QE-gauge_gradient+PC1`. The registered
`QEred=QE-gauge_gradient` is its `C1=0` simplification. This value-level
simplification itself uses no derivative of `C1`.

The retained natural identity is

```
B0(QE)+p Qphi_general e^y+2 deltaB_q(E1)+2 S1 d eta=0.
```

Substitute `QE=QEred_general+gauge_gradient-PC1` and apply the linear `B0`.
The half-normalized symmetrization gives

```
B0(sym nabla Z)_b
 = (Box0 Z_b+[nabla^a,nabla_b]Z_a)/2
 = (Box0 Z_b+Ric0_b{}^c Z_c)/2.
```

The resulting general identity has `+(1/2)L0 QC_actual-B0(PC1)` on its
left side. Solving for `B0(QEred_general)` changes these to
`-(1/2)L0 QC_actual+B0(PC1)` on the right. Neither residual correction
disappears from generic free first jets.

The background is scalar-Einstein, with
`Ric0_ab=p^2 delta_a^y delta_b^y+(2/3)V g0_ab`. Raising the contracted
index gives

```
Ric0_b{}^c Z_c=(2/3)V Z_b+p^2 delta_b^y Z_y.
```

In particular, the temporal lower Ricci entry and the inverse temporal
metric both carry a minus sign; their product is positive. There is no extra
temporal sign in the endomorphism. The radial `p^2` term must be tested on a
covector with nonzero radial component. The warped-frame curvature agrees
with this Einstein expression only after using the algebraic constraint
`V=6/rho^2+p^2/2-6a^2`, as well as `a_y=-rho^-2-p^2/3`.

The exact sufficient premises for specialization are:

| Operation | Sufficient formal-jet premise |
| --- | --- |
| Remove `2 deltaB_q(E1)` | `E1=0` and `nabla E1=0` |
| Remove `2 S1 d eta` | `S1=0` |
| Replace `Qphi_general` by `Qphi` | `C1=0` |
| Replace `QEred_general` by `QEred` at value level | `C1=0` |
| Replace `B0(QEred_general)` by `B0(QEred)` or remove `B0(PC1)` | `C1=0`, `nabla C1=0` |
| Replace `QC_actual` by `QC_background` at value level | `q.C1=0`; `C1=0` suffices |
| Exchange the two QC fields inside `L0` | `C1=0`, `nabla C1=0`, `nabla^2 C1=0` suffice |

These are sufficient conditions, not claimed necessary conditions; special
products can vanish with weaker premises. The first-natural master identity
must hold in a neighbourhood to supply its first prolongation. Harmonicity
is a separate neighbourhood identity. With both premises the conditional
eight-row identity is
`B0(QEred)=-p Qphi e^y-(1/2)L0 QC_actual`.

A coherent counterexample to value-only QC exchange can be constructed
without inventing independent covariant `C1` jets. Set `q=A(y)Y dy^2`, with
all other metric profiles zero. Direct first connection contraction gives

```
C1_y=(A_y/2+4aA)Y,
QC_actual-QC_background=(A A_y+8aA^2)Y^2 dy.
```

At a point choose `A=0`, `A_y=1`, `A_yy=A_yyy=0`, and every ordinary time
jet zero. These are coherent ordinary profile jets. The QC difference then
has value zero, radial derivative `Y^2`, and second radial derivative
`16aY^2`. For a lower radial covector `F dy`, the rough wave radial entry is
`F_yy+4a F_y-4a^2 F+rho^-2 Box_dS4 F`; hence

```
L0(QC_actual-QC_background)_y=20aY^2.
```

For example `rho=p=1`, `a=1/2`, `a_y=-4/3`, `V=5` is compatible with the
retained background algebra, and the displayed nonzero witness is `10Y^2`.
This demonstrates that equality of QC values at one point is insufficient
for exchange inside the second-order operator. This is an analytic witness,
not a claim that a corrupted subprocess has already been run.

The source-map audit concerns the retained seven-profile first metric
ansatz and one scalar profile, with free ordinary jets at the required
orders. It does not establish an arbitrary full spacetime metric ansatz.
The proposed finite binding supports an already retained analytic subsidiary
argument; it proves no new PDE existence or uniqueness theorem, accuracy
improvement, nonlinear family, or Big Bang origin statement. Implementing
the actual numerical third-jet adapter remains a separate task.

Inspected source bytes:

| Source | SHA256 |
| --- | --- |
| `/workspace/hdblast-natural-replay-30n/NEXT_REDUCED_BIANCHI_BINDING.md` | `18407b13d4a27fff9aa63f362416cab0406bb1648f3d734c23601b828a908ed8` |
| `/workspace/hdblast-reduced-publish-31a/research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/geometry/indexed_geometry.py` | `7fff2472bacba8b1a79e8b6f9bf29ef46d58a06d987439e51394584980135bd4` |
| `/workspace/hdblast-reduced-publish-31a/research/HDBLAST_NATURAL_FIELD_BINDING_20261030/theory/NATURAL_BULK_BIANCHI_NOTE.md` | `7a04a0d7e739d83e345c1139bf35256b1a8e2f522e9168260d4c94c00da81a30` |

Read commands were `cat` on the three named files and `sha256sum` on those
same paths. No scientific executable closure was loaded by this note.
