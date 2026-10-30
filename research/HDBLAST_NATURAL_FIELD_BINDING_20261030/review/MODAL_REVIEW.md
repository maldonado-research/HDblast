# Independent implementation review: natural residual-corrected Bianchi

30 October 2026. Internal implementation review only; no map certificate is
claimed in this note.

The retained `geometric_sources.py` frame is radially parallel. Its six slots
are y, one temporal slot, and four redundant ambient tangent slots on S^3.
Their contraction signs are `(1,-1,1,1,1,1)` and the ambient projector trace
is three in the sphere quotient. The frame derivative formula gives zero
connection along y. Along time and sphere directions it must act on every
covariant slot, including a prior derivative slot. Ambient differentiation
must project before adding the sphere slot terms.

An efficient independent linear construction does not need to materialize
all rank-five third derivatives. With every index of P lowered, define

```
P[c,a,b] = (nabla_a q_bc+nabla_b q_ac-nabla_c q_ab)/2,
Ric1[a,b] = sum_c sign[c]*(nabla_c P[c,a,b]-nabla_b P[c,a,c]),
E1[a,b] = Ric1[a,b]-p delta[a,y] deta[b]-p delta[b,y] deta[a]
            -(2/3)*(Vp eta g0[a,b]+V q[a,b]).
```

Differentiate the resulting rank-two E1 with the same covariant operator.
This constructs the first residual from the connection variation, without
using the quadratic map or rearranging the target identity. All ordinary
metric profile derivatives through total degree three remain free; only
scalar jets through degree two are required.

With `C1_lower[d]=sum_a sign[a]*P[d,a,a]`, the operator correction is

```
deltaB[b] = sum_(a,c) sign[a] sign[c] q[a,c]
                 *(-nabla_a E1[c,b]+(1/2)nabla_b E1[a,c])
              -sum_d sign[d] C1_lower[d] E1[d,b].
```

The scalar linear residual offers a short independent diagnostic. For the
first scalar eta=sY,

```
S1/Y = syy+4a sy-(stt+3theta st+8s/c^2)/rho^2
          -py A-ap(-D+3L-8M)-p C1_y-Vpp s.
```

In particular, the temporal contribution in the Hessian contraction has
`-D`, since the physical temporal background Hessian equals `-ap`.
Independent harmonic covector amplitudes are

```
C1_y = (Ay+Dy-3Ly+8My)/2-Bt/rho
         +a(4A+D-3L+8M)-3theta B/rho-8C/(rho c),
C1_0 = By-(At+Dt+3Lt-8Mt)/(2rho)
         +5aB-3theta(D+L-8M/3)/rho-8E/(rho c),
C1_a = Cy-Et/rho+5aC-4theta E/rho
         +((-A+D-L)/2-2M)/(rho c).
```

Before the quadratic join, `B0(E1)+p*S1*e^y=0` is a useful independent linear
Bianchi precheck of the geometry, background chain rules, and coherent
metric third jets. It is a diagnostic, not a replacement for the requested
quadratic identity.

For a canonical background quotient, eliminate `theta^2` by
`1-c^-2`; theta then has exponent zero or one and c remains Laurent.
The de Sitter relation is stable because
`Dt(c^-2+theta^2-1)=-2theta*(c^-2+theta^2-1)`.
The radial background substitution also satisfies

```
Dy(6/rho^2+p^2/2-6a^2)
 = -12a/rho^2+p(Vp-4ap)-12a(-rho^-2-p^2/3)
 = Vp*p.
```

The first metric uses an untraced angular Hessian:
`J=L-8M/3,H=M`. The registered quadratic maps already use J/H and require
no second conversion. Both TT sectors are killed by B0, so an eight-row
check cannot validate those source components. If coefficients are extracted
by angular inner products, retain an exact angular reconstruction check;
projection alone would miss a remainder outside the registered modal span.

No first field equation or harmonicity is imposed by any construction above.
Deleting deltaB requires E1 and its first prolongation to vanish; deleting
`2S1 deta` requires S1=0. C1=0 alone supplies neither fact.
