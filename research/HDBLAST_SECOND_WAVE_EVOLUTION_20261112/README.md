# Evaluated second-wave evolution in X13

This continuation follows merged PR #136. It evaluates the complete mixed
source hierarchy and the forced evolution of the selected second amplitude
derivative at zero. With the retained background, actual completed initial
data, homogeneous four-Dirichlet/twelve-Neumann shell rows and outer wall,

```
v=V2=partial_epsilon^2 V|0,
X13(v;t)=sum_(j=0)^13 ||partial_t^j v(t)||H^(13-j),
sup_(0<=t<=13/8000000) X13(v;t)
    <10^10200000 exp(10^152019).
```

The norms use the fixed product connection, full geometric sphere measure
and positive Frobenius/direct-sum metric on all sixteen fields. The result
includes the field H13 and velocity H12 bounds as well as all displayed
mixed time derivatives.

The first-wave W14 estimate now supplies both parts of the retained m=13
forcing cost: the evolving time/space source hierarchy through total order
twelve and the initial source jets through eleven. The forcing includes
`V1 D²V1` and all its coefficient-variation factors. Its normalized radial
source rows and its contribution to initial second-wave temporal jets are
charged explicitly. No first-wave order-fifteen extension is required.

[SECOND_WAVE.md](SECOND_WAVE.md) gives the forced energy, inhomogeneous
initial recursion, radial recovery and numerical bound.
[MIXED_SOURCE.md](MIXED_SOURCE.md) proves the complete source estimate.
[REPRODUCE.md](REPRODUCE.md) describes 438 exact source/formal/arithmetic
checks, normal/optimized comparison and seven corruption controls.
Eight immutable original public inputs are pinned in [SOURCE_PINS.json](SOURCE_PINS.json).

The actual nonlinear curve is retained as C2 into X10. Uniqueness identifies
its particular second derivative at zero with this stronger X13 solution.
The proof does not assert that the whole parameter curve is C2 into X13,
or that the actual finite-amplitude family stays in a quantified metric patch.
No angular quadrupole restriction is imposed on the second wave.

[LITERATURE.md](LITERATURE.md) compares three primary sources on forced
regularity, corner jets and sphere harmonic products, with coverage limits
recorded. [CONTINUATION.md](CONTINUATION.md) states the remaining work.
The analytic background and linear existence/uniqueness premises remain
conditional inputs. Finite exact checks support the written arguments;
they are not external peer review or physical confirmation. The very large
constants have no SI blast-energy calibration. A nonlinear amplitude interval
and the proper-time actual-area observable remainder remain open.
