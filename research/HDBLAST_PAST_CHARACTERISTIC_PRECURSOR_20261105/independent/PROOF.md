# Independent bounds for a short-past master precursor

5 November 2026 UTC. This companion derivation establishes conditional
coefficient and characteristic bounds for the registered first linear master.
The checker performs rational arithmetic on selected retained inputs; it does
not execute an inherited campaign, numerical evolution or oracle.
`recompute.py` uses only standard-library rational arithmetic and five pinned
retained sources. Its normal/O reports agree byte for byte.

## Static collar and the sign of the potential

The width premise is the later `8ell` static collar, not the older `8e-10`
collar. The checker uses the **raw shell enclosure directly**. On the declared
outer rectangle it freshly encloses the four vector-field components and checks
`I_shell+[0,8ell]F(I_outer)` lies strictly inside every outer interval. It does
not reuse the wider report's stored flow intervals. Smooth branch existence,
the validity of the shell box and the standard continuation argument remain
analytic premises.

The model polynomial gives `U3>0`, so `U2` increases with eta; its
exact endpoint values lie in `(31/10,25/8)`. Therefore U1 increases and remains
negative. Exact endpoint division by the negative p interval gives
`7/4<h=U1/p<7/3`. For

```
Qpot=2h^2-10ah+35a^2/4-U2-p^2/6,
```

the h derivative `4h-10a` is positive and the a derivative `-10h+35a/2` is
negative throughout the box. Joint endpoint evaluation gives

```
17246399951/15000000000 < Qpot < 38208223/7200000,
1 < Qpot < 16/3.
```

The conjugation is independently consistent: `rho_z=-rho^2 a` and the scalar
ODE give `alpha_z/alpha=rho(5a/2-h)=-C/2`, with `C=rho(2h-5a)`. Differentiating
C gives

```
C_z=rho^2(2h^2-10ah+5a^2-2U2-5p^2/3)-5.
```

The constant -5 cancels the temporal constant after inserting C into the
conjugated wave equation, leaving
`W=rho^2 Qpot-(35/4)tanh(t)^2`. On the full **new** rectangle
`-ell<=t<=0,0<=z<=8ell`, the all-real inequality `|tanh(t)|<=|t|` consequently
gives `16000<W<92000`. In all response and
node estimates below we retain the primary safe cap
`mu=ell^2*103000=103/10^9`. The potential is positive; “negative-time bound”
must not be described as a negative-valued potential.

## The two principal L2 channels

Put `s=-x` and `A(s,r)=U(-s,r)`. Its initial velocity is
`-ell exp(-ell*r)Q(r)`, with the minus sign essential. Finite propagation gives
support `[1-s,2+s]`, so for every fixed `tau<1` the whole support cone lies in
`(0,3)`. A bounded smooth coefficient extension used for a whole-line formula
has no physical content: its values outside this positive-depth cone cannot
affect the solution. Nothing about the shell interaction at or before x=-1 is
inferred.

Translations and the wave kernel give, in either L2 or supremum norm,

```
B_tau=(1+ell*tau)/(1-mu*tau^2/2),
sup_s ||A(s)||_2 <= B_tau ||Q||_2,
sup |A| <= B_tau.
```

The denominator comes from the actual Volterra operator norm
`mu*tau^2/2`, not from an assumed past evolution error. Smooth data and
coefficient regularity supply the usual smooth linear solution. The compact
endpoint Cauchy estimate gives `||Q||_2<=||Q'||_2/2`. Transporting
`Jplus=U_x+U_r` and `Jminus=U_x-U_r` along their characteristics needs no
potential derivative. With `D=||Q'||_2>0` and `B=B_1` used only as a limiting
constant, it yields

```
(1-delta)D <= ||Jplus/minus(x)||_2 <= (1+delta)D,
delta=ell/2+mu*B/2=2206000103/3999999794000000 < 10^-6.
```

Each principal quadratic-weight share is thus bounded below by
`(1-delta)^2/((1-delta)^2+(1+delta)^2)>499999/1000000`, and above by its
complement. The optimized denominator is necessary; replacing both denominator
terms by the larger bound would lose the desired margin.

These are the normalized-master weights. A nonconstant positive conjugation
preserves neither their values nor purity conditions. No physical energy,
Noether norm or radiation flux is identified, and this finite preshell result
is not an all-past incoming-state obstruction.

## The past corner and the two face zeros

At `(-1/2,3/2)`, both free displacement endpoint values vanish and

```
U_free=-(ell/2) integral_1^2 exp(-ell*r)Q(r) dr.
```

On `[5/4,7/4]`, `Q>=exp(-4/3)>1/4`; a rational exponential-series upper bound
proves `exp(4/3)<4` without a bump evaluation. Also `exp(-ell*r)>=1-2ell`.
With `B_half=(1+ell/2)/(1-mu/8)`, the common corner is at most

```
-ell(1-2ell)/16+(mu/8)B_half
=-3175991536500103/63999999176000000000000 < -4e-8.
```

At each face midpoint `(-1/4,5/4)` and `(-1/4,7/4)`, the free value is at least
`1/2-ell/4`. A local quarter-time error, stronger than the primary's uniform
half-slab error, gives

```
U >= 1/2-ell/4-(mu/32)B_quarter
  =31999983691/63999999794 > 0.499999.
```

Continuity proves an interior zero on **each** face between the common corner
and its midpoint. It proves neither uniqueness nor simplicity of those zeros.
The positive conjugate longitudinal master has the same trace signs and zeros;
the reconstructed scalar itself contains a radial derivative and does not
inherit that assertion automatically.

## What the characteristic construction establishes

Use `a=r+x-1,b=2-r+x`. The half-diamond is exactly
`a,b>=0,a+b<=1`, and `U_ab=-(Vscaled/4)U`. Every rectangle `[0,a]x[0,b]`
stays in that triangle; `ab<=1/4` gives a genuine characteristic Volterra
contraction of norm at most `mu/16<1`. No future x>0 rectangle is needed.
The two actual traces share the corner value. Their smooth scalar Goursat
solution is unique and reproduces the Cauchy development by restriction and
uniqueness. They are a constructed description of that selected linear state,
not independently prescribed physical preparation data. Their negative common
corner also prevents describing them as two independent nonnegative compact
pulses or extending both by zero across the corner.

For `u=x+r,v=x-r`, the trace derivatives satisfy `2f_v=Jminus` on u=1
and `2g_u=Jplus` on v=-2. Derivatives with respect to x along the faces
equal the corresponding J quantities.

This scalar interior diamond is not the earlier shell-corner null cut. The
result does not solve all Einstein null constraints, scalar/Israel corner data,
a nonlinear Goursat problem or a second constrained preparation. Retained
linear reconstruction and full harmonic one-jet/uniqueness premises are needed
to identify the coupled first perturbation. A momentum-reflected exact-family
existence interpretation is conditional; the contact-free first-master support
does not bound the support or development time of higher-order constraint
corrections. No physical source-selection law, uniform nonlinear amplitude
window or cosmological-origin claim follows.
