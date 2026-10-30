# Local causal bounds for the first actual-bulk update

This supplements the global error majorant with a local enclosure at
`(x,r)=(9/8,0)`. The actual bulk convolution is not quadrature-evaluated.
Instead, local initial-profile caps, the physical causal geometry, actual
coefficient bounds and a genuine uniform affine-source bound yield narrower
point enclosures. This local information is absent from the global identity
`D+E1=E0`; that global identity remains valid.

## Fixed functions and the local domain

The reconstruction `W0`, its smoothed free wave `Bs`, and the update

```
W1=Bs+K(La W0+F)+L(M PN W1+G)
```

use the same fixed actual coefficient functions and the same actual affine
sources. No actual source is replaced by a midpoint or zero. Put

```
a=9/8, d=a-1=1/8,
Omega_a={(s,y):s>=0,y>=0,s+y<=a}.
```

This physical triangle is inside the retained cone. Compatible fields and
sources vanish on the inherited zero wedge `s+y<=1`, and boundary loads
vanish through shell time one. These support facts are separate inherited
premises; a magnitude cap alone would not imply them.

The initial data entering this triangle have arguments in `[-a,a]`. The
checker bounds the unsmoothed nominal displacement and its actual radial
derivative using the continuous P1 interpolation-error certificates on all
cells through `a+max(delta_u,delta_v)`. The P1 velocity uses its nodal convexity
bound. Maxima are taken componentwise before the positive Gram norm.
Positivity of the convolution kernel transfers these local caps to the
smoothed functions and derivative. Write the resulting H caps as `U0,U1,V`.
Then safe bounds are

```
||Bs|| <= U0+a V,
||Bs_x||,||Bs_r|| <= U1+V,
||Bs|shell|| <= U0+d V.
```

The last bound uses the zero velocity collar. The all-triangle value bound
is deliberately conservative; a shorter support-based value bound could
also be proved. The numerical calculation retains the displayed formula.

## Physical reflection and active area

Let `f` be a physical-side source with H bound `Fcap`, vanishing when
`s+y<=1`. For an output `(s,y)` in `Omega_a`, the direct and reflected
physical integration domains both lie in `Omega_a`. A physical source point
has multiplicity at most two; the Duhamel prefactor is one half. Reflection
of the **completed coupled output** is an isometry of the diagonal H norm.
Hence the value bound uses the active physical area

```
A=area{(s,y):s>=0,y>=0,1<s+y<=a}
 =(a^2-1)/2=17/128,
||K f|| <= A Fcap.
```

The two characteristic-edge integrals have length at most `a`, so
`||(K f)_x||,||(K f)_r||<=a Fcap`. The smaller constants `d^2/2` for area or
`d` for derivative length are not valid: a characteristic may stay on the
active line `s+y=a` for its entire length. All these bounds hold continuously
throughout the local triangle, not only at its shell endpoint. No matrices
are moved through a parity projector.

## Bounding the earlier reconstruction on the triangle

For the frozen source `f0=L0 Bs`, the local H bound is

```
f0cap=(nA0+nB0)(U1+V)+nD0(U0+a V).
```

Let `m` bound the actual Robin matrix norm. At the shell the trace `Z` of W0
satisfies the actual Volterra equation. Therefore

```
hcap=U0+d V+A f0cap,
Zcap=(hcap+d Gcap)/(1-m d),
gcap=m Zcap+Gcap=(m hcap+Gcap)/(1-m d).
```

The checker verifies `m d<1`. The actual total boundary load is `MZ+G`;
it is distinct from the homogeneous correction load of `W1-W0`. The lift
then gives, everywhere in `Omega_a`,

```
||W0-Bs|| <= A f0cap+d gcap,
||(W0-Bs)_x||,||(W0-Bs)_r|| <= a f0cap+gcap.
```

## A genuine pointwise affine bulk bound

The earlier `SOURCE` bounds `X(KF+LG)` and cannot be inserted as a pointwise
bound on F. The local checker instead authenticates the actual geometric
source report and its complete source map. Its domain `0<=x<=3,r>=0,x+r<=8`
contains `Omega_a`, and its normalization already includes
`metric_forcing=-2 ell^2 rho^2 QEred` and
`scalar_forcing=ell^2 rho^2 Qphi`.

All 93 source block-Gram entries and field multiplicities are compared with
the current full23 H weights. The source norm square is reconstructed from
the recorded component caps and agrees with the positive bundle certificate.
This establishes the uniform bound in the required normalization,

```
||F||_H < 17387513156573.
```

This is a uniform magnitude bound, not a point evaluation or source-sign
determination. The source report's actual common-gauge evolution and complete
geometric projection remain inherited premises.

## Local residual and correction

For `r0=La W0+F-L0 Bs`, write `deltaA,deltaB,deltaD` for actual-minus-frozen
coefficient norm caps. The checker uses

```
Rlocal=(deltaA+deltaB)(U1+V)+deltaD(U0+a V)
       +(nA+nB)(a f0cap+gcap)
       +nD(A f0cap+d gcap)+FH.
```

Actual coefficient norms, error caps, block closure and every actual entry
interval are checked directly from the raw pinned reports. The exact
arithmetic gives approximately

```
Rlocal = 1.905589378957374e25.
```

The correction `C=W1-W0` obeys `C=K r0+L(M PN C)`. Its boundary condition is
homogeneous because the same actual G cancels. A uniform local value/trace
majorant follows from

```
Dtrace=A Rlocal/(1-m d)
      =A Rlocal+d m Dtrace.
```

At the shell, `PD C=PD C_x=0` and `PN C_r=M PN C`. Thus

```
||C_x||_Robin <= a Rlocal+m Dtrace,
|(C_r)_i| <= rowdual_i(M) Dtrace   for Robin rows,
|(C_r)_i| <= a Rlocal/sqrt(H_i)   for Dirichlet rows,
rowdual_i(M)=sqrt(sum_j max|M_ij|^2/H_j).
```

The resulting local majorants are approximately

```
Dtrace = 2.530947563204937e24,
a Rlocal+m Dtrace = 2.143857386748785e25.
```

## Updated target point information

The global one-update error certificate is independently recomputed with the
retained contraction and compared with the separately authenticated bulk
majorant report. The readable global bound is `E=3.541e25`. At the shell,
the error value is at most `E/8`; its time derivative is at most E. Since U
and W1 share the same actual M and G, their radial error satisfies
`PN(U-W1)_r=M PN(U-W1)`. Thus the total target radii around the old W0 sample
are

```
value: (Dtrace+E/8)/sqrt(H_i),
time:  (a Rlocal+m Dtrace+E)/sqrt(H_i),
Robin radial: rowdual_i(M)(Dtrace+E/8),
Dirichlet radial: (a Rlocal+E)/sqrt(H_i).
```

All 69 best target point boxes are contained in the old PR96 boxes, and 53
are strictly narrower. The unchanged cases are homogeneous Dirichlet
value/time boxes and zero-Robin-row source-only radial boxes. These latter
rows are `u:2,u:4,H:2,H:4,T:+,T:-`; zero M does not mean zero source.

For the registered scalar mean, the exact local enclosures give

```
8.29036e30 < U_w:0 < 8.29042e30,
1.974621e32 < U_w:0,x < 1.974623e32,
2.154067e27 < U_w:0,r < 2.154080e27.
```

The widths improve by factors approximately 1.806, 2.684 and 1.806,
respectively, relative to the earlier best target boxes. The improvement
comes from new local analytic information. It does not contradict the
global norm-only width identity and does not evaluate the actual convolution
or determine its correction sign.

## Reproduction and limits

Run `check_local_bulk_enclosure.py --repo-root PATH`. Optional prior-capsule
arguments select separately located authenticated copies. The default
`--bound-report` is `../reference/BULK_BOUND.json`. Normal and optimized Python
produce identical exact JSON. No new primary bulk code or scientific source
producer is imported by this checker.

The result is conditional on the retained first-field, source projection,
background, operator, zero-wedge and formal integral-equation premises. It
does not complete the finite update's bulk self-consistency, provide a
quantitative C2 error or full Einstein/constraint residual, establish a
nonlinear family, or supply a physical observable. The sample quantities are
registered gauge-dependent modal coefficients and coordinate derivatives.
