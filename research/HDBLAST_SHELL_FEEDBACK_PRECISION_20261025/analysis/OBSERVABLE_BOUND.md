# Conditional anchored-volume second-jet bound

The fixed time is `x*=9/8`, with physical proper time
`tau*=rho_b ell x*` and `ell=1e-6`. The clock and anchor are those in
`FINITE_TUBE_OBSERVABLE.md` and the frozen proper-time specification.
This calculation concerns the parameter derivative `V2/V0`, not a
finite-amplitude fractional volume change or an observed expansion rate.

The inherited determinant/clock formula is

```
v2 = (3/2)[J0 + a IN] + Q_first,
a=ell tanh(ell x*),  IN=integral_0^x* N0(x) dx.
```

The complete candidate interval already contains quadrature uncertainty.
It is widened once by the complete actual-response error (input, omitted
source, coefficient defects and the full third-and-higher feedback tail).
At `x*=9/8`, the resulting total `J0` interval is contained in

```
[-11456254991148748, -10003713776284829].
```

The positive ten-state comparison independently gives a direct component
value bound `sup_C(9/8)|N0| < 1.213e22`. The exact zero history through
`x=1` and `|tanh z|<=|z|` give

```
|(3/2) a IN| <= (3/2) ell^2 (9/8)(1/8) Ncap < 2.56e9.
```

There is no extra factor three here: `Ncap` bounds the component value,
rather than the rescaled global X norm.

For completeness the first-order quadratic bound is reproduced explicitly.
The primary retained shell report supplies strict caps

```
X*=2987/2000000000,          Z*=1119/100000000000,
V*=11117/2500000000000000,  d*=33351/10000000000000000000000,
Zt*=28811/1000000,          dt*=11117/5000000000000000.
```

`V*` bounds the observable shift `V_gauge`, not another source coefficient.
Set `delta=1/8`, `C*=delta X*`, `b*=ell^2(9/8)`, `t*=ell(9/8)`,
`S*=V*+ell C*`, `E*=Z*+b*C*`, `Zx*=ell Zt*`, `dx*=ell dt*`.
The full scaled clock/shift/determinant formula and triangle inequality give

```
|Q_first| <=
 b*[delta X*^2/2 + 4 delta S*^2 + C*X*]
 + C*/3 [3Zx* + 6b*Z* + 8dx* + 16b*d*]
 + ell^2 C*^2 (1+t*^2)/2
 + 8ell C*V*/3 + 4ell^2 C*^2/3
 + E*^2/2 + (8/3)E*d* + (16/3)d*^2
 = 6686222274308036835249874398211279
   /1228800000000000000000000000000000000000000000000
 < 5.5e-15.
```

This uses `cosh(ell x)>=1` and the zero first-order history through `x=1`.
The ordinary-time derivative conversion is essential. No sign assumption
or correlation among the first-order quantities is used.

The exact rational join in `join_results.py` therefore encloses the
conditional anchored-volume coefficient by the integer interval

```
v2 = V2/V0 in [-17184385045350733, -15005568105799633] < 0.
```

The first volume coefficient is zero. Accordingly `R2/R0=v2/3` also has
negative sign for the volume radius `R=(V/(2 pi^2))^(1/3)`. No sign of
the volume-averaged expansion coefficient follows without a separate
time-derivative calculation. Negative second derivative does not assert
negative total volume, a contracting universe, or a physical force.

The finite-tube lemma uses a smooth kinematical completion of the exact
mean histories and the known first-order induced metric. It supplies
an unspecified sufficiently small parameter interval, without solving
the missing angular second-order field equations. All regular-branch,
mean-projection and subsidiary premises of the parent remain conditional.
The full physical second-order construction, an exact nonlinear family,
and observational evidence for the Big Bang hypothesis remain open.
