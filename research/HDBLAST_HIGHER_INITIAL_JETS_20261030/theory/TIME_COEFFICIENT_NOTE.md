# Actual time derivatives of all full23 homogeneous coefficients

This closes the homogeneous coefficient-derivative input of the proposed
differentiated-error gate. It does not provide the differentiated complete
bulk source, higher initial errors, or a quantitative C2 error theorem.

## Fixed-coordinate analytic differentiation

Use scaled time `x=t/ell`, with `ell=10^-6`, and fixed inward scaled radius r.
The background variables `rho,a,p,eta` and their potential/curvature functions
depend only on r. Put

```
tau=tanh(ell*x),  s=sech(ell*x),  h=tau/rho,
tau_x=ell*s^2,    s_x=-ell*s*tau, h_x=ell*s^2/rho.
```

The actual homogeneous equation is
`U_xx-U_rr=A U_x+B U_r+D U+F`, with the same actual coefficient assembly as
the registered full23 response. The matrices below are fixed exact rational
angular/fibre operators. The time-dependent parts of the assembly are

```
A=-3ell*tau*I + radial-background-only terms,
B=radial-background-only terms,
D=ell^2 [s^2 Delta+2rho*s*a Ca+2rho*s*h Ch
          +rho^2 C+rho*p*s Ephi],
C=C_static+a*h(Sah-3O)+h^2 Sh+3p*h W_L.
```

All potential and curvature terms belong to `C_static`. No potential or
curvature term has been omitted: its x derivative is zero because r is held
fixed. Product differentiation therefore gives

```
A_x=-3ell^2 s^2 I,      B_x=0,      M_x=0,
D_x=ell^3 [ -2s^2 tau Delta -2rho*a*s*tau Ca
            +2s(s^2-tau^2) Ch
            +rho*a*s^2(Sah-3O)+2tau*s^2 Sh
            +3rho*p*s^2 W_L-rho*p*s*tau Ephi ].
```

For example, the `Ch` coefficient differentiates both factors of `s*h`:
`2rho(s_x h+s h_x)=2ell*s(s^2-tau^2)`. Omitting either product term gives
an incorrect result. The distinction between partial_x and partial_t supplies
the third power of ell in the D_x formula.

The actual Robin matrix M uses only the shell background and scalar-tension
derivatives. It is time independent in this frame, so its time derivative is
exactly zero; no Robin matrix inverse or spectral sign assumption is needed.

## Continuous interval enclosure and comparison norm

The computation uses the containing rectangle
`0<=x<=3/2`, `0<=r<=7/4`, which covers the retained cone D+. Actual radial
background intervals come from the authenticated response Context. For
nonnegative t in this rectangle, the analytic bounds

```
t-t^3/3 <= tanh(t) <= t,
1-t^2/2 <= sech(t) <= 1
```

give valid continuous hyperbolic ranges. The identities are differentiated
before interval substitution. A derivative of an interval endpoint is not
used. All arithmetic is rational with outward interval rounding. Displayed
`A_x` matrix endpoints are kept as exact Fractions so the sharp lower
endpoint `-3ell^2` is not unnecessarily widened by dyadic rounding.

The full23 positive H Gram retains the scalar, symmetric mixed-component,
gradient, trace-free Hessian and both TT weights. Exact primitive block
closure gives orthogonal blocks of dimensions `5,8,8,1,1`. Each D_x block's
weighted Frobenius cap bounds its H operator norm, and the maximum bounds
the complete operator. No mean-only or TT-only replacement is made.

The resulting whole-rectangle bounds are

```
||A_x||_H <= 3e-12,
B_x=M_x=0 exactly,
||D_x||_H < 2.354e-16.
```

The sharper recorded D_x cap is approximately `2.35313021018e-16`. Its five
block caps are approximately
`[2.005116306e-16,2.339491779e-16,2.353130210e-16,
  2.4000000001e-23,2.4000000001e-23]`.
The A_x bound uses its exact shared scalar-identity form rather than a
dimension-dependent Frobenius bound; at x=0 it attains `3ell^2`.

For an error with inherited zero-wedge value bound
`||e||_H<=(3/4)X(e)`, this supplies the differentiated coefficient term

```
||A_x e_x+B_x e_r+D_x e||_H
 <= [a1+b1+(3/4)d1] X(e),
a1+b1+(3/4)d1 < 3.000177e-12.
```

The numbers refer to the declared scaled-coordinate H norm. They are not
coordinate-invariant physical smallness statements. The complete
differentiated forcing defect and the higher paired initial-data errors
remain separate inputs before the C2 gate can be evaluated.

## Reproduction and independent algebra

`check_time_coefficients.py --repo-root REPOSITORY` authenticates the raw
coefficient driver, primitive operator and profiles before project imports.
The Response Context authenticates all seven scientific upstream roots;
the checker reauthenticates the complete dependencies afterward. It uses
actual background intervals and fixed rational primitives, not the driver's
frozen midpoint coefficient matrices. Normal and optimized outputs match
byte for byte.

The separate `check_time_derivative_identity.py` evaluates the original
primitive assembly with exact rational dual numbers and compares every
entry against the differentiated formula. Its independent product-rule
construction, 17 groups and eight wrong-formula controls are described in
`TIME_DERIVATIVE_AD_NOTE.md`. Manufactured rational fixtures verify algebra;
the displayed continuous caps separately require the authenticated actual
background ranges and interval argument above.

This is conditional homogeneous-operator input preparation. It does not
bound the derivative of the complete actual bulk source, constraint
derivatives, full Einstein/scalar residuals, an exact nonlinear family, or
an observable supporting a higher-dimensional Big Bang origin.
