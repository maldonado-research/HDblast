# Independent scalar observer-tilt reduction

The proposed compact correction is correct for the stated same-event scalar
observer tilt. This note is separate from the frozen metric/coefficient work.
It evaluates no actual-field intervals and establishes no new numerical sign.
It also records an exact cancellation after substitution of the registered
first-metric form; that cancellation uses no field equation.

Use an oriented adapted orthonormal background frame `(U,n,e_a)` and
`v_a=D_aY`, `h_ab=D_aD_bY`,
`A_ab=epsilon_(a^{cd}v_c h_{b)d}=2(Tplus-Tminus)`. Take the relative first frame
variation

```
delta U1=b0 Y n+V v^a e_a,
delta n1=b0 Y U+W v^a e_a,
delta e_a1=V v_a U-W v_a n.
```

Here `b0,V,W` are coefficient functions of the nonangular variables. Since
`C0=0`, direct differentiation of the four Weyl slots gives

```
Delta H''_ab=Sym_ab epsilon_a^{cd}
 [b0 Y C1_bncd+V v^e C1_becd+2V v_c C1_bUUd
  -W v_b C1_nUcd-2W v_c C1_bUnd].                         (1)
```

The symmetrization is normalized. Each paired last-slot contribution uses
last-pair antisymmetry, giving the displayed factor two and its sign.

## Angular algebra

For the complete first scalar metric sector, the mixed blocks have the forms

```
C1_bncd=z(kappa_bc v_d-kappa_bd v_c),    C1_nUcd=0.
```

The symmetric epsilon contraction therefore eliminates the `b0` term and the
first `W` term. The first identity follows from scalar Codazzi curls and the
constant-curvature Hessian commutator; the second follows because a scalar
first sector supplies no angular antisymmetric two-form.

The all-angular algebraic curvature tensor on a three-dimensional screen can
be written uniquely as

```
C1_becd=kappa_bc S_ed-kappa_bd S_ec-kappa_ec S_bd+kappa_ed S_bc.
```

Consequently, for any vector `v` and symmetric `S`,

```
Sym_ab epsilon_a^{cd} v^e C1_becd
 =-2 Sym_ab epsilon_a^{cd}v_c S_bd.
```

Write the Hessian coefficients of `S`, `C1_bUUd`, and `C1_bUnd` as
`S_H`, `E_H`, and `X_H`; their trace terms have zero symmetric dual. Equation
(1) reduces to

```
P=E_H-S_H,       Q=X_H,
Delta H''=2(VP-WQ)A=4(VP-WQ)(Tplus-Tminus).                 (2)
```

This establishes the radial-boost cancellation and all signs and factors
without assuming the desired final cancellation.

## Independent local ADM and Gauss derivation of P and Q

Put `R=rho cosh(t)`, `theta=tanh(t)`, `H=theta/rho`, `a=rho_y/rho`, and use
ordinary physical derivatives. The angular Hessian coefficients of the first
*orthonormal* time and radial second fundamental forms are

```
m=M_t/(2rho)-E/R,
j=M_y/2-C/R.
```

The first mixed time extrinsic curvature has
`K(n,e_a)_1=(eta/R)v_a`, where

```
eta=R(C_t-theta C)/(2rho)-R E_y/2-B/2.
```

It is essential here to use orthonormal-frame curvature variations, including
the frame changes acting on nonzero background Riemann. First Weyl itself is
frame independent because its background vanishes. With these conventions,
the angular first Riemann restricted-Schouten Hessian coefficient is

```
S_R,H=-L/(2R^2)-a j+H m.
```

The first term is the intrinsic three-screen curvature variation. The `M`
Hessian metric contribution is an angular diffeomorphism of a constant-curvature
screen and cancels after orthonormal-frame variation. The last two terms are
the radial and time second-fundamental-form contributions in the screen Gauss
equation.

The time-normal evolution identity gives

```
(R1_bUUd)_H=m_t/rho+2H m+D/(2R^2)-a j.
```

The lapse term uses `N1=-rho D Y/2`: its spatial Hessian contributes
`-D/(2R^2)+a j` before the minus sign in the Riemann evolution identity.
The mixed Codazzi identity gives

```
(R1_bUnd)_H=m_y-eta/R^2+a m+H j.
```

Subtracting the Weyl Schouten terms shifts the electric and restricted
all-angular Hessian coefficients by the same amount, so they cancel in
`E_H-S_H`. The mixed `bUnd` Schouten term is pure screen trace and does not
change its Hessian coefficient. Thus

```
P=m_t/rho+H m+(D+L)/(2R^2)
 =(M_tt+theta M_t)/(2rho^2)-E_t/(rho R)+(D+L)/(2R^2),

Q=m_y-eta/R^2+a m+H j
 =(M_yt+theta M_y)/(2rho)-E_y/(2R)
  -(C_t+theta C)/(2rho R)+B/(2R^2).                         (3)
```

All `a` terms and the intermediate `theta E` terms cancel. The sign of the
mixed normal derivative agrees with
`R_bUnd=D_n K_bd-D_d K_bn` for the stipulated positive time extrinsic curvature.

## Exact registered-amplitude substitution

The registered first metric is of the form

```
B=f_t/rho-rho b_y,       C=f/R+R d_y,
D=-2(u+a f+b_t),         E=cosh(t)d_t-b/cosh(t),
L=2(u+a f+theta b),      M=2d.
```

Here `u,f,b,d` may be arbitrary smooth functions for this algebraic check. In
particular the scalar amplitude `u` and the radial gauge combination `a f`
cancel from `D+L`. Differentiation gives

```
D+L=-2b_t+2theta b,
E_t=cosh(t)(d_tt+theta d_t)+(-b_t+theta b)/cosh(t),
E_y=cosh(t)d_yt-b_y/cosh(t),
C_t+theta C=f_t/R+2theta R d_y+R d_yt.
```

Substituting these identities in (3) yields

```
P=0,    Q=0
```

as formal rational identities, without equations of motion or numerical
solution values. Therefore the coefficient (2) vanishes for this registered
first-metric form and the stated scalar relative tilt. This is a stronger
specialization than the generic seven-amplitude scalar metric: generic
independent amplitudes need not satisfy either identity.

The statement is limited to the same-event first observer/radial tilt and the
registered metric representation used above. This audit does not certify a
new geometric clock transport, point/label anchoring, or a universal observer
independence theorem. Its tensor-algebra step also works for arbitrary vectors,
but any extension beyond the stated scalar-tilt comparison should be separately
specified and reviewed rather than inferred from the control report's scope.

## Exact checks

`check_scalar_tilt.py` imports no frozen project code. It verifies the
all-angular contraction on all 18 independent symmetric-tensor/vector basis
pairs, the scalar radial-boost and trace cancellations, and the substitutions
`P=Q=0` in a sparse rational Laurent polynomial ring with independent gauge
jets. Normal and optimized results agree. These checks support the algebra;
the ADM component derivation above remains an explicit mathematical audit.
