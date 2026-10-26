# Independent audit of the candidate first transport variations

Equations (4)–(7) in the adjacent geometric-observer working note are correct
with the stated signed `p=partial_y phi0`, nonzero-gradient collar, background
frame conventions, and corrected boundary identification.

Choose the scalar-normal orientation as
`n_phi=sign(p) grad(phi)/sqrt(g^-1(dphi,dphi))`, so its background vector is
`ey`. Differentiating the inverse metric, scalar gradient, and normalization
cancels the radial `s_y/p` contributions and gives

```
(n_phi)1=-(A Y/2)ey+(B-s_t/(rho p))Y e0+(s/(pR)-C)D^aY ea.
```

The transport equation `n_phi(tau)=0`, with `tau0=t`, then yields

```
partial_y tau1=(-B/rho+s_t/(rho^2 p))Y.
```

Direct differentiation of the normalized clock gradient gives

```
(U_geo)1=(D Y/2)e0-[B Y+rho partial_y tau1]ey
          -[E D^aY+(rho/R)D^a tau1]ea.
```

Subtracting the registered coordinate-normal pair therefore confirms

```
b=(B-s_t/(rho p))Y=-rho partial_y tau1,
v_a=-(rho/R)D_a tau1,
w_a=s D_aY/(pR).
```

Both boundary clarifications identified in this audit have been applied to
the main note:

- For a perturbed shell embedding with time displacement `zeta^t`, the effective
  Eulerian boundary datum is `beta=beta_pullback-zeta^t`. Thus the integral
  formula has boundary value `tau1(y_s)=beta`. A purely radial displacement
  contributes nothing at this order because `partial_y tau0=0`.
- Scalar-normal transport yields the required observer at the background. Its
  observer at the perturbed shell is not asserted to equal the intrinsic shell
  geodesic congruence without compatibility between the scalar normal and the
  shell normal.

The fixed positive calibration `rho_s` introduces no first variation. On a
connected collar where signed `p` stays bounded away from zero, its orientation
choice is smooth. The note still requires the specified clock/section data,
valid collar, and transport certification before numerical use.

This is an algebraic audit of candidate first variations. No code, actual-field
evaluation, new curvature coefficient, or sign claim is added.
