# Actual canonical seed: a numerical graded-field envelope

This calculation concerns the registered first canonical datum, not an
evolved solution or the completed wave chart. Put `b=1/rho`, `theta=dy`,
`P=gamma0-theta tensor theta`, and let `n_A` be the four ambient coordinate
functions restricted to the unit three-sphere. With `Y=n_1^2-n_2^2`, the
background-orthonormal, fixed-density seed fields are

```
h = 2 u Y gamma0 - 6 u Y theta tensor theta,
psi = v Y,
tau = -(3/2) exp(-z) u b Y theta tensor theta,
omega = v_t b Y.
```

Here `tau=pi_1/rho0` with both indices lowered by the parallel background
metric, and `omega=varpi_1/rho0`; the sign of the displayed tau is immaterial
for its norm. These are the actual nonzero-velocity canonical data pinned
in the preceding source capsule.

## 1. Radial derivatives, with no sampling

The companion `checks/seed_envelope.py` imports the pinned rational
majorant producer `HDBLAST_WEIGHTED_REGULARITY_20261106/derive_source_jets.py`.
It recomputes its polynomial ODE induction, global flat-bump Bell bounds,
and their composition, all with outward-rounded nonnegative rational
Taylor coefficients. Its stored order is sixteen. Differentiating `v`
once consumes one order, so all orders used below are available. It proves

| Profile | Derivatives | Common absolute upper bound |
|---|---:|---:|
| `u` | 0 through 14 | `10^86` |
| `v` | 0 through 14 | `10^97` |
| `(3/2) exp(-z) u b` | 0 through 13 | `10^77` |
| `v_t b` | 0 through 13 | `10^88` |
| `a`, `b` | 0 through 14 | `1` |

The recurrence is uniform on the entire closed seed support, including
its flat endpoints. The program uses `need`, not optimizable assertions.
It checks the exact polynomial coefficients and constants of its producer
before evaluating the bounds. Its import is an explicitly pinned dependency,
not an independent reimplementation of that producer.

## 2. Every covariant tensor slot

Write `v_A=d n_A` on the warped annulus. The exact background identities are

```
nabla theta = a (gamma0 - theta tensor theta),
nabla n_A = v_A,
nabla v_A = -b^2 n_A (gamma0 - theta tensor theta)
            -a (theta tensor v_A + v_A tensor theta),
nabla gamma0 = 0,
nabla f(y) = f_y theta.
```

The first term in the Hessian identity is the sphere Hessian; the two mixed
terms are the warped connection. Thus this calculation differentiates the
original tensor slots and all newly introduced derivative slots. It does
not estimate angular components by a finite harmonic truncation.

In the background Hilbert tensor norm, `|theta|=1`, `|n_A|<=1`,
`|v_A|<=b<1`, and `|gamma0|=2`. Expand the four displayed seed formulas using
the two monomials of `Y`. Each field has at most four initial monomials,
numerical coefficients at most six, and at most five factors when the
radial profile is counted. Each differentiation of a factor gives at most
four terms and increases its factor count by at most four. Treat each
radial derivative of `a` or `b` as a factor, including the two factors in
`b^2`. The product rule then bounds the number of terms at derivative order
`j<=14` by

```
4 product_{i=0}^{j-1} [4(5+4i)].
```

There are at most `j+1` parallel-metric factors: at most one initially, and
each derivative introduces at most one. All radial factors other than the
one seed profile have absolute value at most one by part 1. Consequently
all full covariant jets needed here have pointwise norm at most

```
4 * 6 * 2^15 * product_{i=0}^{13} [4(5+4i)] * 10^97 < 10^140.
```

This bound is deliberately loose, but includes mixed radial/angular
derivatives and every tensor component.

## 3. Graded norms

For `d=0,1,2` use the literal derivative-sum norm

```
E_d(T) = sum_{j=0}^{14-d}
         || w^-1 q^(j+d-2) nabla0^j T ||_{L2(gamma0)}.
```

On the actual seed support `z in [ell,2ell]`, `r>=3/4`,
`q=10^-16 r^2>10^-17`, `q<=1`, and `w^-1=exp(1/r)<4`.
The latter inequality has the same exact rational Taylor upper bound as
the pinned source producer. Also

```
volume(K) <= 2 pi^2 * 131^4 * 10^-6 < 100^2.
```

For the final strict inequality one may use `pi<4`. Thus the sum of the
four seed norms is bounded by

```
E_0(h)+E_0(psi)+E_1(tau)+E_1(omega)
 <= (2*15*4*10^34 + 2*14*4*10^17) * 100 * 10^140
 < 10^180 < B_seed := 10^200.
```

Their zero extensions are smooth, because the bump and every derivative
are flat at the two endpoints. The same fields lie in the corresponding
completed graded spaces. This is an envelope for this seed and this
distinguished weight member only. Its enormous slack is not an estimate
of a physically relevant perturbation amplitude.
