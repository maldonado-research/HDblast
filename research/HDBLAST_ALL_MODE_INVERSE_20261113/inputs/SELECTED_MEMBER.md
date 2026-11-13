# Quantified high-order second initial data

6 November 2026 UTC. This conditional mathematical checkpoint gives
evaluated high-order bounds for the actual second initial data of the
explicitly specified weighted Einstein–scalar family. It retains the
background, first packet, affine scalar/density prescription and selected
correction range. The historical weight shape remains unspecified; this
numerical member is not identified with it.

The full geometric background-covariant Sobolev bounds on one bulk copy are

```
||gamma2||_H14 <10^2131,     ||tau2||_H13 <10^355,
tau2=pi2/rho0,
I2=||gamma2||_H14+||K2||_H13+||phi2||_H14+||Pi2||_H13 <10^2132.
```

Here rho0 is the fixed background volume density, K is physical extrinsic
curvature and Pi is normal scalar momentum. Norms sum the squared
covariant derivatives using gamma0 and its full geometric volume. The
second-amplitude derivatives are evaluated at zero; the e² Taylor
coefficients are half of them. The constants are extremely conservative
upper budgets in model normalization, not estimates of a measured effect.
They have no SI calibration. A reflected two-copy norm has its separate
sqrt(2) factor, covered by the rounded I2 ceiling 10^2133.

The specified member has

```
ell=10^-6, A={ell/2<z<5ell/2}xS3,
r=(z-ell/2)(5ell/2-z)/ell², kappa=ell/100,
d=kappa r, beta=kappa, a_w=2, qweight=d², w=exp(-1/r).
```

The construction remains conditional on the retained analytic background
and weighted inverse/implicit-family theorem. The new work quantifies the
high-order estimate for its actual second source. It does not establish a
uniform nonzero amplitude interval, an evolved norm or the nonlinear
proper-time, actual-area quadrupole remainder.

The main steps are:

- Global source derivatives from exact coupled-ODE Taylor majorants and
  Bell-polynomial bump bounds, including the canonical density convention.
- Constructive weighted coercivity for the scalar and momentum potentials.
- The complete coordinate ds=dy/qweight and unitary factor
  w*sqrt(qweight*rho³), retaining the radial volume Jacobian.
- Global distribution-to-Sobolev estimates with explicit Fourier constants.
- Covariant tensor reconstruction, angular/connection estimates, zero
  extension, and the exact conversion to physical K2 and Pi2.

Read the [physical Cauchy-jet theorem](PHYSICAL_CAUCHY_SECOND_JET.md),
[source proof](SOURCE_JETS.md), [scalar coercivity](SCALAR_COERCIVITY.md),
[scalar high-order estimate](SCALAR_REGULARITY.md),
[momentum coercivity](MOMENTUM_COERCIVITY.md),
[momentum regularity](MOMENTUM_REGULARITY.md) and
[momentum extension](MOMENTUM_EXTENSION.md).
[Source pins](SOURCE_PINS.json), [reproduction](REPRODUCE.md),
[internal review](INTERNAL_REVIEW.md), [primary references](REFERENCES.md)
and [continuation](CONTINUATION.md) delimit the result and its next steps.

No physical Big Bang origin, energy or force, observational confirmation,
external peer review or journal acceptance follows from these bounds.
This GitHub working checkpoint creates no new Zenodo edition or journal
submission. Earlier frozen papers, capsules and releases remain unchanged.
