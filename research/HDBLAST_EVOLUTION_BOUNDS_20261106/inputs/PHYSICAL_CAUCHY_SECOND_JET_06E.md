# From canonical corrections to the physical second Cauchy jet

6 November 2026. This is a conditional initial-data estimate for the
explicit weight member specified in README.md. It retains the fixed model,
first packet, initial identification and canonical affine prescription.
Every subscript 2 denotes the second amplitude derivative at zero, twice
the quadratic Taylor coefficient. All constants are in the fixed model
normalization, without SI calibration.

## Norm and statement

On one bulk copy of the initial hypersurface, use the fixed background
metric gamma0 and its Levi-Civita connection. Define the full geometric
integer Sobolev norm by

```
||T||_Hm(gamma0)^2 = sum_{j=0}^m integral |nabla_gamma0^j T|_gamma0^2 dvol_gamma0.
```

The sphere has its full geometric volume 2pi². This definition has no
uncomputed equivalence with coordinate-component or Fourier norms.
The fields are compared in the retained initial gauge and identification;
these are not invariants under an amplitude-dependent change of that
identification. Corrections extend by zero from the annulus.

Write h2=gamma2, and let K2 and Pi2 be the second derivatives of the
extrinsic curvature and normal scalar momentum. Then

```
I2=||h2||_H14+||K2||_H13+||phi2||_H14+||Pi2||_H13 <10^2132.       (1)
```

For two identical reflected bulk copies with the sum-of-integrals norm,
each norm gains sqrt(2), so the deliberately rounded ceiling is 10^2133.
No double-copy energy interpretation is attached to that bookkeeping.

The premises are the retained background branch and selected weighted
inverse/implicit construction, the actual canonical source, and the
explicitly distinguished weights. Equations (1) and its inputs concern
only derivatives at zero on the initial slice. They do not give a uniform
amplitude interval, a development time or a finite-amplitude remainder.

## Exact canonical conversion

The [canonical constraints](inputs/CANONICAL_CONSTRAINTS.md) use

```
pi^{ij}=sqrt(det gamma)/2 (K^{ij}-(tr_gamma K)gamma^{ij}),
varpi=sqrt(det gamma) Pi,
K_ij=2/sqrt(det gamma) [gamma_ik gamma_jl pi^{kl}
                        -(gamma_kl pi^{kl})gamma_ij/3].
```

At the time-symmetric background K0=Pi0=0. Let rho0 be the fixed
background volume density, not the warp radius rho. The actual first data
in [the binding](inputs/REGISTERED_EXACT_FAMILY_BINDING.md) are

```
h1=uY(-4dy²+2gamma_angular),
pi1/rho0=tY partial_y tensor partial_y, t=-3w_velocity/(2rho),
Pi1=(vt/rho)Y, Y=n1²-n2², gamma_angular=rho²kappa.
```

The canonical second scalar and scalar-density data are phi2=varpi2=0,
whereas pi2=rho0 tau2. Since tr_gamma0 h1=2uY, the first relative volume
derivative is uY. Differentiating the two inverse formulas gives exactly

```
K2=2[tau2_lower-(tr_gamma0 tau2)gamma0/3]
   +(u*w_velocity/rho)Y²(36dy²-6gamma_angular),
Pi2=-2(u*vt/rho)Y², phi2=0.                                  (2)
```

Lowering uses gamma0. No h2 term enters these conversions because the
background momenta vanish. In particular, zero canonical scalar-density
second data do not mean zero physical normal scalar momentum.

One direct check sets tau2=0 and uses the affine metric factors
A=1-4euY and B=1+2euY. The known diagonal terms are
Kyy=(4/3)etY A^(3/2)B^(-3/2) and
Kab=-(2/3)etY rho²kappa_ab A^(1/2)B^(-1/2). Their second derivatives
are -24utY² and 4utY²rho²kappa_ab, yielding 36 and -6 in (2).

The linear tensor map L(tau)=2[tau-(tr tau)gamma0/3] has eigenvalue 2
on trace-free tensors and -2/3 on the trace line in spatial dimension
four. Its operator norm is exactly 2. Since gamma0 is parallel, it
commutes with every covariant derivative, hence

```
||L(tau2)||_H13 <=2||tau2||_H13.                              (3)
```

## Explicit known-product and combined bounds

[SOURCE_JETS.md](SOURCE_JETS.md) proves, over the actual compact packet,
the global proper-distance derivative bounds through order thirteen

```
max_{j<=13}|partial_y^j(u*w_velocity/rho)| <=10^73,
max_{j<=13}|partial_y^j(u*vt/rho)| <=10^84.
```

These use normalized Taylor majorants, not point samples. The flat
packet extension preserves the bounds globally. The radial support
volume is at most 131^4 ell=(17.161)^2. The angular basis Y² and its
needed covariant derivatives admit the same generous 10^200 envelope
as the degree-at-most-four basis in SCALAR_REGULARITY.md.

For completeness, the background y-vector field in (rho^-1,a,p,eta)
has polynomial coefficient l1 norm at most 100 and degree at most 5
on the retained box. Thus its derivatives of a and rho^-1 through
order eighteen have the rational majorant
product_{k<18}100(1+4k)<10^80. In the warped orthonormal splitting,
nabla_gamma0 acts radially as partial_y and angularly as
rho^-1 nabla_kappa plus a times the radial/angular slot-mixing map.
The latter has norm below 100 at the ranks needed here. At each of
thirteen stages the additive branch/norm budget is below 1000,
including derivative slots and directions. There are at most thirteen
background jet factors in each expanded term. This also bounds the
mixed derivatives, rather than just repeated radial derivatives.

The following deliberately oversized number consequently bounds either
known product in (2) in full geometric H13:

```
Bcross=5*200*4*18*14*1000^13*(10^80)^13*10^200*10^84 <10^1600.    (4)
```

The factors allow sphere conversion, all fixed tensor coefficients,
the Sobolev derivative sum, support volume, and the stated connection
expansion. The scalar product has fewer slots and is covered as well.
Compact support and flatness give ordinary zero extension for these
known terms. The selected unknown correction's zero extension is proved
in the regularity/extension notes without assuming a potential trace.

The two high-order estimates are

```
||h2||_H14 <10^2131,
||tau2||_H13 <10^355.                                       (5)
```

The second is a rounded ceiling above the sharper per-copy 10^354.
Combining (2)--(5) gives
I2<10^2131+2*10^355+2*10^1600<10^2132, as claimed.

These very large upper bounds are analytic budgets for initial amplitude
derivatives. They are not estimates of the measured size of a cosmological
effect. Uniform parameter control, coupled evolution, completed physical
observer/shell pullbacks and the nonlinear actual-area quadrupole remain
separate requirements.
