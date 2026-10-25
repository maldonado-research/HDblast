# Fixed finite-time interpretation of the shell volume jet

25 October 2026. A conditional kinematical statement, separate from the
Einstein–scalar–shell existence problem.

The observable is the volume of spatial sections reached at a specified
intrinsic proper time by normal geodesics from the registered physical
initial shell section. Its definition and second-variation conventions are
those of the frozen `HDBLAST_TOTAL_RESPONSE_PREFLIGHT_20261023/observable/`
specification. Parameter coefficients are derivatives: the metric expansion
is `gamma_e = gamma_0 + e m + e^2 n/2 + O(e^3)`.

## Finite tube lemma

Let the background shell be a product neighborhood of
`[-b,T+b] x S^3`, where `b>0`, with metric

```
gamma_0 = -ds^2 + A(s)^2 kappa,
A(s) = rho_b cosh(s/rho_b),
```

and anchor `{0} x S^3`. Assume a smooth family of smooth Lorentzian
metrics `gamma_e` and spacelike anchor embeddings `j_e`, with a common
future orientation, on that neighborhood. Suppose they equal the
registered background data at `e=0`. Then, after shrinking the parameter
interval, the future normal exponential maps are embeddings on
`[0,T] x S^3`. Their constant-proper-time sections are spacelike.
The parameter interval can depend on the chosen finite `T`.

Here is a proof that includes global injectivity. Background unit normal
geodesics are `(s,z)=(tau,z)`. They exist beyond the compact target tube,
and their map is an embedding there. Smooth dependence of the geodesic
ODE on its metric, initial point, and normal vector gives existence on
the whole compact target interval and convergence of the perturbed maps
to the background map in `C^1`, for sufficiently small parameter. The
buffer prevents loss of existence at the endpoint. In particular their
differentials remain invertible uniformly on a slightly larger tube.

If injectivity failed arbitrarily near zero, compactness would give a
sequence of distinct domain pairs with equal images, converging to two
background domain points. Injectivity of the background map makes those
limit points equal. The inverse function theorem applied to
`(e,tau,z) -> (e,F_e(tau,z))` in a neighborhood of that common point then
rules out the collisions. The time endpoints cause no problem because
of the larger open tube. A continuous injection from the compact target
tube into the Hausdorff shell is an embedding. Gauss' lemma gives
`F_e^* gamma_e = -d tau^2 + h_e(tau)`; the spatial metric is positive
by compact `C^0` stability of the positive background spatial metric.

Consequently the specified target
`T = rho_b * ell * 9/8`, `ell=1e-6`, is admissible for a smooth
kinematical realization defined on an open buffer inside the retained
history `x<5/4`. This is an existence statement for some sufficiently
small, unspecified parameter interval. It supplies no numerical
finite-amplitude range, curvature bound, or measured caustic radius.

Smoothness is a sufficient hypothesis. This lemma does not infer a
smooth realization from a merely `C^1` numerical candidate. The
observable jet below assumes smooth exact mean histories in the parent
coefficient problem; the numerical candidate approximates those means.
It does not assert that a full angular Einstein solution has been found.

## What a kinematical realization means

Given prescribed smooth induced metric jets `m,n` and the smooth anchor,
`gamma_e=gamma_0+e m+e^2 n/2` supplies such a metric family on a fixed
compact tube for sufficiently small `e`. The Lorentzian signature and
spacelike anchor persist by openness and compactness. This construction
does not solve nonlinear Einstein equations, scalar equations, or shell
junction conditions. An exact nonlinear physical family with these jets
remains a separate problem. The second volume variation depends only
on these jets, so different smooth kinematical realizations with the same
jets give the same second coefficient.

For this functional, a full angular second-order physical solution is
unnecessary. Given the exact smooth mean histories, one may choose the
explicit kinematical completion
`n = N0(s) ds^2 + J0(s) A(s)^2 kappa`, with zero second-order shift.
The known full first-order tensor `m` is retained. The volume second jet
depends linearly on `n` only through `mean(n_ss)=N0` and
`mean(tr_(h0) n_ij)=3J0`. Adding a second-order tensor with both these
means zero leaves the coefficient unchanged; its second-order shift
does not enter this variation. The completion need not satisfy the
missing second-order physical constraints. Identifying its coefficient
with a physical spacetime response remains conditional on the parent
mean projection and subsidiary assumptions.

## Dependence on the mean coefficients

For the registered quadrupole and anchored proper-time clock, the frozen
determinant and clock calculation gives

```
v2 = V2/V0 = (3/2) [J0(s) + H(s) integral_0^s N0(a) da] + Q_first(s),
H(s) = tanh(s/rho_b)/rho_b.
```

`Q_first` contains the complete quadratic first-order clock, shift and
spatial-metric terms in the parent specification. The first volume
coefficient vanishes. Thus `R2/R0=v2/3`; the second expansion coefficient
also involves derivatives of `J0` and the clock term. A sign of `N0`
alone or of a candidate `J0` alone does not establish a sign of this
observable. The total conditional response, all error terms, the
physical anchor and the proper-time correction must be used together.

This lemma is standard stability of a compact normal tube and smooth
ODE dependence. It is supplied to clarify a finite-time use of the
existing jet observable, with no claim of a new general geometric theorem.
