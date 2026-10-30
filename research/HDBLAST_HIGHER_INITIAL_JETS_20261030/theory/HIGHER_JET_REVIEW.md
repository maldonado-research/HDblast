# Independent higher-initial-jet and coefficient review

This review concerns sufficient inputs for the differentiated error gate.
It does not claim that the complete differentiated bulk source or a spacetime
C2 error certificate has been supplied.

## Initial derivative budget and background recurrence

The registered scalar displacement differentiates the quadratic fs terms
once in the physical radial direction. Since fs already contains the second
packet derivative, the scalar displacement contains packet derivatives
through three. Its second scaled radial derivative therefore requires input
jets through order five. The old order-four paired implementation does not
contain that second scalar derivative; nominal magnitude caps are not a
replacement for actual-minus-nominal higher-jet errors.

The higher initial helper correctly extends the paired input order to five.
The existing compact-packet input already contains its exact derivative
polynomials through six. The new independent endpoint checker verifies their
parity recurrence and continuous polynomial envelope through six. Physical-z
normalized packet coefficients have the factor `ell^(-n)/n!`; a displacement
c2 coefficient is converted to a scaled radial second derivative by
`2 ell^2`, and a velocity c1 coefficient by `ell`.

For the background, if the existing uniform box bounds the actual normalized
fifth derivative, Taylor's theorem applied separately to each derivative
gives the remainder coefficient

```
choose(5,n) c5 z^(5-n),   n=0,...,5.
```

At n=5 this is simply c5. This does not differentiate an unknown remainder
as if it were a polynomial with a fixed coefficient. The helper preserves
that distinction. Its three paired difference channels retain the original
telescoping multiplication and inverse identities. The generalized
exponential/rational derivative identity through five also retains its
factorial and rising-factor conventions.

An independent endpoint route computes an autonomous ODE jet around the
actual state value at R, rather than using the paired Taylor derivative
remainder. All 23 actual-curvature intervals from the two routes overlap.
The scalar-mean affine-endpoint error sign and magnitude agree.

## Nominal second-derivative jumps and mollifier bounds

The cutoff polynomial matches the packet and its first four derivatives,
but not its fifth derivative, at both joins. Registered nonlinear products
can make particular output jumps vanish; jumps should be computed from the
one-sided output jets rather than assumed from input regularity alone.
In general, the scalar nominal displacement is guaranteed C1 with bounded
piecewise second and third derivatives. Its second derivative need not be
globally Lipschitz.

For an even positive kernel, a second-derivative jump contributes at most
half its magnitude when a convolution crosses it from a base point on one
side. If the jump locations are separated by more than twice the kernel
width, at most one is crossed from a fixed base point. Consequently a valid
componentwise uniform bound is

```
|k_delta*f''-f''| <= (35 delta/128) M3
                       +(1/2) max_join |[f'']|.
```

The old affine continuation contributes the jump `-f''(R-)`. A quadratic
continuation matching the curvature removes that endpoint jump, while any
interior cutoff jumps still have to be charged. The half-kernel endpoint
identity and the actual lower bound for the old continuation are proved in
`ENDPOINT_CURVATURE_NOTE.md`.

The quadratic continuation has exterior second derivative equal to its
fixed interior endpoint curvature kappa, and exterior third derivative
zero. Independently, all 23 endpoint kappa interval magnitudes lie below the
existing `nominal_second_derivative_component_caps[0]` entries. Therefore
the quadratic continuation is C1,1 with the same componentwise global M2
bound used by the existing C1 mollifier lemmas. Positivity, evenness and the
kernel moments give the same inequalities

```
|k_delta*f-f| <= M2 delta^2/18,
|(k_delta*f)'-f'| <= 35 M2 delta/128
```

for either extension, on the retained footprint. This supports a separate
re-derivation of C1 bounds for the new candidate with the same smoothing
constants; it does not silently identify the two candidates or transplant a
frozen certificate without checking its remaining inputs and equations.

## Smoothed P1 velocity derivatives

The nominal velocity is C1 with bounded piecewise second derivative. Thus
its first derivative is globally Lipschitz after tangent continuation.
For a P1 cell of width h, the slope error relative to the nominal derivative
is bounded by

```
Dv_cell=(h/2) M2v_cell +(radius_left+radius_right)/h.
```

The first term follows by integrating the pointwise Lipschitz derivative
difference across the cell; the second charges the two fixed nodal center
errors. The last-slope extension of P1v and the tangent extension of the
nominal velocity preserve the last-cell bound beyond R. Convolution then
gives the safe derivative comparison

```
|(k_delta*P1v)'-nominal_v'| <= max_cell Dv_cell
                                +(35 delta/128) M2v.
```

Actual-minus-nominal velocity-gradient errors remain a separate paired
input. No derivative is inferred by differentiating a velocity magnitude
or interpolation error cap.

## Actual coefficient time derivatives

I independently checked the complete time-dependent primitive expansion.
With `s=sech(ell x)` and `tau=tanh(ell x)`, the time-varying part is

```
D=ell^2 [s^2 Delta +2 rho a s Ca +2 s tau Ch
         +rho a tau(Sah-3O)+tau^2 Sh
         +3 rho p tau W_L+rho p s Ephi+time-constant terms].
```

Using `tau_x=ell s^2` and `s_x=-ell s tau` reproduces all eight coefficient
functions in `check_time_coefficients.py`. The radial background is held
fixed under this time derivative. Also
`A_x=-3ell^2 s^2 I`, `B_x=0`, and the actual shell M is time constant.
The positive hyperbolic interval ranges, exact orthogonal degree blocks,
weighted Frobenius bounds and zero-wedge commutator conversion are
consistent. These derivatives do not supply the missing complete F_x.

The resulting data and coefficient inputs remain conditional model
calculations. They are not a full Einstein residual, a nonlinear realization,
an invariant observable or evidence for a higher-dimensional Big Bang.
