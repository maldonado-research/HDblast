# Independent internal review of the natural residual bridge

30 October 2026. The reviewed equations are the trace-reversed metric equation
`E_ab=Ric_ab-dphi_a dphi_b-(2/3)Vg_ab` and the natural scalar equation
`S=Box phi-Vprime`. They are evaluated at the second **parameter derivative**
of the registered expansion. `E` is not the Einstein tensor `G` itself, and
the lambda-squared coefficient is half the stated parameter Hessian.

## Formula review

For `R=Z_xx-Z_rr-La Z-F`, the metric forcing normalization
`F_metric=-2 ell² rho² QEred` gives the reduced metric residual
`R_metric/(2ell²rho²)`. Adding the natural gauge-gradient and the difference
between the actual natural and reduced quadratic sources therefore gives

```
E2[Z]=R_metric/(2ell²rho²)+nabla_(a DC[Z]_b)+(QE-QEred).
```

The parentheses include one half. The exact complete source-map subtraction
is appropriate under the retained actual `C1=0` neighbourhood hypothesis.
It must not be interpreted as an arbitrary off-shell identity; the general
gauge expansion also has the first-gauge-variation coupling.

The primitive homogeneous scalar assembly retains the natural
`-p DC_y[Z]` term. In particular its `Aphi` and `Bphi` rows include the
time `L_x` and radial `K_r+N_r-3J_r` pieces; their algebraic background
pieces combine with the scalar Hessian coupling into `-(py+4ap)K`.
Consequently its source convention `F_scalar=ell²rho² Qphi` gives

```
S2[Z]=-R_scalar/(ell²rho²).
```

No scalar harmonic term should be deleted or added a second time.

## Physical derivative and norm review

The preceding exact harmonic formula supplies
`P=ell rho DC[Z]=Rx Z_x+Rr Z_r+Zmat Z` in the positive S3 probability L2
norm. The separately bounded derivative inputs give the legitimate rational
factor `15/8>sqrt(7/2)`. Independently summing the degree-wise weighted
Frobenius contributions gives `245/2,131/2,68` for the alpha, eta and zeta
value matrices, so their respective inclusive caps `12,9,9` are valid.

The time and radial value-matrix derivatives in the report follow from the
actual-background identities. In particular `rho_r=-alpha rho`, so

```
e_y DC=-(P_r+alpha P)/(ell²rho²).
```

This radial product term has been retained with the correct sign. The
temporal-frame connection contributes at most `|a| ||DC||`. The three
sphere-direction extrinsic connections obey a universal positive-square
identity:

```
3(a²+h²)|C|²-sum_i|Omega_i C|²
 =3(h C_y+a C_0)²+2(a²+h²)|C_sphere|² >=0.
```

The independent checker expands both sides as exact polynomials, rather
than checking sample covectors. The intrinsic sphere derivative has squared
norm at most 24 on the retained scalar/gradient covector harmonics. Combining
the time, radial and angular derivative-direction bounds by the triangle
inequality gives the stated conservative connection factor
`2 alpha_abs+5ell+(7/4)(alpha_abs+eta_abs)`.

Symmetrizing a covariant derivative cannot increase its positive Frobenius
norm. The metric and scalar output projections likewise have norms bounded
by the full positive full23 H norm. The complete quadratic source correction
uses all 20 metric amplitudes and their full Frobenius weights.

## Numerical audit and scope

`check_natural_bound_audit.py` independently reconstructs the actual
background charges, the degree-wise matrix normalization, all 20 quadratic
source correction amplitudes, the smoothed-source discrepancy charge, and
the resulting exact rational residual caps. The 21 named checks include
four wrong-formula controls. Normal and optimized outputs agree byte for
byte. The reviewed readable bounds are

```
||E2[Z_epsilon]||_positive,physical-frame,L2(S3) <2.918e43,
||S2[Z_epsilon]||_physical-frame,L2(S3)          <2.118e30.
```

These are large absolute upper bounds. They establish finiteness under the
retained premises; they establish no relative smallness or useful physical
accuracy. The bound for the actual target's second-derivative error and the
complete actual `F_x` remain open. No exact constraint propagation, exact
nonlinear family, finite-amplitude remainder or cosmological-origin
conclusion follows from this bridge.

**Assessment:** the displayed trace-reversed metric and natural scalar
coefficient residual bridge and numerical upper bounds pass independent
internal review with the above premises and limitations.
