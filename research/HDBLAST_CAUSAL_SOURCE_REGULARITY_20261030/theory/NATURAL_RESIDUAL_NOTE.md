# Direct natural-equation coefficient residual for the filtered candidate

This result bounds the **second parameter derivative** of the retained
natural equations evaluated with the distinct numerical second field
`Z_epsilon`. It is not a second-derivative error bound relative to the exact
target, a relative-accuracy result, a finite-amplitude nonlinear residual,
or evidence for a cosmological origin mechanism.

All statements retain the actual background, actual first field, complete
quadratic bulk sources, affine shell law, and their previously authenticated
conditional premises. In particular, the actual first harmonic identity
`C1=0` holds on a neighbourhood. Merely selecting independent values from
an interval jet box would not establish that identity.

## Equations, parameter convention, and signs

The natural equations considered here are the trace-reversed Einstein–scalar
and scalar equations

```
E_ab = Ric_ab - (d phi)_a (d phi)_b - (2/3) V(phi) g_ab,
S = Box_g phi - Vprime(phi).
```

`E` is not the Einstein tensor `G_ab` itself. Expansions use
`g=g0+lambda q+(lambda^2/2) Z`, with the analogous scalar expansion. Thus
`E2[Z]=DE[Z]+QE` and `S2[Z]=DS[Z]+Qphi` are second parameter derivatives;
the coefficients of `lambda^2` are half these quantities.

Let the complete unit-principal scaled defect be

```
R[Z] = Z_xx-Z_rr-La Z-F.
```

The registered metric forcing is `-2 ell^2 rho^2 QEred`, whereas the scalar
forcing is `+ell^2 rho^2 Qphi`. The checker authenticates the primitive
homogeneous assembly and verifies these signs/scales as exact Laurent
polynomial identities in all 20 metric and three scalar output amplitudes.
The primitive scalar row already retains the natural `-DC[Z].dphi0`
coupling. It must not be removed or supplied a second time. Consequently

```
E2[Z] = R_metric/(2 ell^2 rho^2)
        + nabla_(a DC[Z]_b) + (QE-QEred),
S2[Z] = -R_scalar/(ell^2 rho^2).
```

Parentheses include one half. The actual source identity

```
QE-QEred = nabla_(a QC_b)
```

is verified by exact subtraction of every collected polynomial in all 20
metric output amplitudes, including both TT slots. This gives a route to a
direct finite natural residual bound **without** a new adapter for `QC_x`.
`QC` here is the actual full-metric-lowered gauge-covector Hessian. The
identity is its retained actual harmonic specialization; it is not the
arbitrary off-shell formula, which additionally contains the `2 P C1` term.

## Explicit charge for smoothing the candidate forcing

The actual target still uses `F`; only the auxiliary candidate uses its
one-sided causal average `F_epsilon`. With `W1q=W1_quadratic`,

```
R[Z] = R[W1q] - La(Z-W1q) + F_epsilon-F.
```

The signs above are for the defect defined as principal minus actual RHS.
The prior published forcing-defect cap may use its opposite sign; the norm
cap is unaffected. The exact new majorant is

```
R_H <= old_matched_defect_H
       + (2/3) k_actual X(Z-W1q) + 2 F_H.
```

The `(2/3) k_actual` factor is exactly the actual `La:X->H` majorant
`nA+nB+(3/4)nD`. The value factor `3/4` uses the inherited zero wedge,
which is preserved by past-time averaging. The checker imports neither a
zero-source target nor a filtered target equation. The source discrepancy
is charged in full. The new defect cap is approximately `3.574e22`.

## Coordinate-to-physical derivative bound

The authenticated full23 harmonic formula reads

```
P = ell rho DC[Z] = Rx Z_x + Rr Z_r + Zmat Z,
Zmat = alpha Za + eta Zeta + zeta Zzeta,
alpha=ell rho a, eta=ell tanh(ell x), zeta=ell sech(ell x).
```

The input `H` norm is the probability-normalized S3 L2 positive metric/scalar
bundle norm. The output covector `B` norm uses no extra factor two on mixed
covector components. The checker authenticates the preceding exact formula
report and independently recomputes

```
Rx Rx* + Rr Rr* = diag(7/4,7/4,7/4,7/4,7/4,7/4,1,1).
```

It uses the rational upper bound `15/8 > sqrt(7/2)` when both input jets are
bounded separately. Weighted Frobenius square bounds on `Za,Zeta,Zzeta`
are respectively `245/2,131/2,68`; the inclusive operator caps `12,9,9`
therefore suffice.

The radial background identities give

```
rho_r = -alpha rho,
alpha_r = ell^2 rho^2 b,
b = rho^-2-a^2+p^2/3.
```

The caps use authenticated actual `rho,a,p` intervals and the triangle
bound `|alpha_r| <= ell^2[1+rho_upper^2(a_abs^2+p_abs^2/3)]`.
Also `alpha_x=eta_r=zeta_r=0`, `|eta_x|<=ell^2` and
`|zeta_x|<=(3/2)ell^3`. Thus the checker bounds `Zmat`, `Zmat_x`,
and `Zmat_r` from these actual-background premises.

For candidate caps `V0,Vx,Vr,Vxx,Vxr,Vrr`, define

```
P0 = (15/8)max(Vx,Vr) + z0 V0,
Px = (15/8)max(Vxx,Vxr) + z0 Vx + zx V0,
Pr = (15/8)max(Vxr,Vrr) + z0 Vr + zr V0.
```

These bound `P,P_x,P_r` in the covector positive L2 norm. Since the physical
frame has `e_y=-(ell rho)^-1 partial_r` and
`e_0=(ell rho)^-1 partial_x`, differentiating the **covector components**
`DC=P/(ell rho)` gives

```
e_y DC = -(P_r+alpha P)/(ell^2 rho^2),
e_0 DC = P_x/(ell^2 rho^2).
```

The plus `alpha P` in the radial numerator must not be dropped. The frame
covariant derivative further includes the connection: the time-direction
connection has operator norm at most `|a|`; in the three sphere directions,

```
sum_a |Omega_a C|^2
  = (a^2+h^2)|C_sphere|^2 + 3(a C_y-h C_0)^2
  <= 3(a^2+h^2)|C|^2.
```

Indeed the difference of the right and left expressions equals
`3(h C_y+a C_0)^2+2(a^2+h^2)|C_sphere|^2`, so this is an exact
positive-square inequality, including the Lorentzian time connection signs.
The intrinsic S3 covector space has rough eigenvalues at most 24: scalar
normal/time amplitudes have degrees 0,2,4 and the two spatial gradient
amplitudes have eigenvalues 6,22. Hence
`||D_kappa C||_L2 <= sqrt(24)||C||_L2`.

Combining derivative-direction bounds by the triangle inequality gives

```
||nabla DC[Z]||_L2 <=
 [Px+Pr+{2 alpha_abs+5 ell+(7/4)(alpha_abs+eta_abs)}P0]
 /(ell^2 rho_lower^2).
```

The constants 5 and 7/4 are rigorous rational upper bounds for `sqrt(24)`
and `sqrt(3)`, respectively, with `sech<=1`. The norm of the symmetrized
gradient is no larger in this positive tensor norm.

## Numerical consequence and its limits

Using the separately authenticated candidate report with `epsilon=1/1024`,
the exact rational calculation gives readable strict upper bounds

```
||E2[Z_epsilon]||_positive,physical-frame,L2(S3) < 2.918e43,
||S2[Z_epsilon]||_physical-frame,L2(S3)          < 2.118e30.
```

These are continuous upper bounds on the declared retained domain. They are
**large absolute bounds**. No comparison scale or relative smallness is
established, and they do not show the candidate satisfies the equations
accurately. The Einstein bound is dominated by the conservative second-jet
bound in `nabla DC`; the scalar bound follows directly from the reduced
unit defect, because the scalar assembly is already natural. These are
positive comparison norms in a declared background frame, not
Lorentz-invariant scalar observables or physical forces.

The full actual `F_x` bound and C2 closeness to the unfiltered exact target
remain open. The filter also changes the candidate's interior initial
acceleration forcing; the previous vanishing initial defect cannot be
inherited. None of these estimates proves exact harmonic constraint
preservation, an exact nonlinear family, a finite-amplitude remainder,
a cosmological prediction, external originality, or a Big Bang origin.
