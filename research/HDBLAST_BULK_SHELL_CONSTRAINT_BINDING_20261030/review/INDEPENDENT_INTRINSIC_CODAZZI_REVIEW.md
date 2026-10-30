# Independent intrinsic-coordinate review of shell Codazzi source binding

30 October 2026. Internal mathematical and implementation review. This is an
exact source-to-proof binding for the unchanged registered model. It is not
external specialist assessment or a new Codazzi or nonlinear existence theorem.

## Result

The independently constructed Lorentzian shell divergence agrees with all five
normal-tangential **natural** Einstein-source coefficients of the authenticated
exported geometric map. The same five rows also obey the off-shell scalar law
before the first scalar physical residual is set to zero. No mismatch was found.

The independent checker imports neither the producer, its arithmetic helper nor
any project scientific module. It contains a separately written exact rational
Laurent ring and a direct decoder for the exported JSON monomials. It
authenticates the original map and analytic/convention sources before use and
reads them back afterward. All collected residual coefficients vanish exactly;
no point values, midpoint substitutions or zero-containing intervals are used.

## Direct intrinsic coordinate calculation

On the fixed shell the background radius `rho` is constant in tangential time
and angular coordinates. Write

```
gamma = rho^2[-dt^2 + c(t)^2 kappa],
c_t = c theta,       Hubble = theta/rho,       j = 1/(rho c).
```

At an arbitrary sphere normal-coordinate point the program builds the metric,
its first derivatives and every Christoffel symbol by the coordinate formula.
In particular `Gamma^t_ij=c^2 theta delta_ij` and
`Gamma^i_tj=theta delta^i_j`. The temporal inverse metric has the negative sign.
The scalar angular jet `Y`, all three `Y_i`, and the Hessian entries remain
formal. The only angular trace/divergence restrictions used are

```
tr_kappa(Hess Y) = -8Y,
div_kappa(Hess Y) = -6 dY.
```

The second identity follows from `Delta Y=-8Y` and the unit three-sphere Ricci
tensor `Ric=2kappa`. The third angular jets are kept free subject to these
contractions; this is an all-point tensor identity, not evaluation at one
selected angular direction.

For primitive physical-frame shell-source profiles `q00,q0,qa,qb`, the
coordinate covariant tensor is assembled independently as

```
QR_tt = rho^2 q00 Y^2,
QR_ti = rho^2 c q0 Y Y_i,
QR_ij = rho^2 c^2[qa Y^2 kappa_ij + qb Y Hess_ij(Y)].
```

The program forms every coordinate first derivative of this tensor, explicitly
computes `D^b QR_ba` using the Christoffel symbols, and subtracts the coordinate
derivative of the **Lorentzian** trace. Differentiation of the inverse metric
in the time derivative of the trace is retained. Division by `rho` for the
time covector and by `rho c` for angular covectors converts the results to the
background orthonormal frame.

With `T=3qa-8qb`, the resulting temporal primitive coefficients in
`(Y^2,|DY|^2)` are exactly

```
(-T_t/rho - (theta/rho)(3q00+T) - 8q0/(rho c),
 q0/(rho c)).
```

The angular primitive coefficients in `(Y DY, HessY.DY)` are exactly

```
((2q00-4qa+10qb)/(rho c) - q0_t/rho - 4theta q0/rho,
 qb/(rho c)).
```

These agree with the compact expressions proposed in the private planning
note. In the independent checker, however, the expected modal coefficients
are **extracted from the explicitly contracted coordinate tensor**, rather
than produced from those compact formulas. Scalar product projection uses
`Y^2=(1/6,1/2,1)` and `|DY|^2=(4/3,2,-4)`; gradient projection uses
`Y DY=(1/4,1/2)` and `HessY.DY=(1,-2)`.

## Physical shell and natural geometric map

The expected physical affine Hessian is independently reconstructed from the
existing outward one-copy shell formula:

```
beta = a h + sigma1 v/3,       sigma1=-2p,
A0 = a h^2/4 - sigma1 h v/6 - sigma2 v^2/6,
q00=-A0+2beta X, q0=-beta V, qa=A0-2beta L, qb=-2beta d.
```

Its tangential derivative is calculated with `a,p,sigma2` constant on the
background shell. The natural exported `QE/y0` and `QE/ya` maps are decoded
independently on the first shell jet class

```
A=h, D=-2X, E=V, L_frame=2L, M_frame=2d,
B=C=0 and their tangential time derivatives vanish,
D_y=-beta, L_frame,y=beta, E_y=M_frame,y=0,
s=v, s_y=p h/2-sigma2 v/2,
```

with tangential derivatives of the four metric normal rows only. Every unused
normal jet, including `A_y,B_y,C_y`, remains free and cancels. In particular
the normal harmonic row for `A_y`, any first bulk PDE or reduced acceleration,
the background bulk ODE and `c^-2+theta^2=1` are not needed. No identity valid
only at the shell is differentiated in the normal direction.

All three scalar temporal and both angular-gradient natural coefficients
equal the corresponding intrinsic shell divergence. This is a finite
implementation join to an existing geometric identity; the geometric map's
all-point export theorem and the analytic physical-shell derivation are
retained inputs.

## Off-shell scalar cancellation and second variations

The exact physical identity is obtained before setting the scalar shell
residual to zero. In four shell dimensions, put

```
K=R+sigma gamma/6,       n phi=-sigma'/2+P.
```

The geometric Codazzi contraction contributes
`div R-d tr R-(sigma'/2)d phi`; subtracting the Einstein-scalar mixed matter
term `(n phi)d phi` gives

```
E_na = div R-d tr R-P d phi.
```

Deleting that matter cancellation is a genuine rejected control. At second
parameter order, `P*dphi` contributes `P2*dphi0+2P1*dphi1+P0*dphi2`.
With the actual background `P0=0,dphi0=0` and first metric shell residual zero,
the five exported natural-map coefficients therefore satisfy

```
QE_y0 = div QR_0 - 2P1 v_t Y^2/rho,
QE_ya = div QR_a - 2P1 v Y DY/(rho c),
P1=s_y-p h/2+sigma2 v/2.
```

The checker verifies these five **actual map** equalities with `P1` free, and
only then separately verifies the physical `P1=0` specialization. It does not
pretend that altering `QP` after the surviving scalar term has vanished is an
effective corruption test.

The parameter convention is `g=g0+epsilon q+epsilon^2 k/2`. For the affine
boundary metric with zero mixed normal trace, the raised unit normal has
coefficient series `1-epsilon h/2+3epsilon^2 h^2/8`. The checker verifies its
normalization and multiplies it by the ordinary natural residual series. The
second pullback component is

```
QE_na - h E1_na + (3/4)h^2 E0_na.
```

The extra terms vanish with the background and first physical mixed residuals;
the latter follow from the same first Codazzi identity. A perturbed tangential
metric, raising or shell connection similarly multiplies `R0` or `R1`, which
are zero. These are explicit lower-residual premises, not omitted terms for
arbitrary nonzero first residuals. No new independent normal derivative of a
shell identity is introduced.

## Review of the main producer

The main checker `check_natural_shell_codazzi.py` and its arithmetic helper
were read independently. Their natural-versus-reduced source distinction,
ordinary time/frame factors, scalar off-shell `-2P1 D eta` term, Lorentzian QR
trace, minimal coherent normal-jet restrictions and source authentication are
consistent with the direct coordinate calculation above. The registered
`source_coefficients.evaluate_sources` call is a third implementation join;
it is not used as the independent expected divergence. The two imported
scientific modules have only standard-library dependencies on the exercised
paths and are authenticated before import.

The main producer's 35 named rows include scalar/tensor QR adapter joins,
five off-shell and five physical Codazzi modal rows, and independent-decoder
versus registered-evaluator comparisons. These are verification bookkeeping,
not 35 new theorems. No mathematical mismatch was found in the reviewed code.

The final main producer reviewed here has checker SHA256
`c3d260c951d54a3e47bf2b6a62d096d49e5dff2cb8811c0fed04d69546754b31`,
helper SHA256
`178b0a4e58db63eae6c24112ed9d91bf430f70dca9f6f8126c5e352e429f8fc2`,
and complete normal/optimized report SHA256
`e1a44b5c3125b9f15f7dc708a97ac0940189c3be393af25c8e1f7a9f9e82514b`.
Its final positive-time-trace corruption genuinely replaces
`tr QR=-q00+T` by `+q00+T`; the resulting temporal expression changes by
`-2q00_t/rho` and is correctly rejected. Its report explicitly retains
`R1=0` in the scalar-off-shell specialization; leaving only `P1` free must
not be described as leaving every first shell residual free.

## Verification and limits

The independent checker has 23 named exact identities and rejects eight
corrupted fresh subprocesses per execution mode. The normal and optimized
complete reports are byte-identical. Controls cover a positive temporal
divergence sign, omitted shell connection terms, a wrong background scalar
tension sign, omitted anisotropic source, an extra warp term in the frame
normal derivative, reduced instead of natural source, halved second variation
and omitted off-shell scalar cancellation. Every corrupted run must fail at
its intended nonzero identity and create no successful report.

The audit binds sources on the coherent first physical shell-jet class. It
does not independently prove that the actual evolved first wave belongs to
that class, the global Bianchi identity, smooth boundary regularity or
subsidiary uniqueness. Existing proofs of those statements remain analytic
premises. The induced subsidiary Robin conversion is a separate convention
check; a new incoming primal constraint datum is not added. No finite-amplitude
nonlinear spacetime, small physical error, novelty priority or higher-dimensional
Big Bang origin follows from these identities.
