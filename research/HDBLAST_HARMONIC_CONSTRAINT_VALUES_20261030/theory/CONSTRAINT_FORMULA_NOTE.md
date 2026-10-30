# Linear harmonic constraint, frame normalization and exact norm identities

This note concerns the first-order operator `DC` in the second metric field.
It does not identify harmonic-gauge constraint values with the Hamiltonian or
momentum constraints, nor with a full natural Einstein residual. The same
actual first field and quadratic source must be retained when comparing two
second fields.

## Conventions and the full covector

Use the retained orthonormal parallel background frame

```
e_y=partial_y,
e_0=rho^(-1) partial_t,
e_a=(rho cosh(t))^(-1) bar_e_a,
signature=(+,-,+,+,+).
```

Here `bar_e_a` is a unit-S3 frame, `a=rho_y/rho`, `c=cosh(t)` and
`h=tanh(t)/rho`. The connection on fibre indices agrees with
`projection/inputs/linear_operator.py:omega`: along `e_0`, the normal/time
mixing coefficients are both `a`; along `e_a`,

```
nabla_a e_y=a e_a,       nabla_a e_0=h e_a,
nabla_a e_b= -a delta_ab e_y+h delta_ab e_0
              + intrinsic-sphere connection.
```

The relative harmonic covector linearization is

```
DC[k]_i=nabla^j k_ij-(1/2) nabla_i tr_g k.
```

The complete selected metric has scalar coefficients `K,L,N,J` in
`k_yy,k_y0,k_00,k_ab` at degrees `l=0,2,4`; scalar-derived mixed coefficients
`v,u` in `k_ya,k_0a` at degrees `l=2,4`; trace-free Hessian coefficients `H`
at degrees `l=2,4`; and two transverse trace-free tensor coefficients `T+/-`.
Write `lambda=l(l+2)` and

```
k_ab=J_l Y_l delta_ab + H_l S_ab[Y_l] + transverse-TT terms,
S_ab[Y_l]=D_a D_b Y_l+(lambda/3) gamma_ab Y_l.
```

The physical trace is `K-N+3J`. Unit-S3 identities are

```
Delta Y_l=-lambda Y_l,
Div(dY_l)=-lambda Y_l,
Div S[Y_l]=(2-2 lambda/3) dY_l.
```

For the last identity, commute the gradient with the scalar Laplacian using
`Ric(S3)=2 gamma`, then add `(lambda/3)dY`. Put
`beta=2-2 lambda/3`. Direct contraction of the connection gives

```
DC_y,l = (K_y+N_y-3J_y)/2 - L_t/rho
         +a(4K+N-3J)-3hL-lambda v/(rho c),
DC_0,l = L_y-(K_t+N_t+3J_t)/(2rho)
         +5aL-3h(N+J)-lambda u/(rho c),
DC_a   = sum_l [v_y-u_t/rho+5av-4hu
                +(-K/2+N/2-J/2+beta H)/(rho c)] D_a Y_l.
```

At degree zero, `v,u,H` and the angular output are absent. The scalar field
coefficient `w` does not enter a metric harmonic constraint. Both TT channels
are annihilated because their trace and divergence vanish. Thus the five
physical covector components have eight registered amplitudes: three normal
scalars, three time scalars, and two scalar-derived spatial covectors.

As a convention crosscheck, the longitudinal first metric
`K=-4u, N=-2u, J=2u, L=v=u_mixed=H=0` gives
`DC_y=-6(u_y+4au)Y`, matching the retained first-harmonic construction.

## Scaled coordinates and positive comparison norms

The retained scaled coordinates satisfy

```
x=t/ell,       r=z/ell,       partial_y=-(ell rho)^(-1) partial_r.
```

Multiplication by `ell rho` is pointwise multiplication of the covector;
it is not differentiation of a rescaled unknown. With

```
alpha=ell rho a,      eta=ell tanh(t),      zeta=ell/c,
P[k]=ell rho DC[k]=Rx k_x+Rr k_r+Z(alpha,eta,zeta)k,
```

the normal/time/gradient output rows are

```
P_R,l=-(K_r+N_r-3J_r)/2-L_x
       +alpha(4K+N-3J)-3eta L-lambda zeta v,
P_T,l=-L_r-(K_x+N_x+3J_x)/2
       +5alpha L-3eta(N+J)-lambda zeta u,
P_A,l=-v_r-u_x+5alpha v-4eta u
       +zeta(-K/2+N/2-J/2+beta H).
```

The output comparison norm is the normalized S3 average of
`C_y^2+C_0^2+sum_a C_a^2`. This is a positive norm in the declared background
frame, not a Lorentz-invariant scalar observable. If `S_l=average(Y_l^2)`,
its eight diagonal weights are `S_l` for each normal/time scalar and
`lambda S_l` for each gradient. The input H weights retain the symmetric
metric multiplicities: `K,N:S_l; L:2S_l; J:3S_l;
v,u:2lambda S_l; H:(2/3)lambda(lambda-3)S_l`, plus the registered scalar/TT
weights. In particular a mixed metric component's factor two is not present
in the corresponding output covector norm.

The code verifies the exact weighted matrix identity

```
Rx Rx* + Rr Rr* = diag(7/4 on six normal/time outputs,
                       1 on two gradient outputs),
```

where adjoints use the input H and output B Gram matrices. Thus

```
||Rx v+Rr w||_B^2 <= (7/4)(||v||_H^2+||w||_H^2).
```

When each first jet has H norm at most E, this yields the sharper derivative
bound `sqrt(7/2) E`, rather than the separate-operator triangle bound
`sqrt(5) E`. An exact degree-zero balanced normal/time adjoint witness
attains the squared factor `7/2`, so this derivative-only bound is sharp for
two independently bounded first jets. With the inherited zero-wedge estimate
`||e||_H<=(3/4)E`, any verified `z0>=||Z||` gives

```
||ell rho DC[e]||_B <= [sqrt(7/2)+(3/4)z0] E.
```

For a physical-frame covector bound, divide by `ell rho` using a certified
positive lower bound on rho. Bounds on scaled amplitudes alone must not be
reported as bounds on the unscaled covector.

## Exact normal constraint at the shell

Homogeneous scalar-degree shell rows in the physical y orientation are

```
By K=-4aK-2aN+6aJ+(4sigma1/3)w,
By N=-aK-(sigma1/3)w,
By J= aK+(sigma1/3)w,
M=-ell rho By.
```

Substitution into the normal DC row cancels the homogeneous metric and
scalar terms identically. Dirichlet `L=v=0` also gives `L_x=0` on the shell.
No assertion that `L_r` or `v_r` vanishes follows; those normal derivatives
remain in the time and spatial constraint rows.

For actual one-copy shell Hessians write `tr QR=-QR00+3QRJ`. Their source
packing is

```
GN=-2 tr QR-2 QC_y,
G_K=-ell rho GN,    G_N=2ell rho QR00,    G_J=2ell rho QRJ.
```

Consequently `G_K+G_N-3G_J=2ell rho QC_y`, and any integral-defined candidate
with these exact actual shell rows satisfies

```
DC[k]_y+QC_y=0  on the shell.
```

This direct source-packing identity does not require knowledge of the bulk
waveform or an assumption that a separate formal target has zero constraint.
It remains conditional on the inherited actual source and boundary identities.
It proves the normal harmonic constraint component, not the other four
components, physical Hamiltonian/momentum constraints or an Einstein residual.

## From sphere-averaged to all-angular-point bounds

The selected scalar harmonics lie in the entire degree `0,2,4` irreducible
hull on S3, of dimensions `1,9,25`. For normalized sphere measure the scalar
addition kernel is constant with value `1+9+25=35`. The exact and coexact
sectors are not being mixed here: the spatial constraint lies in the gradient
hull of degrees `2,4`. Dividing an L2-normalized scalar gradient by
`sqrt(lambda)` gives its normalized covector basis. Its evaluation kernel is
isotropic under the point stabilizer SO(3), has constant trace `9+25=34`, and
therefore is `(34/3)I_3`.

The full five-covector evaluation kernel is block diagonal

```
diag(35,35,(34/3)I_3).
```

Its largest eigenvalue is 35. The selected modal subspace is an orthogonal
subspace of this hull, so restricting it cannot increase the evaluation
operator norm. Uniformly at every angular point,

```
|DC[k](n)|_positive_background_frame <= sqrt(35) ||DC[k]||_B.
```

The more conservative trace-of-kernel factor `sqrt(104)` is unnecessary.
This finite-harmonic norm conversion is classical addition-theorem algebra;
it is not a new physical observable or a general infinite-mode estimate.

## Conditional constraint-value transfer and verification

For the same fixed quadratic source `QC`, if the retained formal completion
has `DC[U]+QC=0`, then

```
DC[W1]+QC=DC[W1-U].
```

The available C1 error can therefore bound the harmonic constraint's value
after inserting explicit coefficient norms and frame factors. This identity
does not certify derivatives of the constraint or a quantitative C2 error.
Those still require the differentiated source, higher initial-data errors and
the covariant conversion ledger. Smoothed nominal initial data are not assumed
to satisfy the physical initial constraints exactly.

`check_constraint_formulas.py` imports only the Python standard library.
It contracts manufactured normal-coordinate metric jets using the full
five-frame connection and Lorentz signs, separately from the modal formulas;
checks weighted adjoint identities and their sharp witness; records the shell
source cancellation and the angular hull dimensions; and rejects explicit
wrong-formula controls. The manufactured jets validate algebra and conventions,
not the physical premises or actual second field. No project scientific input
is mutated or imported by this checker.
