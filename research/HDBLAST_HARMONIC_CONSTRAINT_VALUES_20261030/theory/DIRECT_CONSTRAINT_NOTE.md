# Direct harmonic-constraint values at the selected shell event

The candidate is the same integral-defined `W1` as the preceding causal-bulk
checkpoint. At `(x,r)=(9/8,0)`, its complete 23-component first-jet intervals
are authenticated from that checkpoint's target-transfer report. This audit
assembles the second-parameter relative-harmonic constraint directly:

```
mathfrak_C^(2)[W1] = DC[W1] + Q_C[q,q].
```

The actual first field q and actual quadratic source are fixed. The source
is retained using its genuine uniform physical component intervals. It is
not assigned a point value, sign, or zero. The direct audit does not assume
that a separate formal target U has vanishing constraint. It remains
conditional on the inherited actual background, first-field source map,
finite angular basis, candidate definition and boundary-source consistency.

## Five-dimensional frame derivation

Use an orthonormal background frame with radial index y, timelike index 0,
and three sphere-tangent indices a. Its signature is `(+,−,+,+,+)`. The
registered metric amplitudes are

```
k_yy=K_l Y_l, k_y0=L_l Y_l, k_00=N_l Y_l,
k_ya=v_l D_a Y_l, k_0a=u_l D_a Y_l,
k_ab=J_l Y_l gamma_ab+H_l HessTF_ab(Y_l)+TT_ab.
```

The scalar field amplitudes w do not enter the metric-only gauge functional.
The two TT slots have zero trace and divergence and therefore contribute
zero to DC, without being removed from the authenticated candidate input.
Put `a=rho_y/rho`, `h=tanh(t)/rho`, `c=cosh(t)`, and
`lambda_l=l(l+2)`. The scalar harmonics have l=0,2,4; the gradient and HessTF
slots exist only at l=2,4. The identities on the unit three-sphere are

```
Delta Y_l=−lambda_l Y_l,
div HessTF(Y_l)=(2−2lambda_l/3) dY_l.
```

At a sphere-normal frame, contract the complete covariant first jet:

```
DC[k]_B = sum_A sigma_A (nabla_A k)_AB − (1/2) E_B tr_g k.
```

The nonzero extra warped-frame connections are
`nabla_0 E_y=a E_0`, `nabla_0 E_0=a E_y`,
`nabla_a E_y=a E_a`, `nabla_a E_0=h E_a`, and
`nabla_a E_b=delta_ab(−a E_y+h E_0)` at that event. This gives, degree by degree,

```
DC_y = (K_y+N_y−3J_y)/2 − L_t/rho − lambda v/(rho c)
       +4aK+aN−3aJ−3hL,
DC_0 = L_y − (K_t+N_t+3J_t)/(2rho) − lambda u/(rho c)
       +5aL−3h(N+J),
DC_a = [v_y−u_t/rho+5av−4hu
        +((-K+N−J)/2+(2−2lambda/3)H)/(rho c)] D_a Y.
```

The radial coordinate is inward conformal distance, so the authenticated
coordinate convention is `dy=−rho ell dr`, while `t=ell x`. Thus every
physical covector row equals `1/(ell rho)` times its scaled row. With
`alpha=ell rho a`, `eta=ell tanh(ell x)` and `zeta=ell sech(ell x)`, the
scaled rows are

```
R = −(K_r+N_r−3J_r)/2 − L_x
    +alpha(4K+N−3J)−3eta L−lambda zeta v,
T = −L_r−(K_x+N_x+3J_x)/2
    +5alpha L−3eta(N+J)−lambda zeta u,
A = −v_r−u_x+5alpha v−4eta u
    +zeta[(-K+N−J)/2+(2−2lambda/3)H].
```

The physical source Q_C is added after this conversion. It is not multiplied
by the wave-forcing normalization `−2 ell^2 rho^2`; that normalization belongs
to the reduced Einstein evolution forcing, not the gauge covector.

## Full lowering and source cancellation

The projected source map states explicitly that Q_C is the second variation
of the **full metric-lowered** gauge covector. For a path
`g(eps)=g0+eps q+eps^2 k/2`, its second derivative is `DC[k]+Q_C[q,q]`.
Along an affine first-pair path the lowered functional can be written
`C_B(g)=g^(CD) S_BCD(g)`, where S is linear in the background-covariant metric
first jet. Its affine quadratic term is `−2 q^(CD) S_BCD(q)`.

If a quadratic contravariant source is lowered only with g0, the missing
second-variation lowering term is `2 q_BA DC[q]^A`. The actual first pair is
relative harmonic, so that term vanishes. The registered projection instead
already supplies the full lowered source, making its convention explicit.
For any two second fields with the same actual first field, the same Q_C
cancels in their constraint difference. That comparison identity is useful,
but is not required to form the direct candidate intervals in this audit.

## Shell normal component and Dirichlet derivatives

The inherited exact shell-source identity is
`GN=−2 tr(QR)−2 Q_C,y`. Together with the same actual mixed shell law used to
define W1, it gives `mathfrak_C^(2)[W1]_y=0` at the shell, degree by degree.
The homogeneous normal identity also follows directly from the three
metric Robin rows: the combination `−(K_r+N_r−3J_r)/2` cancels
`alpha(4K+N−3J)`, including the scalar-tension coupling. This is a value of
the normal covector component; it is not a normal derivative of the complete
constraint.

The five Dirichlet field traces L and v and their tangential derivatives
vanish. Their **normal derivatives do not vanish**. The time and angular
constraint rows retain `−L_r` and `−v_r`. The independent exact control
explicitly rejects formulas that omit those normal terms.

## Numerical enclosure and norm

All eight direct modal intervals contain zero. Consequently this calculation
proves no nonzero harmonic-defect sign. The raw normal intervals are also
retained in the report, then intersected with the independently exact normal
boundary identity. The physical covector norm is the positive frame norm,
integrated against the unit-S3 probability measure. Its weights are

```
R_l,T_l: g_l,     A_l: lambda_l g_l,
g_0=1, g_2=1/12, g_4=19/720.
```

The gradient output is a covector rather than a symmetric mixed metric slot;
there is no extra factor two in its Gram weight. Direct exact endpoint sums
give the integer upper bounds

```
raw norm       <= 1318623212395798958550689,
normal-refined <= 1318560930352165554331780,
both           < 1.319e24.
```

These bound a first-order gauge-constraint **value** at one event. They do
not bound its derivatives, all Hamiltonian or momentum constraints, the full
natural Einstein equations, or a C2 field error. Zero inclusion is not proof
that the constraint vanishes in the remaining components. The bound is in
registered physical-frame units and is not a physical-force measurement,
nonlinear amplitude limit, or observational result.

## Independent exact tests

`check_direct_formula_audit.py` imports no project module. It compares the
modal formulas with a separate five-index connection contraction for 54
individual field-jet inputs. The sphere jets are those of explicit degree-two
and degree-four harmonic polynomials restricted from R4 at the north pole;
their Ricci commutation yields the HessTF divergence coefficient independently.
It also checks the exact homogeneous normal boundary cancellation and the
full lowered second-variation source convention. Three wrong-formula controls
reject dropped Dirichlet normal terms, reversed radial orientation, and an
`eps^2` coefficient substituted for a second parameter derivative.

The direct numerical checker authenticates the upstream driver and source
reports before imports and reauthenticates them afterward. Its normal and
optimized reports are byte identical. Independent formula review reproduced
all eight physical interval endpoint sums and both weighted norm squares.
