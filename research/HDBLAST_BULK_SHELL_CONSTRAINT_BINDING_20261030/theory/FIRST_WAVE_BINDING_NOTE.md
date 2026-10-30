# Exact binding of the registered first gauge wave to the harmonic operator

30 October 2026. This is a finite implementation-binding calculation for the
unchanged actual first field. It does not solve a new PDE, alter an actual
source or numerical candidate, or establish a nonlinear cosmological model.

The first-wave existence/regularity argument, the coordinate Lie derivative
and the raw gauge-vector operator were already derived in the retained
October 16 capsules. The new calculation connects the current
`numeric/actual_first_jets.py::frame_amplitudes` implementation to those
formulas by a separate coordinate-tensor and frame-connection construction.
It also verifies the same registered operator's initial specialization.

## Exact input and independence boundary

The checker uses standard-library sparse Laurent polynomials over the rational
numbers. Formal scaled derivative jets of the raw functions `f,b,d,W` commute;
there is an explicit guard against differentiating past their declared order.
The background identities used in the calculation are

```text
rho_r = -ell rho^2 a,        a_r = -ell rho a_y,
c_x = ell H c,              H_x = ell(1-H^2),
theta_x = -(3/2)ell H theta, R_r = Rr,
theta_r=R_x=rho_x=0,        c^(-2)=1-H^2.
```

Here `c=cosh(t)`, `H=tanh(t)`, `theta=c^(-3/2)` and `R=alpha/rho^2`.
The identity remains algebraic in the static derivative `a_y`; no actual
background field equation is substituted to manufacture an operator
cancellation. Nonzero `rho`, `ell` and `c` justify the Laurent divisions.

Before source AST execution, the full 66-member first-wave capsule is
authenticated using its external manifest SHA256 and every member's exact
size and SHA256, including closed-tree and symlink checks. Additional pins
authenticate the current actual reconstruction and the retained pullback
proof note. The same bytes are checked after use. No project module is
imported: only authenticated function/assignment ASTs are compiled in the
formal arithmetic namespace.

The implementation side supplies:

- the entire `frame_amplitudes` function;
- the entire `physical_derivatives` function;
- the registered `A,B,D,Ps,P0` matrix/source assignments in
  `check_global_gauge_coefficients.py`;
- the original `A,B,D,force` assignments in `check_gauge_operator.py`.

The independent side supplies coordinate metric/Lie tensor construction and
the complete five-dimensional parallel-frame connection contraction. It
does not use a retained printed `DC` row or copied matrix formula to produce
the expected harmonic covector.

## Coordinate pullback and the actual frame sectors

For raw contravariant amplitudes

```text
xi^y=fY, xi^t=bY, xi^i=d D^iY,
u=R theta W, v=B theta(ell^(-1)W_r+kW),
h'=h_long+Lie_xi(g0), scalar'=vY+p fY,
```

the checker independently forms all 15 symmetric coordinate metric
components from

```text
(Lie_xi g0)_ij = xi^k partial_k g0_ij
                +g0_kj partial_i xi^k+g0_ik partial_j xi^k.
```

It converts them to the background frame
`e_y=partial_y`, `e_0=rho^(-1)partial_t`,
`e_a=(rho c)^(-1)bar_e_a` and compares them to the compiled implementation.
This recovers every sector

```text
A=-4u-2f_r/(ell rho),       Bmix=f_x/(ell rho)+b_r/ell,
Cmix=f/(rho c)-c d_r/ell,   D=-2(u+af+b_x/ell),
E=c d_x/ell-b/c,            L=2(u+af+Hb), M=2d,
scalar=v+pf.
```

The radial gauge amplitude has not been rescaled by `rho`. The `b_x/ell`
clock derivative in `D` and the `H b` clock value in `L` are distinct.
The source reconstruction coefficient named `B` in the scalar formula is
not the mixed metric sector `Bmix`.

The independently differentiated ordinary radial conversion also agrees
with the compiled function in all six returned slots. In particular

```text
partial_y F = -(ell rho)^(-1)F_r,
partial_yy F = (ell rho)^(-2)F_rr + a(ell rho)^(-1)F_r.
```

Only exact formal smooth functions are differentiated. No Taylor remainder
or interval value is treated as a differentiable function.

## Independent contraction of the first harmonic covector

The checker constructs `DC_i=nabla^j h'_ij-(1/2)nabla_i tr(h')` using the
timelike inverse sign `(+,-,+,+,+)` and the complete parallel-frame
connection. At an arbitrary sphere-normal frame point it retains free
`Y,Y_a,Hess_ab(Y)` generators with their exact trace and third derivative:

```text
sum_a Hess_aa(Y)=-8Y,
D_c Hess_ab(Y)=-delta_ca Y_b-delta_cb Y_a-2delta_ab Y_c.
```

These are identities for the actual quadrupole on the unit three-sphere,
not a finite angular point sample. The calculation retains the angular
Hessian as an untraced Hessian, matching the seven-sector first-field
adapter, and includes the time/angular FRW and radial warp connection terms.

Let the actual registered scaled matrices be `P,Q,Rmat` and define the
source-free residual

```text
Res = U_xx-U_rr-P U_x-Q U_r-Rmat U, U=(f,b,d),
s_f=-6ell^2 rho^2(u_y+4au), s_b=s_d=0.
```

The complete collected coefficients satisfy

```text
ell^2 rho^2 DC_y = (s_f-Res_f)Y,
ell^2 rho   DC_0 = Res_b Y,
ell^2 rho/c DC_a = -Res_d D_aY.
```

All 27 matrix entries are also compared between the authentic global
coefficient evaluator and the authentic original raw operator after the
required `ell,ell,ell^2` scalings. The actual `Ps W_r+P0 W` source matches
the same `s_f`, retaining its negative ordinary-y source sign and both
`ell^2 rho^2` factors.

Consequently the exact actual solution of the registered PDE `Res=s`
has `DC[h']=0` on its physical domain. The existence, smoothness and exact
PDE identity of that actual solution remain the retained analytic premises.
This finite calculation binds the implementation to that implication; it
does not infer a PDE identity from interval values containing zero.

## The same operator supplies the corrected initial one-jet

Zero gauge Cauchy functions mean `f=b=d=f_x=b_x=d_x=0` throughout the initial
radial interval. Their spatial derivatives and spatial derivatives of the
zero velocities also vanish. The same audited operator therefore fixes

```text
f_tt=-6rho^2(u_y+4au), b_tt=d_tt=0.
```

Specialization of the compiled reconstruction then gives

```text
A=-4u, D=-2u, L=2u, Bmix=Cmix=E=M=0,
A_t=-4u_t, D_t=-2u_t, L_t=2u_t,
(Bmix)_t=-6rho(u_y+4au), (Cmix)_t=E_t=M_t=0,
scalar=v, scalar_t=v_t.
```

Since the value equalities hold as functions of the initial radial
coordinate, their spatial derivative equalities follow under the same
retained smoothness. This is the corrected first one-jet used by the
initial physical/ADM audit; setting `f_tt=0` merely because `f=f_t=0`
would lose its mixed time metric coefficient. The independent companion
review additionally binds this table directly to the earlier initial
`QC` checker's source AST; that separate join is not falsely counted among
this producer's 77 rows.

The physical shell mixed rows are checked separately by substituting
`f=f_x=b_r=d_r=0`. Both `Bmix` and `Cmix` vanish and the scalar pullback is
unchanged. The tangential gauge values `b,d` remain free. No zero values
for their Neumann components, no reflected smooth background coefficients,
and no extra incoming constraint rows have been imposed.

## Verification and remaining scope

Both normal and optimized fresh processes pass 77 named coefficient/array
joins with identical complete report bytes. Each mode rejects 15 actual
corrupted subprocesses at its expected residual, without a successful
output report. Controls cover the scalar pullback sign, missing radial
frame factor, rescaling raw `f`, interchanging clock terms, the Lorentz
sign, the angular Hessian divergence, both raw radial drift factors,
the time `ell` factor, source sign and separate missing `ell^2`/`rho^2`
factors, the angular frame factor, the variable-rho second derivative,
and zeroing the initial acceleration. Counts are bookkeeping, not
independent theorems.

The source-role ledger relating general/background-lowered and harmonic
quadratic source formulas is not newly calculated here. The full-metric-
lowered `QC_actual` formula is distinct from background lowering; the
actual harmonic simplifications may now use this audited operator
implication together with the retained actual-solution premise. The shell
Codazzi and subsidiary uniqueness arguments are separate. No new PDE solve,
finite-amplitude nonlinear family, small physical error, cosmological
origin, novelty priority or external specialist approval is established.
