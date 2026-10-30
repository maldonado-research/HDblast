# A differentiated error gate for the unchanged original-source candidate

This calculation concerns `e=U-W1_quadratic`, where `U` is the retained
second-coefficient target with its **original full actual source F**.
`W1_quadratic` is the PR100 finite candidate, also called `W1_dagger`.
It does not use the causally averaged candidate `Z_epsilon` introduced in
PR101. The latter report is used here only for authenticated initial
absolute data and the frozen-stage forcing ledger, which predate its
source averaging and are identical for both candidates.

Without a complete actual `F_x` report, the checker produces an
affine input gate rather than a numerical C2 error certificate. A supplied
report must explicitly refer to original full23 `F_x` in the scaled
coordinate `x=t/ell`; a derivative bound for averaged forcing, a response
norm, or a single source summand cannot fill this role.

## Domain, norms and fixed operators

Use the retained physical half-domain
`D+={0<=x<=3/2,max(0,x-5/4)<=r<=7/4-x}` and positive full23 norm H.
The norm is `X(V)=max(sup||V||/3,sup||V_x||,sup||V_r||)`.
Fields vanish on the inherited exact wedge `x+r<=1`; hence
`sup||V||<=3X(V)/4`. Let

```
La=A partial_x+B partial_r+D,
L0=A0 partial_x+B0 partial_r+D0,
P=(I-L M PN)^(-1),
b=||M||/4,  ||P||_X<=1/(1-b),
q=3409/8000000,
c=a+b_bulk+(3/4)d,
comm=a1+b1+(3/4)d1.
```

Here `a,b_bulk,d` denote actual bulk operator caps; the script uses separate
variable names for `b_bulk` and the boundary inverse constant. Derivative
caps `a1,b1,d1` are the complete actual coefficient time derivatives.
The actual `B_x=M_x=0`. No coefficient is commuted through parity or the
causal propagator. All physical source coupling occurs before output
reflection.

Write `Bs` for the same free wave and

```
f0=L0 Bs,
W0=P[Bs+Kf0+LG],
W1=P[Bs+K(La W0+F)+LG].
```

The targets, background, first field, all sources, initial data and actual
Robin matrix are unchanged. The candidate is a finite update, not an
asserted exact solution of the self-coupled equation.

## A cancellation-preserving bound for the differentiated correction

Set `E0field=W0-Bs` and `C=W1-W0`. The already verified quantities
`E0X=X(E0field)` and `C1=X(C)` are used as upper bounds. Both differences
have zero Cauchy displacement and velocity. The free wave has homogeneous
Dirichlet/Neumann reflection, so its PN normal derivative is zero on the
physical shell. Therefore

```
(E0field)_xx-(E0field)_rr=f0,
PN (E0field)_r=M PN E0field+M PN Bs+G,
PD E0field=0.
```

Let `U2,V,V1` bound `||u_s,rr||,||v_s||,||v_s,r||`, and let `B_s` bound
`X(Bs)`. The free-wave second-jet cap is `Qfree=U2+V1`.
Constant frozen coefficients give

```
f00=a0 V+b0 ||u_s,r||+d0 ||u_s||,
f0x=(a0+b0)Qfree+d0 B_s.
```

For `Y0=(E0field)_x`, differentiation in time gives zero initial
displacement but initial velocity `f0(0)`. Its shell source is
`M PN Bs_x+G_x`. Thus

```
X(Y0) <= [f00+(3/2)f0x+||M|| B_s+Gx]/(1-b).       (1)
```

Dropping `f0(0)` because the original difference has zero Cauchy data would
be incorrect. For example, a zero-data wave driven by constant forcing
has nonzero initial acceleration. Equation (1) uses the complete actual
`G_x` packing; it is not a boundary-error term in the final target-error
equation.

The bulk correction obeys `C=P K r0`, where

```
r0=La(W0-Bs)+(La-L0)Bs+F.
```

Let `delta_a,delta_b,delta_d` bound the actual-minus-frozen operator
differences, and let `Fx` denote the still-to-be-supplied genuine full
pointwise `sup||F_x||_H`. Differentiating this difference form yields

```
r0x <= (a+b_bulk)X(Y0)+(d+comm)E0X
       +(delta_a+delta_b)Qfree+(delta_d+comm)B_s+Fx.  (2)
```

The complete time derivative of `La-L0` is the actual time derivative of
`La`, because `L0` is constant. This form preserves the large free-wave
cancellation; summing unrelated absolute W0 and Bs second-jet bounds is
unnecessarily loose.

On the initial face `W0` and `Bs` have equal first jets, so

```
r00 <= delta_a V+delta_b ||u_s,r||+delta_d ||u_s||+FH.
```

The full actual source cap `FH=17387513156573` is retained here: `F(0,r)`
is generally nonzero inside the initial packet. For `Y1=C_x`, the initial
displacement is zero and initial velocity is `r0(0)`. Its boundary source
is zero, because W1 and W0 use the same actual affine G. Hence

```
X(Y1) <= [r00+(3/2)r0x]/(1-b),                    (3)
F1 >= sup||(r1)_x||
   <= (a+b_bulk)X(Y1)+(d+comm)C1,                 (4)
r1=La C.
```

There is no `G_x` term in (3) or in the final target error boundary law.
It already appears where needed in (1).

## Original-target differentiated error and initial acceleration

The exact error equation is

```
e_xx-e_rr=La e+r1,
PN e_r=M PN e,  PD e=0.
```

Although `C_xx(0)=r0(0)` need not vanish, C and all its initial first jets
do vanish. Since La is first order, this proves the exact identity
`r1(0)=0`. It is distinct from the filtered candidate's generally nonzero
initial error forcing; no filtering identity is imported here.

For initial error caps `Eu,Eur,Ev,Eurr,Evr`, the correct acceleration and
free-wave inputs for `v=e_x` are

```
W=Eurr+a Ev+b_bulk Eur+d Eu,
J=max((Ev+(3/2)W)/3,Evr+W).
```

The three lower initial error caps are exactly the individual caps
reconstructed for PR100's matched extension. The continuous PR100 higher
errors supply `Eurr<6.219e33` and `Evr<9.672e24`. Actual-minus-smoothed
errors are used, not nominal absolute second derivatives.

Let the retained complete C1 error upper bound be E. The differentiated
causal inverse gives

```
Z=[J+(3/2)(F1+comm E)]/(1-q),
sup||e_xx||,sup||e_xr|| <= Z,
sup||e_rr|| <= Z+c E+F0,
```

where `F0=c C1` is the already verified complete original reduced-defect
magnitude. The affine shell source cancels exactly, so boundary derivative
error is zero. The constant actual M introduces no `M_x` term.

The physical-side derivative construction retains the smooth actual
source/background, initial smoothing and flat-corner premises. A doubled
odd output source can jump at the shell, so no C2 extension across that
artificial interface is asserted. These are uniform one-sided
scaled-coordinate estimates on the original physical half-domain.

## Exact affine consequence before actual Fx is supplied

The exact rational output is of the form

```
sup||e_xx||,sup||e_xr|| <= Z_const+Z_Fx Fx,
sup||e_rr|| <= R_const+Z_Fx Fx,
Z_Fx=(3/2)^2(a+b_bulk)/[(1-b)(1-q)].
```

For display only, the verified constants are approximately

```
Z_const=6.221556148097404e33,
R_const=6.221556148140732e33,
Z_Fx=0.000482444763679284.
```

Thus an actual Fx below approximately `9.20e32` would retain a common
readable strict bound below `6.222e33`. This statement is a threshold,
not a claim that Fx has that bound. The exact threshold is recorded as a
rational in the report.

The baseline is dominated by the initial displacement second-derivative
error. These are conservative absolute error bounds, not lower bounds on
the true error or demonstrated relative physical accuracy. Conversion to
physical/covariant derivatives requires ell, rho, product-rule and frame
connection factors. A coordinate C2 error bound alone does not prove all
initial constraints, exact subsidiary constraint propagation, a compatible
nonlinear family, a finite-amplitude remainder or a Big Bang origin.

## Applying the complete original-source derivative bound

The new authenticated complete actual-source derivative report supplies

```
Fx = sup||F_x||_H <=905405545954284 <9.055e14.
```

It explicitly concerns the same original full23 source, not the averaged
candidate. Substitution in the exact affine formulas gives

```
sup||U_xx-(W1_quadratic)_xx||_H <6.222e33,
sup||U_xr-(W1_quadratic)_xr||_H <6.222e33,
sup||U_rr-(W1_quadratic)_rr||_H <6.222e33.
```

This closes a conditional, continuous, physical-side scaled-coordinate C2
**absolute error upper bound for the unchanged original-source candidate**.
It does not demonstrate that this is an acceptably small physical error.
The exact majorants are approximately `6.2215561481e33`; their affine Fx
contribution is only about `4.37e11`. In this upper-bound ledger the
initial free-wave term accounts for more than `0.999999999` of the numerator.
This fraction describes the proof majorant, not a decomposition of the
unknown true error. A zero denominator or comparison to a nonzero physical
scale is not supplied, so no relative accuracy claim follows.

The inner cutoff joins and coarse initial interval dependencies therefore
merit attention before more bulk iterations. The existing exact jump table
has zero second-derivative jump in all components at `65/64`; that join must
not be counted as a nonzero curvature jump merely because it is a branch
boundary.
Reducing a majorant alone is not proof of reducing the actual error; a
changed higher-order cutoff would define another initial candidate and
would need all corresponding estimates re-established.

## Reproduction

Run `check_original_target_c2.py` with a full checkout as `--repo-root`.
Override the prior capsule paths only for sparse local worktrees. Without
`--fx-report` it explicitly reports the open affine gate. With an actual
derivative report, provide its exact SHA-256 in `--fx-sha256`; the program
requires its full original-source scope and emits numerical error caps.
Normal and optimized Python must produce identical output. All arithmetic
and acceptance decisions use exact rationals and explicit checks.
