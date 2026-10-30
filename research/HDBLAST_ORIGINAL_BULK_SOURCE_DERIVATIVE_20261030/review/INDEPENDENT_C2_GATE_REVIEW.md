# Independent review of the original-target second-coordinate error gate

30 October 2026. This is an independent internal analytic and arithmetic
review. It does not review novelty, supply outside specialist approval, or
prove the retained background and first-field premises independently.

The candidate is the unchanged PR100 `W1_quadratic`, not PR101's causally
averaged auxiliary candidate. The actual target uses the original unsmoothed
complete `F` and the same actual shell source `G` and constant Robin matrix
`M`. The positive full23 norm and scaled coordinates are unchanged.

## First differentiated principal construction

Write `Bs` for the smoothed quadratic-continued free wave, `La` for the
actual homogeneous first-order operator, `L0` for its constant frozen
counterpart, and `P=(I-L M PN)^(-1)`. Set

```
W0=P[Bs+K L0 Bs+LG],
E0=W0-Bs,   f0=L0 Bs.
```

`Bs` satisfies the reflected homogeneous principal boundary conditions.
Consequently, on the physical shell,

```
PN E0_r=M PN E0+M PN Bs+G,   PD E0=0.
```

The Cauchy data of E0 vanish. Differentiating the principal equation in
time gives E0_x initial displacement zero and initial velocity `f0(0)`.
Its affine boundary source is `M PN Bs_x+G_x`; neither term vanishes merely
because E0's initial data vanish. If `B=X(Bs)`, `f00>=||f0(0)||`, and
`f0x>=sup||f0_x||`, the original principal-wave estimates imply

```
X(E0_x) <= [f00+(3/2)f0x+||M|| B+Gx]/(1-b).
```

The `f00` free-wave bound follows from zero displacement, velocity f0(0),
and T=3/2: its derivative branch is f00 and its value branch is f00/2.
The same quarter-length active shell history gives b=||M||/4. No bulk
feedback inverse is required at this stage, because f0 is prescribed.

## Second differentiated principal construction

Define the actual update correction C and source r0 by

```
C=W1_quadratic-W0=P K r0,
r0=La E0+(La-L0)Bs+F.
```

This form preserves the cancellation between the large absolute W0 and Bs
second jets. Let `a,b0,d` be actual A,B,D norm caps, `da,db,dd` their
discrepancy caps relative to L0, and `m1=a1+b1+3d1/4` the actual time
commutator. If Q bounds the free-wave second jets, then

```
||r0_x|| <= (a+b0)X(E0_x)+(d+m1)X(E0)
            +(da+db)Q+(dd+m1)B+||F_x||.
```

This is a product differentiation before bounding, with coefficients left
in their original order. In particular, there is no parity/matrix or
boundary-inverse/operator commutation. Since W0 and Bs have the same
initial first jets,

```
||r0(0)|| <= da||v_s||+db||u_s,r||+dd||u_s||+FH.
```

C has zero Cauchy data and homogeneous actual Robin boundary. Thus C_x
has zero displacement and initial velocity r0(0), with zero affine shell
source. Its first-jet bound is

```
X(C_x) <= [||r0(0)||+(3/2)sup||r0_x||]/(1-b).
```

The original reduced error forcing is `r1=La C`, so

```
||r1_x|| <= (a+b0)X(C_x)+(d+m1)X(C).
```

The initial value of r1 is exactly zero because C,C_x,C_r vanish initially.
C_xx need not vanish initially; its generally nonzero r0(0) contribution is
already present in the C_x estimate. The final cancellation r1(0)=0 does
not license dropping either earlier differentiated forcing acceleration.

## Final target-error gate

For `e=U_target-W1_quadratic`, the actual sources and boundary laws cancel
to give

```
e_xx-e_rr=La e+r1,
PN e_r=M PN e,   PD e=0.
```

Its time derivative obeys the same complete actual linear feedback with
source `r1_x+A_x e_x+B_x e_r+D_x e`. With the retained C1 error E and
full contraction `q=3409/8000000`,

```
w0=e0_rr+A0 v0+B0 e0_r+D0 e0,
J=max((||v0||+(3/2)||w0||)/3,||v0_r||+||w0||),
Z=[J+(3/2)(sup||r1_x||+m1 E)]/(1-q).
```

This controls both e_xx and e_xr. The undifferentiated principal equation
then gives `||e_rr||<=Z+(a+b0+3d/4)E+sup||r1||`. No affine G_x
approximation error is added: it cancels exactly in this error problem,
although G_x was correctly needed in the earlier W0 derivative estimate.

The individual initial C1 errors are independently reconstructed from the
retained interior actual-minus-nominal bounds, P1 velocity interpolation
error and PR100 smoothing allowances. They agree exactly with all three
producer inputs. PR100's paired higher-initial errors provide e0_rr and
v0_r separately; a nominal second-derivative magnitude cap is not substituted
for their actual errors.

## Exact numerical reconstruction and scope

`check_original_c2_audit.py` imports no producer code. It authenticates the
retained reports, reconstructs both differentiated principal estimates,
all initial errors, the complete feedback gate and the exact affine
sensitivity to the original scaled source derivative. It separately
expands a polynomial forcing before differentiating to test all product
terms. Seven deliberately omitted terms or wrong inverse choices fail the
checks. The affine audit has 36 named checks.

With the independently authenticated complete original-source derivative
report, `||F_x||_H<=905405545954284` in scaled x, the numerical
audit has 44 named checks. All three error majorants satisfy

```
sup||e_xx||_H <6.222e33,
sup||e_xr||_H <6.222e33,
sup||e_rr||_H <6.222e33.
```

Their size is dominated by the actual-minus-smoothed initial radial
second-derivative bound, not by the new source-derivative contribution.
The initial-free-wave term contributes more than 0.999999999 of the
differentiated error majorant's numerator. This is a share of an upper
majorant, not a lower bound on the unknown true error.
The normal and optimized numerical audit reports are byte-identical.

These are conditional continuous absolute error bounds in the declared
scaled-coordinate positive norm. The source/initial fields retain the
required qualitative physical-side regularity and open flat corner. Odd
reflected Dirichlet forcing may jump across the artificial interface;
there is no assertion of classical C2 regularity across the doubled shell.
No relative physical smallness, invariant observable, full natural
constraint propagation, nonlinear realization, cosmological prediction or
evidence of a higher-dimensional Big Bang follows from this gate.

**Assessment: PASS**, conditional on the separately reviewed actual full
source derivative certificate and the retained PDE/background premises.
