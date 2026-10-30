# Absolute second coordinate jets of a causally filtered finite candidate

This calculation defines a new numerical candidate, retaining the PR100
quadratic-continued smoothed initial displacement, the same smoothed P1
velocity, actual coefficient functions, actual Robin matrix M, and actual
shell source G. The actual target problem continues to use its complete
original bulk source F. Only the auxiliary finite candidate uses the
one-sided time average F_epsilon, with epsilon=1/1024 in scaled time x.

Write `La=A partial_x+B partial_r+D`, `L0=A0 partial_x+B0 partial_r+D0`,
`P=(I-L M PN)^(-1)` and Bs for the PR100 free wave. The definitions are

```
f0=L0 Bs,
W0=P[Bs+K f0+L G],
f1=La W0+F_epsilon,
Z_epsilon=P[Bs+K f1+L G].
```

The complete physical coefficient/source vector is assembled before output
reflection, with the same 23-mode positive Gram H and the same parity law.
No coefficient is commuted through parity, K or the boundary inverse.
These are integral-defined finite updates, not a nonlinear solution or a
fixed point with La Z in its own forcing.

The physical domain is
`D+={0<=x<=3/2,max(0,x-5/4)<=r<=7/4-x}`. The norm is
`X(V)=max(sup||V||_H/3,sup||V_x||_H,sup||V_r||_H)`.
All fields in this argument vanish on `x+r<=1`, so
`sup||V||_H<=(3/4)X(V)`. The same causal inverse obeys
`||P||_X<=1/(1-b)`, `b=||M||_H/4<1`, and
`X(Kf)<=(3/2)sup||f||_H`, `X(Lg)<=sup||g||_H`.

## The differentiated principal-wave construction

For a prescribed forcing f with a bounded time derivative, let
`W=P[Bs+Kf+LG]`. Under the inherited smooth physical-side and flat-corner
premises, `Y=W_x` satisfies

```
Y_xx-Y_rr=f_x,
PN Y_r=M PN Y+G_x,  PD Y=0,
Y(0,r)=v_s(r),
Y_x(0,r)=u_s,rr(r)+f(0,r).
```

The initial acceleration includes the initial forcing. In the physical
reflected propagator one can differentiate in time distributionally:
`partial_x Kf=K f_x+Bfree[0,f(0)]`. This identity is the source of the
acceleration term; dropping it would produce the wrong initial problem.
Because M is time-independent, and G and the initial trace vanish in an
open corner collar, differentiating the causal boundary inverse introduces
neither M_x nor an initial shell delta. Thus

```
Y=P[Bfree[v_s,u_s,rr+f(0)]+L G_x+K f_x].
```

If `U2>=sup||u_s,rr||`, `V>=sup||v_s||`, `V1>=sup||v_s,r||`,
`f_init>=sup||f(0)||`, and `f_time>=sup||f_x||`, define

```
acceleration=U2+f_init,
J=max((V+(3/2)acceleration)/3,V1+acceleration),
Y_X<=[J+(3/2)f_time+Gx]/(1-b).
```

It follows that `||W_xx||,||W_xr||<=Y_X` and
`||W_rr||<=Y_X+sup||f||`. Only prescribed-forcing/boundary inversion is
needed here, not an inverse of a self-coupled differentiated bulk equation.

The Dirichlet output source can have a nonzero one-sided shell limit. Its
odd reflected extension can therefore jump at r=0. The conclusion concerns
continuous derivatives on the physical side, with one-sided derivatives at
the shell. It does not assert C2 regularity across the artificial doubled
shell. Bounded reflected forcing and its time derivative still give the
required physical-side first-jet estimate for Y; the principal equation
recovers the physical-side W_rr. Qualitative continuity of the physical
source remains an inherited premise; a magnitude cap alone would not prove
it.

## Authenticated initial and coefficient inputs

The separate absolute-data checker authenticates the exact raw profile
table, initial report and PR100 transfer report. It bounds the complete P1
velocity secants in the positive H norm before taking the maximum, and
uses convex convolution and the same final-slope continuation. For the
quadratic-continued displacement, the same nominal component curvature
caps hold globally before convolution. The resulting readable bounds are

```
||u_s||<2.48e31, ||u_s,r||<1.25e33, ||u_s,rr||<1.32e35,
||v_s||<2.71e23, ||v_s,r||<8.66e24.
```

These are absolute initial candidate jets, distinct from the initial error
bounds in PR100. The raw complete forcing cap is `FH=17387513156573`.
The quantity X(KF+LG) is not substituted for FH. The actual shell source
has `G_H<=487122`, and the independent PR91 full23 packing/Gram calculation
gives `Gx_H<58807538`. Its executable output is authenticated here.
The actual A_x,B_x,D_x,M_x ledger from PR100 is also authenticated.

## Applying the construction twice

Let a0,b0,d0 be frozen operator caps, a,b,d actual caps, and a1,b1,d1
actual time derivative caps. Let B be the authenticated Bs X bound.
For the free wave,

```
Q=U2+V1 >= sup max(||Bs_xx||,||Bs_xr||),
f0_H<=(a0+b0+3d0/4)B,
(f0_x)_H<=(a0+b0)Q+d0 B,
f0_init<=a0 V+b0 sup||u_s,r||+d0 sup||u_s||.
```

The previous absolute W0 bound is still valid for the PR100 quadratic
candidate: `W0_X<=B+(W0-Bs)_X`. Applying the differentiated principal
construction gives W0_x_X. For f1, the product rule gives

```
(f1)_x=A W0_xx+B W0_xr+D W0_x
        +A_x W0_x+B_x W0_r+D_x W0+(F_epsilon)_x,
f1_H<=(a+b+3d/4)W0_X+FH,
(f1_x)_H<=(a+b)W0_x_X+(a1+b1+3d1/4+d)W0_X
          +15 FH/(4 epsilon).
```

No derivative of actual F is used. Since `F_epsilon(0)=0`, the candidate
initial forcing is bounded by `a V+b sup||u_s,r||+d sup||u_s||`.
This is not the original target's initial acceleration. The target-error
forcing at x=0 is generally F(0,r), not zero. The open zero collar at the
physical initial/shell corner is preserved, so these distinct interior
accelerations do not introduce a new shell compatibility failure.

For epsilon=1/1024 the exact ledger gives these strict readable bounds:

| Candidate quantity | Uniform positive H-norm upper bound |
| --- | ---: |
| Z_epsilon | 9.37e32 |
| Z_epsilon,x and Z_epsilon,r | 1.250e33 |
| Z_epsilon,xx, Z_epsilon,xr and Z_epsilon,rr | 1.315e35 |

The rough absolute derivative scale is set mainly by the chosen initial
profile. These large coordinate numbers are not dimensionless physical
errors. No inference of small physical Einstein residuals is made from
their finiteness.

## Explicit change of candidate and scope

The source discrepancy is charged, rather than silently changed in the
target. With the same initial data/M/G,

```
Z_epsilon-W1_quadratic=P K(F_epsilon-F),
X(Z_epsilon-W1_quadratic)<=3 FH/(1-b),
X(U-Z_epsilon)<=E_PR100+3 FH/(1-b)<3.541e25.
```

The last strict comparison is checked by exact rational arithmetic. The
source-averaging charge is about 5.22e13 and does not consume the slack in
the displayed PR100 readable cap. Previously reported selected-point
values are not asserted unchanged; they need this perturbation charge or
a new direct evaluation. Reducing epsilon increases the kernel derivative
cost, and the magnitude-only source discrepancy supplies no convergence
rate at the initial face.

The result closes an **absolute physical-side second-coordinate-jet bound
for a distinct numerical candidate**. It does not close C2 accuracy to the
unsmoothed target: that comparison still needs actual F_x or an equivalent
proved estimate. A natural-equation residual requires its own geometric
conversion and source comparison. Tangential constraint propagation, a
full second-order compatible solution, a finite-amplitude nonlinear family
and cosmological or observational conclusions remain unproved.
