# Independent geometry and continuous-bound review

The proposed datum is well defined as a relative initial geometric coefficient
under the following explicit restrictions. The initial Cauchy section is a
physically identified part of the comparison, and must be transported if its
spacetime representative is changed. On a compact annular tube the initial
scalar has nonvanishing spatial differential, so its nearby level sets form a
regular compact foliation. Scalar labels and their reference normalization are
fixed from the retained background: the label `r` means the scalar value
`phi0(y(r))`, with the retained `ell=10^-6`; it is not an independently chosen
coordinate of the perturbed metric.

For a retained C2 family in the stated lower C4 initial topology, define
`A_e(y)` to be the induced three-area of the component of
`{phi_e=phi0(y)}` in that tube. The scalar implicit-function theorem and the
usual induced-area formula give the required C2 parameter dependence and the
label differentiability used below, after family-dependent shrinking of the
parameter interval. No exact SO(4) symmetry of the perturbed fields or their
level surfaces is assumed. This construction uses neither spacetime observer
orbits nor an asserted four-dimensional areal radius.

Compare two families with identical first metric and scalar tangent, and put
`D=(A''_new-A''_base)/A0`. The Hessian of the area functional on their common
first tangent cancels. At the background, a second scalar difference `w`
displaces the level by `-w/p`. A second angular metric coefficient `J` changes
the area density by `3J/2`; background radial displacement changes it by
`-3a w/p`. Thus

```
D = 3J/2 - 3a w/p.
```

For a general nonsymmetric second datum this formula uses the angular mean of
the indicated linear terms. On the restricted preparation difference they are
already radial. In the raw representative `J=h+2a zeta,w=0`; in the conditioned
representative `J=h,w=s=-p zeta`. Both give the same D. A pure radial second
gauge direction `(J,w)=(2a xi,p xi)` is annihilated. The equality concerns
second parameter derivatives, so it contains no extra factor of one half.

The claim must remain restricted to compact metric-mean preparation freedom
with the other canonical sectors fixed as in the published construction.
This single scalar datum is not claimed to control arbitrary second metric,
momentum or scalar preparations.

Write `chi_y=rho^2 h`, `z=chi/(ell rho^3)`, `b=ell rho a` and
`g=ell rho p^2/(3a)`. Substitution into the exact local preparation equation
gives

```
z_r = 3bz-h,
h_r = 3ell^2 z+(-b+g)h-2gD/3,
(z,h)(1)=(0,0).
```

Consequently `D=0` on the causal initial segment implies `h=chi=0` there by
uniqueness for this ordinary differential system. Data supported beyond the
segment need not vanish, and remain outside this local injectivity statement.

The fresh Fraction checker in `recompute.py` pins and reads the two old
background/response artifacts. It imports no producer and invokes no oracle.
The matrix infinity norm is at most `mu=1+3b_max`; its second row sum is
explicitly checked. A positive-term Taylor bound through degree six, with a
geometric bound for the remaining tail, proves `exp(mu/4)<4/3`. The integral
Gronwall bound therefore controls both z and h continuously, including the
flat endpoint. It requires no derivative bound on the ODE matrix.

The remaining free initial data are recovered through

```
beta=p/(3a),
s=beta(3h/2-D),
s_r=beta_r(3h/2-D)+beta(3h_r/2-D_r),
beta_r=-ell rho (p_y/(3a)-p a_y/(3a^2)),
V_L=h_r+b h.
```

All terms, including `beta_r`, enter the bound. With Gram matrix
`diag(1,2,1,3,1)`, the vectors `(h,0,0,h,s)`, `(h_r,0,0,h_r,s_r)` and
`(0,V_L,0,0,0)` are bounded respectively by `2|h|+|s|`,
`2|h_r|+|s_r|`, and `2|V_L|`. Their complete free-wave norm satisfies
`N/R <0.000429098849392<1/2000`, where
`R=max(sup|D|,sup|D_r|)` on `[1,5/4]` refers to the complete added preparation.

For `R<=2e33`, one obtains `N<=1e30`. The exact retained resolvent factor
`1/(1-q)` and the zero wedge imply
`|m(x)|<=(x-1)N/(1-q)`. Integrating against the actual proper-clock kernels
`k<=T-x` and `0<=k_T<=1`, including their distinct time factors, gives

```
Cv=17161/767711760000000000,
Ct=131/63975980000.
```

These bounds preserve strict negativity of both base coefficients. This is a
coefficient theorem conditional on the existing exact-family, C2 development,
and coupled uniqueness premises. It supplies no physical process selecting
the class, no uniform nonlinear epsilon radius, and no new proof of the
physical origin of the initial data. Each fixed admitted exact family retains
only the inherited qualitative, family-dependent small-parameter conclusion.
