# A conditional second-derivative error gate for the coupled response

This note gives a sufficient **input gate**, not a computed numerical C² error
for HDBLAST. It reuses the retained local response equation and the zero-wedge
contraction in `HDBLAST_TOTAL_RESPONSE_PREFLIGHT_20261023/analytic/LOCAL_TOTAL_REVIEW.md`.
It does not re-prove the regular-branch theorem, compute actual high-jet
source bounds, certify an Einstein residual, or establish a nonlinear solution.

## Assumptions

Let `D+ = {0<=x<=3/2, max(0,x-5/4)<=r<=7/4-x}` and let `H` be the retained
finite-dimensional Hilbert space of amplitudes. Set
`X(w)=max(sup_D|w|/3,sup_D|w_x|,sup_D|w_r|)`. Let `P_D,P_N` be the frozen
orthogonal Dirichlet/Robin projections and `M` the time-independent frozen
Robin map. Assume two sufficiently smooth functions `u, uhat` satisfy the
**same** coefficient matrices A,B,D0 and their mixed boundary equations on D+.
Write `e=u-uhat`, `f=F-Fhat`, `g=G-Ghat`, and assume

```
e_xx-e_rr=A e_x+B e_r+D0 e+f,
P_D e|r=0=0,  P_N e_r|r=0=M P_N e|r=0+g.
```

All derivatives, boundary traces, and one-sided corner values below exist and
are continuous. Assume `e`, `e_x`, the forcing and boundary inputs are in the
appropriate closed zero-wedge class (in particular they vanish for `x+r<1`),
including the compatible initial shell collar. The note applies either to two
smooth exact solutions with different data/sources, or to a smooth candidate
whose pointwise bulk and boundary defects are known and assigned to f,g with
consistent sign. A merely C¹ or piecewise-linear candidate does **not** by
itself meet these hypotheses: its classical second/third derivatives and
derivative of its pointwise defect may not exist. Such candidates need a
separate weak/interface treatment or a smooth reconstruction.

Suppose a valid inherited estimate gives `X(e)<=E`. In the same operator norms
assume, throughout the retained domain,

```
||A||<=a0,  ||B||<=b0,  ||D0||<=d0,
||A_x||<=a1, ||B_x||<=b1, ||(D0)_x||<=d1,
sup|f|<=F0,  sup|f_x|<=F1,
sup_shell|g_x|<=G1.
```

The old C¹ error estimate alone does not provide `F1`, `G1`, coefficient
derivative bounds, or the higher initial jets below. If approximate coefficients
are used, their discrepancies and time derivatives must be incorporated into
f and f_x; one may not silently treat the two operators as equal. If M varies
in time, the differentiated boundary has an additional `M_x P_N e` term.

## Reduction to the same contraction

Let `v=e_x`. Classical differentiation yields

```
v_xx-v_rr=A v_x+B v_r+D0 v + f_x+A_x e_x+B_x e_r+(D0)_x e,
P_D v|r=0=0,  P_N v_r|r=0=M P_N v|r=0+g_x.
```

Consequently, using the inherited zero-wedge bound `sup|e|<=3E/4`, the new inhomogeneous
bulk source has norm at most `F1+(a1+b1+(3/4)d1)E`. The existing finite-domain Duhamel estimates are
`X(Kh)<=(3/2)sup|h|` and `X(Lk)<=sup|k|`, and the same linear feedback
operator satisfies `||T||<=q_w=3409/8000000` on the zero-wedge subspace.
The derivative field and its initial/free wave must belong to that subspace.
The scalar inequality for the fixed point gives

```
Z = [J + (3/2)(F1+(a1+b1+(3/4)d1)E) + G1] / (1-q_w),
X(v)<=Z.                                                     (1)
```

Here J is any rigorous X bound for v's reflected free wave. It is not the
old free-wave gate for e: it requires higher initial data. Let
`e0=e(0,r)`, `v0=e_x(0,r)`, and on the relevant reflected initial footprint
`0<=r<=7/4` define

```
w0=v_x(0,r)=e0''+A(0,r)v0+B(0,r)e0'+D0(0,r)e0+f(0,r).
J=max[(sup|v0|+(3/2)sup|w0|)/3, sup|v0'|+sup|w0|].       (2)
```

Equation (2) follows directly from differentiated d'Alembert endpoints and
mixed odd/even isometric reflection, assuming compatible smooth zero collar.
The radial-velocity derivative and acceleration, including `e0''` and `f(0)`,
must be bounded; no radial derivative of `w0` is needed for J. If the first
initial data are exact, these differences may vanish, but that must be checked.

## Recovering the full coordinate second jet

By definition, `e_xx=v_x`, `e_xr=v_r`. The undifferentiated equation gives

```
sup|e_xx| <= Z,
sup|e_xr| <= Z,
sup|e_rr| <= Z+(a0+b0+(3/4)d0)E+F0.                           (3)
```

Thus (1)–(3) provide an explicit C² gate once every listed input has a
validated bound. The actual HDBLAST numerical values of F1, G1, a1,b1,d1 and
J have not been supplied by this note. In particular, the old large E and the
required data/source derivative errors can make the gate too loose for a
physical calculation. A small C¹ bound on e alone gives none of them.

For a scalar modal coefficient, with inherited `x=t/ell`,
`partial_y=-(ell rho)^(-1)partial_r`,

```
e_tt=ell^(-2)e_xx,
e_ty=-(ell² rho)^(-1)e_xr,
e_yy=(ell² rho²)^(-1)e_rr+rho_y/(ell rho²)e_r.
```

The `rho_y` term is essential. Tensor frame/covariant derivatives also need
the inherited connection and angular conversions; (3) alone is not a
pointwise Einstein-equation residual bound.

## Check and precise scope

`check_gate.py` uses exact rational polynomial algebra for a smooth-enough
zero-wedge synthetic scalar system with nonzero time-dependent A,B,D0 and
Robin M. It differentiates the bulk and boundary equations, checks the
nonzero coefficient derivative terms, verifies the reflected free-wave input
gate and rational inequality (1), verifies odd Dirichlet/even Robin-base
reflection, and rejects omitted coefficient-derivative and initial-forcing
terms and swapped reflection parity. It is an internal toy consistency check, not an HDBLAST
numerical source evaluation or an independent PDE existence proof. Run it in
normal and optimized Python. The analytic result depends on the inherited
mixed reflected contraction and classical corner assumptions stated above.
