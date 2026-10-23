# Independent analytic audit of the shell-reaching point certificate

23 October 2026. **Accepted conditionally.** The exact point allowance
`E_trace = 14866325944167915/3992782` is valid for the retained auxiliary
forced increment at `(a,b)=(3/2,1/4)`. This audit independently reconstructs
the geometry, operator estimates, shell-source conversion and rational
arithmetic. It does not substitute for review of an evaluated source ledger.

The accepted smooth background, actual first master and gauge solutions,
their authenticated continuous error bounds, and the retained second-response
construction remain premises. Their expensive PDE and source certificates
have not been independently replayed here. In particular a source box is an
enclosure of one fixed actual solution, not a license to select unrelated
values at different points.

## 1. Dependency closure and exact half-line kernel

Let `d=a-b=5/4` and

```
D+ = {(s,u): 0<=s<=a, u>=0, |u-b|<=a-s},
D- = {(s,u): 0<=s<=d, 0<=u<=d-s}.
```

For a physical backward dependency `(t,v)` of `(s,u)` in D+,
`|v-b|<=|v-u|+|u-b|<=a-t`. Reflected dependencies also satisfy this bound:
their defining condition `u+v<=s-t` implies `|v-u|<=s-t`. Every boundary
time used at `(s,u)` is at most `s-u<=d`, and the backward cone of each
`(tau,0)`, `tau<=d`, is contained in D+. Thus all iterations and all Robin
traces stay in this domain. Its largest `s+u` is `a+b=7/4<8`; no outer
boundary data or extension is being supplied.

Splitting the original signed spatial integral at zero gives

```
Kf(a,b) = (1/2) integral_D+ f + (1/2) integral_D- Jf
         = (1/2) integral_(D+\D-) f + integral_D- P_N f,
J=P_N-P_D,                 Lg(a,b)=-integral_0^d g(s) ds.
```

The yy row is Robin, so its scalar weights are respectively `1/2` and `1`.
This does not justify even reflection of the other amplitudes. Independent
Cartesian polygons give

| Region | Exact physical `(s,u)` area |
|---|---:|
| D+ | 47/32 |
| D- | 25/32 |
| D+ minus D- | 11/16 |
| `s+u<1` in each unfused cone | 1/2 |

The active D+ and D- areas are `31/32` and `9/32`. Consequently the
remaining Robin kernel weight is `(31/32+9/32)/2=5/8`. Removing only one
of the zero triangles would be wrong. These are physical areas, so no
characteristic-coordinate Jacobian is additionally applied.

The sign of L follows without convention guessing: for `x>r`,
`partial_r[-integral_0^(x-r)g]=+g(x-r)`. On the shell this supplies exactly
the prescribed radial boundary derivative. No Robin exponential belongs in
L; Robin feedback is already in T.

## 2. Full-state contraction and point feedback bound

Use the retained norm
`||V||_X=max(sup|V|_H/3,sup|V_x|_H,sup|V_r|_H)` on D+.
Endpoint differentiation of K gives derivative bound `a sup|f|_H` and
value/3 bound `a^2 sup|f|_H/6`. Therefore `||K||<=3/2`. For L the
corresponding bounds are `sup|g|_H` and `(d/3)sup|g|_H`, so `||L||<=1`.
The reflection J is an isometry in H. No extra reflection or component
multiplicity is inserted.

The retained coefficient recipe bounds the complete tensor/scalar operator.
Recalculation of its displayed rational majorants gives

| Majorant | Exact value | Bound used |
|---|---:|---:|
| A* | 58816389/1000000000000 | 1/10000 |
| B* | 2207219/50000000000 | 1/10000 |
| D* | 299136118485884125939/1386750000000000000000000000 | 1/1000000 |
| M* | 40880253/100000000000 | 1/2000 |

Thus `sigma=1/10000+1/10000+3/1000000=203/1000000` bounds the full bulk
combination, and `m=1/2000` bounds the full Robin operator. These are H
operator bounds, not entrywise maxima. Their restriction to the selected
23 amplitudes retains all couplings, including the scalar derivative
coupling and both TT sectors. The retained matrices have five Dirichlet
amplitudes `L:0,L:2,L:4,v:2,v:4` and eighteen Robin amplitudes. The full
287-dimensional space instead has Dirichlet/Robin dimensions 69/218;
those dimensions must not be confused.

For `W=S+TW`, `S=KF+LG`, the actual-source caps are
`F_H=17387513156573`, `G_H=3851400`. Backward closure allows the same
equation and actual source caps to be restricted to D+. Hence

```
q_D = a sigma+3m = 3609/2000000 < 1,
S_D_cap = a F_H+G_H = 52162547172519/2,
B_D = S_D_cap/(1-q_D) = 52162547172519000000/1996391.
```

The retained exact zero wedge implies W and its first derivatives vanish
for `s+u<1`; continuity supplies the limiting values at the support edge.
It follows either by the retained invariant closed C1 subspace or directly
by causality and uniqueness. This is an exact premise, not a threshold
imposed on small numerical values. Consequently bulk feedback is zero
there, and the active bulk point allowance is `(5/8)sigma B_D`.

On the complete shell segment `[1,d]`, the actual W is C1,
`W(1,0)=0`, and `|W_x|_H<=B_D`. The fundamental theorem of calculus gives
`|W(s,0)|_H <= (s-1)B_D`. Therefore

```
|L(M P_N W)(a,b)|_H <= m integral_1^d (s-1)B_D ds
                     = m(d-1)^2 B_D/2.
|W(a,b)-S(a,b)|_H <= [(5/8)sigma+m(d-1)^2/2]B_D
                   = 14866325944167915/3992782
                   = 3723300181.2189884... .
```

No derivative of an actual source or numerical reconstruction of W is
needed for this estimate. The existing C1 actual solution is enough. A
subsequent evaluated feedback field or continuous residual would require
its own C1 error control. C0 forcing alone would not justify a claim of
classical C2 regularity, but this point certificate makes no such new claim.

The yy probability-sphere mean has H-dual norm one: its constant scalar
basis mode has Gram norm one and the other scalar modes have mean zero.
Thus an honest interval `[Slo,Shi]` implies exactly
`mean W_yy(a,b) in [Slo-E_trace,Shi+E_trace]`. Neither a factor three nor
the angular point-evaluation factor 17 is added to this value estimate.

## 3. Shell conventions and continuous source enclosures

The outward direction is `partial_y=-(ell rho)^(-1)partial_r`. With the
one-copy outward affine Hessian QK, the shell normal harmonic source is

```
GN=-2 tr_eta QK-2 QC_y,
G_yy=-ell rho GN=2 ell rho(-QK00+3 QKab+QC_y).
mean LG=-2 ell rho integral_1^(5/4)(-QK00_0+3 QKab_0+QC_y_0) ds.
```

The negative time trace is Lorentzian; the factor three is the tangent
dimension of S3. The four remaining angular tensor slots have zero trace.
An independent trace calculation from the shell adapter's raw scalars A,B
gives `[4A+B(D-3L+8M)]Y^2`, using `Delta Y=-8Y`; its mean is the bracket
divided by six. This agrees with `-QK00_0+3 QKab_0` from the projection.
The scalar and mixed shell sources do not directly occur in this mean
integral, but their feedback is retained in B_D.

Use `one_copy_outward_hessian`, never a doubled two-copy value or an already
converted `coordinate_physical_forcing` with the same conversion factors.
Use the full varying-metric-lowered `QC['y']`, not `QC_background`.
They coincide on the actual harmonic first solution; a control designed
to distinguish them must use off-harmonic manufactured jets. The parallel
frame Robin matrix has no extra coordinate `2a` tangential diagonal.

The inspected oracle `actual_jet(I(slo,shi),I(ulo,uhi),mode='actual')`
encloses derivatives of the actual first fields on the **whole** rectangle.
On `r=0`, `s in [1,5/4]`, it retains the global C3 errors and substitutes
only valid tangential derivatives of the first boundary identities.
Shell-eliminated identities must never be normally differentiated.
The exterior error improvement is permitted only when `ulo>=shi` for
the entire enclosing box, even if the clipped polygon itself lies outside.
All current point-cone boxes have ample margin inside `x+r<=8`.

Exact cell area times a whole-box interval is a valid integral enclosure;
exact segment length times a whole-interval shell enclosure has the same
property. This argument uses inclusion for every point, not sampled values.
It needs no derivative-based quadrature remainder. Shared source uncertainty
and within-cell variation already included in those intervals must not be
added a second time. Repeating the same input interval in an expression may
widen it; it does not invalidate inclusion. The source-map loader checks a
status field rather than a hash, so the independently pinned map hash is
required: `34b10b62d0b36a4506a35407a591a9fdff319971994f25e0b5f031f48eb4f615`.

## 4. Executable receipt and scope

`check_analytic_certificate.py` independently checks polygon areas,
all displayed contraction/feedback fractions, retained source-cap joins,
the selected boundary projectors, and the yy mean Gram normalization.
Eleven exact synthetic false-formula witnesses distinguish omitted
reflection, all-even Dirichlet reflection, a missing wave half, one-sided
zero-wedge removal, wrong L/radial signs, Euclidean time trace, dimension
four, omitted QC, doubled QK, and an extra factor three. These are analytic
controls, not claims that malformed producers were run or actual waveform
values sampled. Normal and optimized runs agree byte-for-byte in
`ANALYTIC_REVIEW_NORMAL.json` and `ANALYTIC_REVIEW_OPTIMIZED.json`.
Their `input_pins` record every read mathematical/source input.

No hidden regularity, copy, orientation, trace, or domain mismatch was found
under the stated retained premises. Acceptance of an evaluated numerical
interval additionally requires its independent geometry/source/ledger
review. An interval containing zero is valid progress and resolves no sign.
The retained conditional analytic completion of total U already exists;
this calculation does not produce a full evaluated numerical waveform.
W is an auxiliary reduced-operator summand and is not separately a sourced
physical-constraint solution. No invariant force, nonlinear family,
finite-amplitude validity, or external novelty conclusion follows.
