# Propagating the refined initial curvature estimate

The candidate is the existing `W1_quadratic` (also called `W1_dagger`) and
the target is the same original formal second-coefficient solution. The
actual quadratic source `F`, shell data `G/M`, first-order fields, candidate
initial functions and candidate bulk iteration are unchanged. This note
changes one certified upper bound, not the functions or their true error.

The physical-side regularity and comparison argument is the one proved in
the original-source derivative checkpoint. Let `e=U_target-W1_quadratic`,
`T=3/2`, and use its positive full23 `H` norm. Write the retained actual
coefficient caps as `a,b,d,m` for the coefficients of `e_x,e_r,e` and the
Robin boundary term. The retained causal contraction bound is
`q=3409/8000000`. Define

```
c=a+b+3d/4,
comm=cap(A_x)+cap(B_x)+3cap(D_x)/4.
```

The original complete unsmoothed `F_x` bound, all absolute candidate
derivative bounds, the cancellation-preserving differentiated finite-bulk
defect bound `F1`, the C1 error cap `E`, and its reduced-defect cap `F0` are
retained. In particular, a bound on the *initial error* cannot replace the
absolute candidate initial curvature cap inside the differentiated
candidate-forcing ledger. The checker independently rebuilds that entire
ledger from the authenticated prior coefficient, source, candidate and
initial reports and reproduces each PR102 rational output exactly.

Only the initial displacement-curvature error cap `eurr` is replaced by
the finer complete-domain bound. The other four initial errors, denoted
`eu,eur,ev,evr`, remain fixed. For `v=e_x`, the initial velocity is bounded by

```
w0=eurr+a*ev+b*eur+d*eu,
J=max((ev+T*w0)/3, evr+w0).
```

The original candidate retains `r1(0,r)=0`; the generally nonzero initial
accelerations of its intermediate finite corrections remain in `F1`.
The same-source error shell source and its time derivative are zero, and
`M_x=0`. Applying the differentiated comparison theorem gives

```
cap(e_xx)=cap(e_xr)=(J+T*(F1+comm*E))/(1-q),
cap(e_rr)=cap(e_xx)+c*E+F0.
```

The new checker evaluates these formulas directly using exact rationals,
then compares their old-input specialization with the authenticated
original C2 report. It does not subtract rounded published bounds. Since
the new `eurr` is smaller and all coefficients in these inequalities are
nonnegative, this substitution preserves the theorem and tightens every
displayed second-derivative cap. The old and new upper-bound ratios compare
certificates; they do not claim the actual reconstruction has changed or
its true error has improved.

The supplied refinement report must identify the unchanged candidate, the
complete positive Gram norm, a *uniform upper bound* for the specified
actual-minus-smoothed initial `u_rr`, and a successful status. Controls
reject mutated lower-bound/point-value roles, a `NOT_PASS` status, changed
candidate identities, incomplete labels, an altered Gram matrix and a
negative upper bound. Formula controls also reject substituting an error
cap for an absolute candidate cap, omitting the initial velocity-gradient
error, dropping the commutator, the feedback inverse, or the original
source derivative.

These remain conditional continuous scaled-coordinate bounds on
`D+={0<=x<=3/2,max(0,x-5/4)<=r<=7/4-x}`, including physical-side one-sided
traces at the shell and outer boundary. They do not establish a C2
extension across the artificial reflected shell, relative physical
accuracy, complete subsidiary-constraint propagation, a finite-amplitude
nonlinear family, or a cosmological prediction. The optional
zero-`eurr` sensitivity table holds all other upper bounds fixed; it is
neither an attainable error nor a lower bound on the true error.
