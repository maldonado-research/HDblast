# Acceptance of the piecewise precontact support correction

14 October 2026. Independent internal mathematical review; not external
peer review. Reviewed the updated `PRECONTACT_SECOND_ORDER_DYNAMICS.md`
and `verify_precontact_dynamics.py` in `continuation-math-20261014`.

The correction is accepted. On the causal support enclosure at
`t=|T|`,

```
max(0,R_-−t) ≤ R ≤ R_++t,
s=R²−t².
```

The lower endpoint is therefore

```
min_s(t) = R_-²−2R_-t   for 0≤t≤R_-,
           −t²         for t≥R_-.
```

Both branches decrease, and they agree at `t=R_-`. The upper
endpoint `R_+²+2R_+t` increases. Thus the entire slab `|T|≤T_*`
has the conservative support range

```
[max(0,R_-−T_*)²−T_*², R_+²+2R_+T_*].
```

Requiring this whole interval to lie in the controlled smooth background
range repairs the global statement. At `R_-=1,T_*=2`, the correct
lower bound is `−4`; extrapolating the pre-origin formula incorrectly
gives `−3`. The verifier now includes the exact piecewise identity and
this counterexample as a wrong-formula control.

The separate safe positive-`s` restriction `T_*<R_-/2` remains
sufficient for a strictly positive lower endpoint. It does not replace
the requirement that the upper and lower endpoints stay inside the
inherited controlled interval.

The corrected geometry does not alter the finite-order wave-equation
induction, the compact-support initial-data reflection theorem, the
Gaussian-time acceleration formulas, or their limited scope. No
shell-contact or exact nonlinear existence claim is added.
