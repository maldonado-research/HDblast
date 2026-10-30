# Transferring the causal correction to the same target

The quadrature encloses the auxiliary correction `C=W1-W0` at `(9/8,0)`.
The previous point table encloses `W0`. Adding these boxes encloses `W1`;
discarding interval correlations can broaden the sum but cannot exclude it.

The inherited global result remains `X(U-W1)<E`, with the same readable
`E=3.541e25` used in the preceding local checkpoint. For each positive Gram
weight `H_i`, the zero wedge and the interval from shell time one to `9/8`
give the pointwise bounds

```
|(U-W1)_i| < (E/8)/sqrt(H_i),
|(U-W1)_{x,i}| < E/sqrt(H_i).
```

For a Robin row, `U` and `W1` have the same fixed actual `M` and `G`, hence
`PN(U-W1)_r=M PN(U-W1)`. Its component radius is `(E/8)` times the weighted
dual norm of that row of `M`. The candidate's radial box is also intersected
with the independently valid affine boundary law `M PN W1+G`. For Dirichlet
rows, the value and tangential derivative vanish exactly; the normal error
retains the generic first-derivative radius `E/sqrt(H_i)`.

The resulting target boxes are intersected with the earlier best target
boxes. Both contain the same conditional `U`, so this intersection preserves
the target. Exact overlap and containment are checked for all 69 entries.
Any narrower target boxes arise from local candidate information, with the
same global error theorem and the same underlying initial prescription.

This comparison distinguishes three quantities: the numerical integration
width enclosing `C`, the inherited error bound for `U-W1`, and their resulting
pointwise enclosure of `U`. A sign of `C` establishes neither a sign of
`U-W1` nor a physical force or observable. The inherited initial-data limiting
majorant is an upper bound; it does not give a lower bound on the true error.
