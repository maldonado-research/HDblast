# Direct preservation of the global C1 upper certificate

Reviewed the frozen comparison report
`4d0b83994270f6e2d3aa5433887a0d4b978acbac372bf28197050345de2881c3`
and its implementation against the original initial-data producer and the
PR95–97 direct error proofs. The sharper direct reproof is valid. It does
not require global boundedness of the quadratic exterior value or slope.

## Why the same smoothing majorants apply

The original `nominal_H_caps` are produced from continuous nominal ranges
on the physical cells of [0,R] in `full23_initial.py`, not from the affine
extension outside that interval. The old global free-wave size proof in
`check_robin_direct_error.py`, `check_refined_response.py`, and
`check_bulk_update_bound.py` uses these interior nominal value/derivative
caps plus the uniform smoothing remainders

```
alpha=M2_H delta^2/18,       beta=35 M2_H delta/128.
```

For each component the original interior nominal derivative is globally
Lipschitz on [0,R], including the cutoff joins, with the certified M2_i.
The new quadratic extension is C1 at R and has constant curvature kappa_i
thereafter. The verified inequality `|kappa_i|<=M2_i` therefore preserves
the same global Lipschitz bound for the derivative. Reflection preserves
the bounds at the negative endpoint; the exact initial zero collar avoids
any join issue at the origin. The checker verifies the vector comparison
`sum H_i M2_i^2<=M2_H^2`.

The even convolution's zero first moment gives the value remainder alpha
by the C1,1 Taylor estimate. The first absolute kernel moment gives the
derivative remainder beta. This proof needs bounded second derivatives and
the first derivative's Lipschitz continuity. It never takes a supremum of
the unbounded exterior quadratic value or slope. On the retained cone all
unsmoothed free-wave argument centers lie in [-R,R], where the original
nominal value and derivative are unchanged.

The actual-minus-interior-nominal initial bounds are unchanged, as are the
P1 velocity and its smoothing/error inputs. Consequently the same numerical
initial mismatch and smoothed free-wave size bounds remain valid for the
new candidate.

## Boundary history and the bulk residual majorant

The old global boundary correction uses shell history only through x=5/4.
The changed initial data begin at R-delta, beyond that history's entire
causal initial footprint, including the kernel width. Thus the shell traces
of the free wave, frozen first feedback and prior finite candidate are
unchanged exactly. The same actual M,G therefore give the same global
boundary-correction majorant.

The global first-feedback estimate is the norm estimate `q0*Bs_cap`; it
does not use an affine-extension-specific sign or exterior value bound.
It therefore applies unchanged to the modified free wave. The new
`W0_dagger-Bs_dagger` bound is the same first-feedback majorant plus the
unchanged boundary correction. Actual/frozen coefficient differences,
affine sources and operator norms are unchanged. Reconstructing the
coefficient-residual and feedback-residual terms therefore yields the same
numerical R bounding the new residual response. This is a new proof of the
same upper bound for a changed residual field, not an assertion that the
global residual field is unchanged.

The retained contraction q and boundary-inverted bulk factor beta also
remain unchanged. Thus `(initial_error+beta*R)/(1-q)` is exactly the same
rational upper certificate and is strictly below `3.541e25`. The direct
reduced forcing-defect majorant follows from the same re-proved R. The
weaker explicit candidate-comparison triangle certificate `<3.588e25`
remains valid independently; its additional charge is not needed in the
sharper direct route.

## Independent numerical audit

The accompanying audit JSON independently reconstructs every smoothing
term, the initial free-wave error, smooth free-wave size, unchanged
boundary correction, global feedback cap, both residual terms, R, and the
one-update C1 certificate using exact rational arithmetic from authenticated
earlier reports. All values agree exactly with the new report. All 23
curvature component caps lie inside the old M2_i bounds.

This preserves a numerical upper certificate for a new auxiliary candidate.
It does not prove equality of unknown actual errors, repair every higher
initial jet, supply the complete bulk derivative F_x, or close a quantitative
C2/Einstein residual theorem. No conceptual or arithmetic gap was found in
the direct C1 reproof.
