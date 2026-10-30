# Independent review of the curvature-matched repair

The reviewed repair defines a distinct auxiliary candidate with the same
nominal initial displacement on the retained footprint, the same convolution
width, the same P1 velocity and the same actual coefficients and fixed affine
sources. Outside `R=7/4`, it adds the quadratic Taylor term with coefficient
equal to the **exact nominal left curvature**. The pinned interval encloses
that expression-defined coefficient. No midpoint or independently chosen
interval representative is substituted for it. The old affine candidate is
preserved as a separate mathematical object.

## Endpoint identity and real error improvement

I independently verified the half-kernel moments and the three exact endpoint
changes: `kappa delta^2/36`, `35kappa delta/256`, and `kappa/2` for value,
first radial derivative and second radial derivative. The added half-curvature
removes the old affine continuation's half-curvature term algebraically.
Consequently the repaired actual endpoint error is the pinned paired
actual-minus-nominal curvature error plus the same half-kernel Lipschitz
remainder. This retains correlations algebraically rather than subtracting
two separately widened enormous curvature intervals.

All 23 new endpoint intervals in `REPAIR_NORMAL.json` were independently
recomputed. For the scalar mean, the old actual error magnitude is proved
strictly greater than `1.6245e33`, while the new actual error magnitude is
strictly below `1.820e27`. Exact cross multiplication gives

```
1.6245e33 > 800000 * 1.820e27.
```

Thus this is a genuine improvement of that one endpoint initial curvature
error, rather than a comparison of two upper certificates. The repaired
interval contains zero and does not prove a new nonzero error sign. The
independent actual ODE endpoint route overlaps every paired actual-curvature
interval and confirms the old nonzero error lower bound.

## C1 price and retained comparison bounds

I independently reconstructed the weighted H norm upper bound on the fixed
endpoint curvature vector and its input value/derivative charges. They match
the reported kernel factors. The free-wave X charge is the larger of the
value charge divided by three and the derivative charge; here it is the
derivative charge. This is only an initial/free-wave perturbation input.
Propagation through the W0/W1 constructions requires the separate shared
actual-Robin comparison formulas and their bounded operators.

All 23 endpoint curvature magnitudes also fit within the existing nominal
M2 component caps. The new quadratic continuation therefore has the same
C1,1 curvature bound on the whole extended displacement. The old C1
mollifier remainder inequalities remain available with those same constants
for a separate proof for the new candidate. This does not identify old and
new candidate functions or automatically transfer a frozen old certificate;
the full reconstruction/error identities must still be checked.

## Uniform initial curvature and causality

The new continuation matches the second derivative at R and has zero third
derivative outside. It removes only the exterior curvature jump. The two
inner cutoff joins are retained in the jump-aware convolution bound. I
independently reconstructed the resulting componentwise smoothing charges
and weighted norm square. Together with the same three continuous paired
channels they give the reported strict uniform **initial-data** upper
certificate `<6.219e33`. This remains a potentially loose upper bound, not
a claim that all initial errors improve pointwise.

The extension difference is supported beyond R before convolution. After
convolution its initial restriction is zero for `r<=R-delta`. Mixed reflected
finite propagation and the same Volterra boundary reconstruction therefore
give zero candidate difference when `x+r<R-delta`. Both the selected point
`(9/8,0)` and the retained shell segment `x<=5/4,r=0` lie strictly inside
that protected region. Their previously verified values are unchanged; no
new point value is inferred from a smallness estimate alone.

The repair is about an auxiliary numerical continuation. It does not remove
the remaining actual-minus-nominal error, the interior cutoff charges or the
complete differentiated-source obligation. It proves neither a global evolved
C2 error improvement nor a full Einstein/constraint residual, a nonlinear
family, an invariant physical quantity or evidence for a higher-dimensional
Big Bang. No mathematical or scope issue was found in the reviewed repair.

Exact recomputation details and the reviewed checker/report hashes are saved
in `REPAIR_REVIEW.json`.
