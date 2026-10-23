# A certified auxiliary response after shell contact

Ricardo Maldonado  
Independent researcher, Henderson, Nevada, US  
23 October 2026

**Review status: independent internal mathematical and numerical checks passed.** These checks were AI-assisted; this status does not represent external expert review. The result below is conditional on the retained smooth actual first fields, full linear response operator, source bounds, and separately certified source integral. The current exact certificate gives a strictly positive auxiliary response mean at one point after first shell contact:

\[
55{,}843{,}882
< \operatorname{mean}_{S^3} W_{yy}(3/2,1/4)
< 737{,}593{,}144{,}521.
\]

The mean uses probability-normalized sphere measure, and yy denotes the covariant normal-normal component in the registered parallel frame. The integer endpoints are rounded outward from exact rational bounds. This statement concerns the single scaled point `(x,r)=(3/2,1/4)`. Its boundary history extends to `x=5/4`, beyond the exact first-contact time `x=1`; no continuous neighborhood or region is certified here.

The quantity `W` is the auxiliary forced increment in the retained 23-amplitude formal second response. It has zero initial data and satisfies `W=S+TW`, where `S=KF+LG` contains the actual bulk and shell sources and `T` contains the complete linear bulk and Robin couplings. The registered common gauge and source conventions are part of the statement.

The source integral uses the selected n16 rational partition and continuous enclosures over complete bulk cells and boundary intervals. These retain uncertainty in the actual first solution. The reflected wave kernel has weight one on the overlap cone and one-half on the remaining direct cone for this Robin component. The boundary lift has the negative sign `Lg=-integral(g)`. Its normal source uses the one-copy outward shell Hessian QK, the Lorentz trace `-QK00+3QKab`, and the full varying-metric-lowered gauge Hessian QC_y. Copy doubling, time-index signs, and the gauge source therefore enter the evaluated integral explicitly.

The sharper feedback estimate follows from exact SO(4) averaging. The five invariant means are the yy, y0, 00 and scalar amplitudes, together with the coefficient of the sphere metric. The retained primitive matrices have zero coupling in both directions between these five amplitudes and their complement. This exact closure preserves every possible linear contribution to the mean, including connection, curvature, metric/scalar and Robin couplings inside the five-dimensional subsystem. Nonzero first-order harmonics still contribute to its quadratic sources before projection.

A positive matrix comparison bounds all five coupled amplitudes and all feedback iterations. It also uses the exact zero wedge `x+r<1` and the actual solution's continuous first derivatives to bound the boundary trace after contact. The total feedback allowance is less than `531,816,344.612908`. Combining it with the evaluated source interval gives a positive lower margin exceeding `55,843,882.796267`. All decisive comparisons use exact fractions; no numerical feedback field is required.

The [proof](../analytic-review/MEAN_BLOCK_POINT_CERTIFICATE_PROOF_FINAL.md) and [exact certificate](../analytic-review/MEAN_BLOCK_CERTIFICATE_NORMAL.json) record the operator and source bounds. Normal and optimized checker reports agree byte for byte. The checker authenticates its retained inputs and the n16 source reports; it does not rerun the inherited PDE or source proofs.

The retained total response U already has a conditional analytic construction and physical completion. This auxiliary W need not separately satisfy the complete sourced physical constraints. The result establishes neither a physical force nor an invariant backreaction observable, and supplies no full numerical waveform, nonlinear family, or finite-amplitude remainder estimate.
