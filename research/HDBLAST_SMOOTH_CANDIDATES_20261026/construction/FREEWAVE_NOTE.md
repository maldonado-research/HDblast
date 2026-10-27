# One-point interval check of the tight candidate's free term

This is a single conditional numerical enclosure of the nominal free wave
B0_s, not the feedback B1_s or total candidate, actual Einstein-scalar field,
pointwise PDE residual, or a physical observable. The analytic candidate and
its global C1 transfer are specified separately.

Use x=1/1000, r=1499/1000, the location of the old candidate's kink, and the
`w:0` modal coefficient. Both endpoints `a=r+x=3/2`, `b=r-x=749/500` are
positive and strictly inside the same nominal analytic branch. Take
`u_s=k_{2^-30}*u0` and `v_s=k_{1/2048}*vP1`, the latter extended affinely
beyond its last node and by the frozen radial parity near the origin.

The selected free-wave coefficient and derivatives are

    B0 = (u_s(a)+u_s(b)+integral_b^a v_s)/2,
    B0_x = (u_s'(a)-u_s'(b)+v_s(a)+v_s(b))/2,
    B0_r = (u_s'(a)+u_s'(b)+v_s(a)-v_s(b))/2.

`check_smoothed_freewave.py` evaluates the nominal displacement and radial
derivative under the convolution at each endpoint with 32 weighted midpoint
panels. The normalized kernel panel masses are integrated exactly and are
positive rationals summing to one. The existing interval evaluator supplies
outward boxes for every sampled u0 and u0' value after authenticating all
seven source roots. The inherited global H-bounds U1 on u0' and U2 on u0''
then bound each component midpoint quadrature error by U1 δu/32 and U2 δu/32,
respectively. This component `w:0` has unit Gram weight. The sample intervals
are away from branch joins; the global midpoint error is conditional on the
cited inherited derivative caps.

Independently, the program reconstructs the `w:0` P1 velocity as a finite
sum of exact rational hinges and uses the closed polynomial kernel convolution
for both v_s and its normalized primitive. It then combines the two endpoint
intervals with the exact velocity terms using the displayed free-wave signs.
It re-authenticates upstream inputs after computing the result.

The complete rational intervals are emitted by the program; outward decimal
illustrations (not independent interval inputs) are approximately:

    B0   in [5.2659269239498705e28, 5.265934190685025e28],
    B0_x in [-5.489139622196741e27, -5.481497878599616e27],
    B0_r in [-5.285605543301422e31, -5.285604779127063e31].

These decimal displays use ordinary binary floating-point conversion and are
for illustration only; use the rational endpoints as the certified result.
The global C1 candidate bound remains far larger than observationally useful
scales. The values have no interpretation as detector signals or forces.

Run from any working directory with a complete published repository checkout:

    python ../checks/check_smoothed_freewave.py --repo-root /path/to/HDblast
    python -O ../checks/check_smoothed_freewave.py --repo-root /path/to/HDblast

The output must agree byte-for-byte. This program is more expensive than the
exact profile-ledger checker because it reauthenticates and evaluates the
inherited interval construction; it does not replay the original PDE solve.
