# Independent scalar TT feedback bound

This note supplies a scalar resolvent bound and an independent exact check of the fixed exterior-cone cell geometry. It evaluates no forcing oracle and grants no new background or source premise.

Let `T=1/4`, `R=5/4`, and `C={(x,r):0<=x<=T, |r-R|<=T-x}`. For continuous `f`, define the unreflected retarded wave map

`Kf(x,r)=(1/2) integral_0^x integral_(r-x+s)^(r+x-s) f(s,v) dv ds`.

The cone is backward closed by the triangle inequality. Its lower radial boundary is `r=1+x`, so neither affine shell forcing nor Robin feedback can enter it. The complete registered TT response on this cone is the scalar fixed point `u=KF+K(Au_x+Bu_r+Du)` provided the structural certificate proves the TT row and column decouple and both registered TT initial values and velocities vanish. The identical forcing and operator for the two orientations then imply identical responses by uniqueness.

For a scalar C1 function on C set

`||u||X=max(2 sup_C|u|/T, sup_C|u_x|, sup_C|u_r|)`.

Endpoint differentiation of K gives

`(Kf)_x=(1/2) integral_0^x [f(s,r+x-s)+f(s,r-x+s)] ds`,

`(Kf)_r=(1/2) integral_0^x [f(s,r+x-s)-f(s,r-x+s)] ds`.

Thus each derivative is bounded by `T sup|f|`; the value is bounded by `T^2 sup|f|/2`. Consequently `||Kf||X<=T sup|f|`. This uses only continuity of f and no derivatives of F, A, B, or D.

Suppose independently authenticated whole-cone bounds are `|A|<=a`, `|B|<=b`, `|D|<=d`, and `|F|<=f`. Then

`sup|Au_x+Bu_r+Du| <= sigma ||u||X`,

where `sigma=a+b+Td/2`. If `q=T sigma<1`, the contraction/resolvent estimate gives

`||u||X <= T f/(1-q)`.

At the apex the complete feedback tail obeys

`|u(T,R)-KF(T,R)| <= (T^2/2) sigma ||u||X <= T^2 f q/[2(1-q)] = E`.

Therefore any exact source-cell Duhamel enclosure `[L,U]` gives the full scalar response enclosure `[L-E,U+E]`. This sums every feedback iterate; it is not a first-feedback truncation. A looser but also safe choice is `sigma=a+b+Td`. The existing full-state norm with value weight 1/3 instead produces `sigma=a+b+3d`; that larger weight is unnecessary for this restricted scalar calculation.

The parent should provide the actual numerical caps and the structural zero-data/decoupling certificate. In particular a full-state probability-sphere norm cannot be reused as a unit TT coefficient cap without its tensor basis Gram factor. Direct maximum of the actual TT coefficient intervals is the preferred scalar forcing cap.

The geometry checker clips the triangular cone by each cell's four coordinate halfplanes and computes its area with the rational shoelace formula. This independently verifies 20 positive-area cells for N=4, 72 for N=8, exact total area 1/16, edge areas h^2/2, and interior areas h^2. Rectangular interiors are disjoint. For each fixed x-row the consecutive retained radial intervals cover the complete cone section, so the union is exact. Wrong edge area and dropped-cell controls fail the area condition.

The original TT D expression can be bounded directly using independent interval background quantities. An optional simplification needs a separately verified background constraint: with `b=-(a_y+a^2)`, `d=rho^-2-a^2`, and `a_y=-rho^-2-p^2/3`, its square bracket equals `8a^2-10rho^-2-2p^2/3+2h^2+4V/3`. If the registered constraint is exactly `V=p^2/2-6a^2+6rho^-2`, this becomes `-2rho^-2+2h^2`, leading to `D=-8ell^2/c^2`. This identity is presented as a candidate algebraic simplification, not an authenticated premise of this review.

The certified quantity remains a registered-gauge conditional metric coefficient. No invariant force, physical gravitational-wave observable, nonlinear family, or cosmological interpretation follows from this bound.
