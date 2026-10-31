# Completed exact Cauchy packet

The new name is `W1_exact_Cauchy_packet`. All older candidate functions,
published formulas and certificates remain unchanged. Its displacement is
exactly the componentwise completed displacement of `W1_exact_packet` in
PR109. Its velocity on `0<=r<=R=7/4` is the full 23-component nominal
analytic output of the authenticated `Context._registered`: the same nominal
background coefficients, exact `E=exp(-ell*r)`, exact flat packet
`q=exp(4-1/((r-1)(2-r)))` for `r>1`, and `q=0` for `r<=1`.
The five mean assignments are included. There is no velocity interpolation
or convolution in this new candidate.

For each *completed velocity component* V_i, define its exterior by

    V_i(R+h) = V_i(R-) + h V_i,r(R-),  h>=0.

The coefficients here are the exact analytic traces, not the endpoints or
midpoints of their enclosures. Consequently the exterior value and first
radial derivative equal the exact interior traces component by component.
Extend across the origin with the original authenticated parity. The parity
matrix has diagonal entries +/-1 and is an isometry for the diagonal full23
Gram norm. Its differentiated action also preserves that norm.

The packet and all its derivatives vanish to every order at r=1: the
exponential of a negative reciprocal dominates each finite pole from
differentiation. Every completed component is built from a finite number
of these derivatives and smooth background functions with uniformly
separated rho and p denominators, which are explicitly enclosed on [0,R].
Thus all completed velocity derivatives needed here vanish at r=1; velocity
is identically zero on [0,1] and on [-1,1] after parity. There is no origin
corner. The former cutoff/interpolation nodes do not introduce joins into
the new analytic formula. The global completed velocity is C1 across +/-R;
no C2 velocity claim is made there. Together with the unchanged globally
C2 completed displacement, these are the initial regularity orders required
by the C2 finite response construction. The flat collar supplies the
homogeneous corner jets at x=r=0. This proof does not change the forcing,
boundary operator or response formula.

The derivative inventory for the velocity value is:

| Components | Highest q derivative in value | Highest in radial derivative |
| --- | ---: | ---: |
| K:0, K:2, K:4 | 1 | 2 |
| L:0, L:2, L:4 | 3 | 4 |
| N:0, N:2, N:4 | 1 | 2 |
| J:0, J:2, J:4 | 1 | 2 |
| w:0, w:2, w:4 | 2 | 3 |
| u:2, u:4 | 2 | 3 |
| v:2, v:4, H:2, H:4, T:+, T:- | exact zero | exact zero |

This is also visible directly from registered formulas: the y0/L velocity
contains f_y, where f is quadratic in the seed and its second derivative,
so its value needs q through order3 and its radial derivative through4.
The authenticated exact-packet adapter supplies q through5. Physical
radial differentiation multiplies normalized jet coefficient1 by ell.

All complete-domain error bounds are paired actual-minus-nominal bounds.
The exact q and E have identically zero functional differences, including
their derivative jets. The background coefficient-difference and actual
ODE remainder channels are retained separately. The ODE remainder bounds
for each derivative are inherited as independently validated derivative
remainder bounds, never obtained by differentiating an unspecified interval
remainder. Complete absolute candidate bounds are distinct from target
error bounds. The seed cover is [0,1], [1,65/64], then successive intervals
of width1/64 through R; the leaf cover is checked exactly. Relative
quotients at the flat support edge are undefined; absolute bounds apply.

The new local calculation concerns v_r at (x,r)=(0,65/64), not the older
u_rr point certificate. Let r*=65/64 and delta=1/2048. The old complete
P1 velocity has constant adjacent slopes s_- and s_+ on intervals wider
than delta on either side of r*. Differentiating its convolution with the
even normalized compact kernel gives

    old_v_r(r*) = s_- integral_0^delta k(s) ds
                 + s_+ integral_-delta^0 k(s) ds
               = (s_- + s_+)/2.

Continuity of P1 means no delta mass enters its first weak derivative.
The exact all23 rational profile supplies these two slopes. A fresh direct
registered actual analytic reconstruction gives each target interval A_i.
For a real interval B, let dist0(B)=0 if it contains0, otherwise the
smaller absolute endpoint. Then, with the exact Gram weights g_i,

    D^2 <= sum_i g_i dist0(A_i)^2,
    L^2 <= sum_i g_i dist0(A_i-(s_-,i+s_+,i)/2)^2,

where D and L are outward rational lower square roots. The new paired
errors supply an outward H-norm upper U. The fresh predicates are D>0,
U<=D/10, and U<L/2. The independent evaluator separately reconstructs
the same target and old P1 derivative. No old v_r point table or comparison
of two unrelated upper certificates is used.

These are conditional formal coefficient Cauchy-data statements. They
do not establish a nonlinear shell-contact solution, a natural residual,
an invariant physical observable, an external specialist assessment,
or a novelty determination.
