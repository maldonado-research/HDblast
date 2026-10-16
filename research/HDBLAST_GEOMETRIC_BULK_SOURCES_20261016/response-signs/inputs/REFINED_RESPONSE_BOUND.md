# Refined bound for the selected formal second response

The retained full causal equation and contraction `q=2109/1000000` are unchanged.
The new interval calculation evaluates the same selected second initial data;
it does not choose different data to obtain a smaller bound. Its complete
cell enclosures include actual background uncertainty, flat packet endpoints,
the scalar compensator, lapse square and mixed shift-velocity conversion.

For initial parallel-frame metric component caps `g0,gr,gx` and scalar caps
`w0,wr,wx`, use

    U0_H <= 5 g0 + w0,
    U0_r_H <= 5 gr + wr,
    U0_x_H <= 5 gx + wx.

The positive Frobenius norm of a symmetric5x5 tensor is at most five times
its largest component; the angular measure has total mass one. The metric
and scalar triangle inequality is conservative. These are partial frame
component derivatives, with the required frame-connection conversion already
included in the initial certificate.

Odd/even reflection is an isometry on H. The mixed d'Alembert formula therefore
gives

    free_X <= max(U0_H/3 + U0_x_H, U0_r_H + U0_x_H).

Keeping the authenticated parent F_H and G_H bounds,

    b_X <= free_X + 3 F_H + G_H,
    ||U||_X <= b_X/(1-q) < 1.03e33.

`check_refined_majorant.py` evaluates every comparison with exact fractions,
authenticated by `PINS.json`. The earlier rounded bound was1.02e36; the ratio
of old to new rounded majorants is102000/103, approximately990.29. This is a
tightening of an upper bound, not a measurement or a reduction of the actual
field amplitude. The independent numerical review checks the norm join.

The separate angular closure proof puts this selected solution and its scaled
x/r derivatives in E_v. The independent sharp evaluation certificate gives
`C*=sqrt(321/19)<4.111`. Hence

    sup |U|_F < 1.270299e34,
    sup |U_x|_F, sup |U_r|_F < 4.23433e33.

These estimates do not control arbitrary higher derivatives or nonlinear
remainders. The universal full287-space factor17 is unchanged. Applying the
sharper factor to a numerical error requires its reconstructed candidate also
to lie in E_v.

For exact partial sums with N retained terms, the tail is bounded by
`factor*b_X*q^N/(1-q)`. Least counts, with both N and N−1 checked exactly, are:

| Absolute target | X norm | Selected pointwise value | Selected scaled first derivative |
| --- | ---: | ---: | ---: |
| 1e30 | 2 | 2 | 2 |
| 1e12 | 8 | 9 | 9 |
| 1 | 13 | 13 | 13 |
| 1e-6 | 15 | 15 | 15 |
| 1e-12 | 17 | 18 | 18 |

These counts concern exact terms only. They contain no quadrature, input,
modal-normalization, reconstruction or floating-point error certificate.
An independently bounded full residual already contains its truncation error;
do not add the same tail a second time.

The signed-profile consequence is separate: the authenticated initial scalar
mean has245 positive and245 negative complete cells. Six alternating signed
cells give five disjoint open zero brackets by continuity; two complete
witness cells have magnitude exceeding7e30. These belong to the selected
initial prescription at t=0. They are not evolved shell zeros, forces or
gauge-invariant observations. The chosen second initial jet is not unique.

Still required: evaluate actual common-gauge quadratic sources and their
approximation errors, complete the metric-source amplitude adapter, compute
a continuous numerical second response, and certify its complete residual.
