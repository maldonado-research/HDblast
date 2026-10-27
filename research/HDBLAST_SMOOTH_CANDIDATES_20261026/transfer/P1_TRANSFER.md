# A conservative C¹ error transfer to a candidate with C⁴ profiles

This note gives an exact arithmetic **C¹** error ledger for a particular
replacement of the frozen `V=B0+B1` approximation. It is conditional on the
prior full23 C¹ certificate, its operator norms and interpolation bounds.
It does not establish a numerical second-derivative error or Einstein residual.
The bound is very loose and is much larger than the prior bound.

## Profile smoothing

Let `p_u,p_v,p_du` be the authenticated continuous P1 displacement, velocity,
and separate radial-derivative feedback profiles on `[0,7/4]`. Let
`δ=1/2048` and use the even nonnegative compact kernel

    kδ(s) = 35/(32δ) [1-(s/δ)^2]^3  for |s|<=δ, and 0 otherwise.

It has integral one and is C². Extend the P1 profiles linearly beyond their
last node, and by their inherited odd/even parity on the negative half-line.
Their original zero collar about r=0 makes the reflected extension smooth
there. Let `u_s=kδ*p_u` and `v_s=kδ*p_v`; set `du_s=(u_s)'`, with opposite
radial parity. Thus the displacement and its derivative are now consistent.
The kernel convolution gives C⁴ profiles. The reflected free wave is built
from `u_s,v_s`; the first feedback uses these same profiles and `du_s`. The
registered exact rational displacement is no longer inserted in the free
wave. All reflected profiles remain zero for |r|<=1, since the original profiles vanish
through 65/64 and `δ<1/64`.

The minimum mesh width is 1/256, so distinct δ-bands are disjoint. On a band
around a knot, write the local P1 profile as `m_left z+Δm z_+` plus a constant.
Evenness and positivity imply

    sup |kδ*z_+ - z_+| = 35δ/256,
    sup |(kδ*z_+)' - (z_+)'| <= 1/2

away from the knot, with the same one-sided bounds at it. Therefore, in the
frozen positive Gram H norm, write

    Δu = (35δ/256) max_knots ||Δm_u||_H,
    Δv = (35δ/256) max_knots ||Δm_v||_H,
    Δdu = max_cells,endpoints ||m_u - p_du||_H
          + (1/2) max_knots ||Δm_u||_H.

Then `||u_s-p_u||∞<=Δu`, `||v_s-p_v||∞<=Δv`, and
`||u_s'-p_du||∞<=Δdu`. The cell endpoint maximum suffices because norm of
a vector affine in r is convex on each cell. All these norms use the exact
23-dimensional H diagonal and the stored node vectors. `../checks/check_transfer.py`
uses rational squared norms, outward integer square-root upper bounds, and
all 188 interior hinges; no floating-point decision is used.

## C¹ difference from the frozen candidate

Let `epsU,epsUr` be the previously authenticated P1 interpolation bounds
between `p_u,p_du` and the frozen exact rational initial displacement and
its radial derivative. The reflected free-wave formula and the half-domain
footprint give

    free = max((epsU+Δu+(3/2)Δv)/3, epsUr+Δdu+Δv).

For the first feedback, the frozen constant matrices have H-operator caps
`nA,nB,nD,nM`. Because reflection is an H isometry, its changed free field
satisfies

    sup |δB0| <= Δu+(3/2)Δv,
    sup |δ(B0)_x|, sup |δ(B0)_r| <= Δdu+Δv.

The shell trace is no larger than `Δu+Δv/4`: the history lies in x<=5/4,
and the difference is zero for x<1. The inherited Duhamel and boundary lift
bounds give

    feedback = (3/2)[(nA+nB)(Δdu+Δv)
                            +nD(Δu+(3/2)Δv)] + nM(Δu+Δv/4).

The triangle inequality yields, for exact total formal response U and the
new smooth candidate V_s,

    X(U-V_s) <= E_old + free + feedback < 9.711e31.                (1)

Here `E_old` is the prior exact rational upper bound, approximately
`1.1787044291e28`. Numerically the displayed upper caps are roughly
`Δu=4.411e27`, `Δv=1.059e20`, `Δdu=6.613e31`,
`free=9.708e31`, `feedback=2.126e28`. They are upper bounds, not measured
errors or evidence that the new candidate is inaccurate by those amounts.
The increased bound reflects the loose inherited independent P1 derivative
and mesh-slope estimates. In particular (1) does **not** preserve the old
`10^29` target. Comparisons of signs, magnitudes or physical effects would
need a far tighter independently verified ledger.

`../checks/check_transfer.py --repo-root /path/to/HDblast` pins the prior profile,
initial report and response report, checks nonoverlapping bands, support,
labels and Gram weights, then calculates (1) with exact fractions. Run normal
and optimized Python; all validation uses explicit exceptions. The script
verifies the algebraic input ledger and bounds, not the upstream PDE or this
note's analytic smoothness proof. Even with C⁴ smoothed profiles, the actual
full response's quantitative C² error, source derivatives, subsidiary
residuals and a nonlinear spacetime remain unestablished.
