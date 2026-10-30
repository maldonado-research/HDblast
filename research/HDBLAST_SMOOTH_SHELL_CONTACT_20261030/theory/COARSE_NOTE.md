# Independent shell enclosure by cone norms

At the active-shell point `(x,r)=(9/8,0)`, this checker encloses the same tighter
integral-defined smooth candidate `V=B0+B1` as the detailed calculation. It does
not import the new detailed feedback checker, its formulas, or its output.
Instead it encloses `B0` directly and bounds the whole `B1` by positive cone
integration and inherited Gram operator norms. Its exact rational report is
identical in normal and optimized Python execution.

## Independent arithmetic

The three old profile/report files and initial producer are pinned by SHA-256.
The old producer's `Context` authenticates its seven complete upstream roots
before scientific imports and again after the evaluation. The checker verifies
common labels, Gram weights, all frozen operator block norm bounds, and all
actual-minus-frozen shell coefficient error bounds.

Local nominal `u0,u0'` caps follow from the continuous P1 interpolation-error
certificates on every cell intersecting `[0,x+max(delta_u,delta_v)]`: maximum
absolute nodal values plus that cell's certified interpolation error. These are
caps on the expression-defined nominal functions, although the P1 interpolants
are only used to derive the caps. Positivity of the smoothing kernel transfers
the local caps to the smoothed profiles and derivatives. The P1 velocity cap
uses convexity of each affine segment directly.

At `x=9/8`, nominal endpoint interval evaluation plus the published `C1,1`
convolution remainders encloses `u_s(x),u_s'(x)`. The velocity primitive is
integrated as a continuous P1 function by exact trapezoids. Its smoothing error
is bounded by `max|v| E|S|`, with `E|S|=35 delta_v/128`; the velocity smoothing
error uses the maximum adjacent P1 slope times the same first absolute moment.
This is a distinct, wider enclosure from the detailed exact hinge convolution.
The even free trace and its time derivative are

`B0_N=u_s(x)+V1_s(x)`, `B0_N,x=u_s'(x)+v_s(x)`.

Odd free value/time derivative vanish, and its radial derivative is the last
expression instead. Even free radial derivative vanishes.

## Physical reflection and the cone bound

The complete physical source is `f=A B0_x+B B0_r+D B0` on `z>=0`.
Only after the coupled matrix action is its output extended with component
parity. This reflection is an isometry of the diagonal Gram norm. No matrix
commutation with parity is assumed.

The local source bound is

`Fcap=(||A||+||B||)(U1+Vcap)+||D||(U0+x Vcap)`.

The initial profiles and the candidate are zero on `x+r<=1`. At the shell,
the physical Duhamel triangle has `s,z>=0`, `s+z<=x`; the active portion has
`s+z>1` and area `(x^2-1)/2`. Therefore

`||Kf|| <= (x^2-1) Fcap/2`,
`||(Kf)_x||,||(Kf)_r|| <= x Fcap`.

The area is **not** `(x-1)^2/2`: a source on the outgoing cone edge can remain
active for the full time range. The implemented control rejects that smaller
area. The Robin trace has magnitude at most `U0+(x-1)Vcap`, and vanishes for
`x<=1`. Consequently its lift is bounded by `(x-1)||M0||` times this trace cap;
its derivatives are bounded by `||M0||` times the trace cap.

This gives, approximately,

- bulk source cap `6.921758e28` in the Gram norm;
- full first-feedback value cap `9.550688e27`;
- full first-feedback derivative cap `8.073161e28`.

These are arithmetic upper bounds, rather than the values of the feedback.

## Conditional actual response and the shell equation

The published tighter-candidate transfer is independently recomputed with
exact outward arithmetic and satisfies `E<2.184e28` in the inherited C1 norm.
On the shell, the error vanishes through `x=1`, so the value error at `9/8` is
at most `E/8`. The primary scalar-mean value and time derivative are positive
even under this wider independent enclosure:

`8.278352e30 < U_w:0 < 8.302914e30`,
`1.973665e32 < U_w:0,x < 1.975718e32`.

The generic derivative error alone does not fix the sign of `U_w:0,r`.
The exact actual Robin law supplies the sharper relation

`PN U_r = Mactual PN U + G`.

The candidate uses the frozen `M0`, so its own boundary identity is
`PN V_r=M0 PN B0`; the actual law uses a different enclosed coefficient.
Writing `deltaM >= ||Mactual-M0||`, `nM0 >= ||M0||`, and
`tracecap >= ||B0||`, an independent enclosure centered on `M0 PN B0` has
Gram-norm radius

`deltaM tracecap + (nM0+deltaM)(B1cap+E/8)`.

The separate pinned affine-source checker supplies component intervals for
`G`. The added coefficient uncertainty is small but is included explicitly.
The resulting scalar-mean radial derivative interval is approximately

`2.150771690613386e27 < U_w:0,r < 2.157500049364436e27`.

All three signs refer to a registered formal scalar-mean coefficient or its
coordinate derivatives. The report gives intervals for all 23 slots. The
actual Robin reconstruction is only emitted for the 18 even/Robin rows.

The candidate's defect against the actual boundary law is

`PN V_r - Mactual PN V - G`
`= -M0 PN B1 - (Mactual-M0) PN V - G`.

This note establishes the candidate identity and the independent conditional
actual radial bounds. The detailed checker separately encloses the candidate
boundary defect. Neither construction makes the defect vanish.

## Reproduction and limits

Run `check_coarse_shell.py --repo-root PATH`; by default the pinned
`check_affine_shell_source.py` must be a sibling. The optional `--source-checker`
argument permits a separate local copy. Normal and `-O` execution use explicit
checks, with no bare assertions.

These results inherit the full23 formal-response, background, coefficient and
first-field premises. They do not provide an actual-field C2 error estimate,
pointwise Einstein residual, completed nonlinear spacetime, invariant physical
observable, or evidence of a higher-dimensional Big Bang origin. The method
comparison is an internal mathematical/code review, not an external specialist
assessment or a priority claim.
