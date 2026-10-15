# Independent coordinate boundary-source audit

This is an internal exact-arithmetic audit of the new coordinate-source
checkpoint. It imports no producer, planning, prior scientific, or other
reviewer program. The actual-H and scalar whole-cell certificates and the
analytic shell box remain retained conditional inputs. Authenticating their
bytes and using their intervals does not independently reproduce their PDE
error proof. The calculation uses the fixed first longitudinal representative
and the background proper-clock slab. It supplies no evaluated second solution
`U2`, invariant response sign, nonlinear remainder, or external peer review.

## Finite audit groups

The standalone program has five finite groups: input identity and continuous
coverage; signed Volterra primitives and contraction; clock/radius and
derivative reconstruction; complete coefficient and harmonic rows; and source
norms. All decisions use reduced `fractions.Fraction` values and explicit
exceptions. Python optimization does not remove a decision. Input paths are
relative to the selected dataset root or the checker's enclosing capsule.

## Boundary reconstruction and residual

The pinned physical reconstruction has
`u_b=rho_b^-2 cosh(t)^(-3/2) w_b(t)` and `v_b=S(t)` because the normalized
radial factor is exactly one at the shell. The dynamic boundary row and the
actual trace identity are

\[
-w_z+\mu w_{tt}+\beta(t)w_b=0,\qquad
H=\ell(w_z-\gamma w_b),
\quad
\beta=\gamma+\mu\left(\frac{35}{4}\operatorname{sech}^2t-\frac{25}{4}\right).
\]

Cancellation of `gamma`, followed by `t=ell*x` and `y(x)=w_b(ell*x)`, gives

\[
y''+Ay=(\ell/\mu)H,\quad
A=\ell^2\left(\frac52-\frac{35}{4}\tanh^2(\ell x)\right),
\quad y(0)=y'(0)=0.
\]

Thus `y=F-K_A y`, with the positive triangular kernel
`F=(ell/mu) integral_0^x (x-s)H(s) ds`. On `[0,3]`,
`M=5ell²/2` bounds `|A|` and `q=ML²/2=9/800000000000<1`.

For each endpoint of the actual-H enclosure, completed-cell prefixes are
`p0=sum H_j(b_j-a_j)` and `p1=sum H_j(b_j²-a_j²)/2`. On the current
cell, the integral endpoint polynomial is

\[
I_\pm(x)=p_{0,\pm}x-p_{1,\pm}
              +\frac12H_{i,\pm}(x-a_i)^2.
\]

Its exact range uses both endpoints and its vertex when the vertex lies
strictly inside the closed interval. Multiplication by the positive interval
`[ell/mu_max,ell/mu_min]` uses all four endpoint products, including when
the primitive is negative. All 4,353 inherited closed cells are used.

Let `Fsup` be the largest absolute endpoint of these signed whole-cell
enclosures. For the approximation `p=F`,
`r=F-p-K_Ap=-K_AF`, so `||r||<=q Fsup`. The a posteriori contraction
bound is

\[
B_y=\frac{F_{\rm sup}}{1-q},\qquad
E=\frac{qF_{\rm sup}}{1-q},\qquad
|y-F|\le E,\quad |y'-F'|\le ML B_y.
\]

The last bound follows by differentiating the Volterra integral, without
differentiating a waveform cell. This yields independently
`Fsup≈1.8888428293538235e-4`, `E≈2.124948183046957e-15`,
`Eprime≈1.4166321220313047e-15`, and
`||u_b||_infinity<=1.1189581281054103e-8`.
The coarser absolute-H convolution independently recovers
`8.643987698350998e-8` and its previously retained exact fraction.

## Clock, radius, and first derivatives

For a cell `[a,b]`, use
`C=cosh(ell*x)^(-3/2) in [1-(ell*b)^2,1]` and
`tanh(ell*x) in [0,ell*b]`. The lower clock bound follows from
`cosh t <= (1-t²/2)^-1` and convexity of `(1-z)^(3/2)`.
The shell JSON's `H_squared` is the positive inverse-radius box `rho_b^-2`,
which is distinct from the waveform square. Then

\[
u=\rho_b^{-2}Cy,\qquad
\partial_\tau u=\frac{\rho_b^{-2}C}{\rho_b\ell}
 \left(y'-\frac32\ell\tanh(\ell x)y\right),
\quad\tau=\rho_b\ell x.
\]

The first shell rows give exact radial reconstruction
`u_y=-2a u-pv/3` and `v_y=-2p u-sigma2*v/2`. Marginal interval
products discard correlations but enclose all admissible shared parameters.

## Complete sources and sign witnesses

The pinned scalar boundary note gives one-copy coordinate coefficients

\[
S_p=-4pu^2+2\sigma_2uv-2v^2,\qquad
S_K=\frac{\sigma_2}{6}v^2-12au^2,
\]

and the degree-four Robin reduction gives

\[
E_N=-16au^2-\frac{16}{3}puv+\frac23\sigma_2v^2.
\]

These are coefficients in the stated representative, not independent
observables. Every square is evaluated with zero as its lower bound if its
interval contains zero. No endpoint-only squaring is used.

On the entire retained central bracket `[1.49998,1.49999]`, independent
quadratic optimization and the new residual give
`u in [4.994210294811386e-9,5.529396634158031e-9]`.
The retained opposite-sign actual-H witnesses and continuity supply a scalar
zero inside the bracket. At such a zero, `v=0`, and hence

\[
S_p\in[1.2405064085171157,1.6769780335004868]\,10^{-20}>0,
\qquad
S_K\in[-4.086332462113927,-3.333587487171056]\,10^{-17}<0.
\]

The full inherited tail `[2,3]` instead gives
`S_p<=-2.062886649452619e-7` and `S_K>=6.87484485838952e-8`.
These establish both signs for each continuous coordinate coefficient.
They imply no count or uniqueness of source zeros and no sign for the
second solution. Multiplication by `Y²` retains angular zeros.

## Harmonics and norms

With `Y=n1²-n2²`, `J=n1²+n2²`, `Z2=J-1/2`, and
`H4=Y²-J/2+1/12`, the exact unnormalized projection is
`Y²=1/6+Z2/2+H4`. The three scalar amplitudes are therefore the
coefficient divided by six, divided by two, and unchanged.

On the unit three-sphere, squared-coordinate moments are
`E[n1^(2a)n2^(2b)]=(1/2)_a(1/2)_b/(2)_(a+b)`.
Direct finite sums give

\[
E[Y^4]=\frac3{40},\qquad E[Z_2^2]=\frac1{12},\qquad
E[H_4^2]=\frac{19}{720}.
\]

Thus the complete `Y²` factor has squared angular norm `3pi²/20`,
whereas the isolated degree-four factor has `19pi²/360`.
These factors must not be interchanged for a degree-four-only row.

For the complete scalar coefficient sources,

\[
\|S Y^2\|^2_{L^2(d\mu_0)}
=\frac{3\pi^2}{20}\rho_b^4\ell
 \int_0^3\cosh^3(\ell x)S(x)^2\,dx.
\]

The coefficient intervals are rounded outward once to denominator `2^144`
for the final CSV. Norm quadrature squares those final intervals over every
whole cell. With `z=(ell*x)^2<=9e-12`, the clock measure has
`1<=cosh^3(ell*x)<=1+2z`. Indeed,

\[
(1+2z)(1-z/2)^3-1
=\frac z8(4-18z+11z^2-2z^3)\ge0
\]

on this window, while `cosh t<=(1-t²/2)^-1` follows termwise from the
positive series. Rational Machin bounds use 64 terms for `atan(1/5)` and
16 terms for `atan(1/239)`, paired with the next alternating partial sum,
and are then rounded outward to the common denominator.

Scalar-amplitude norms do not establish a tensor contraction convention,
an Einstein energy, or a nonlinear response estimate. Separate marginal
norms alone do not determine source cross terms.

The metric correction coefficients are explicitly evaluated before taking
norms: `dSp=-4pu²+2sigma2*uv`, `dSK=-12au²`, and
`dEN=-16au²-(16/3)puv`. They are compared with the corresponding scalar
terms evaluated under the same clock measure. The independently verified
upper squared-norm ratios are approximately `7.406446493642206e-25`,
`9.74640488921516e-50`, and `3.482728988064778e-33`. In particular
`||dSp Y²||/||-2v²Y²||<9e-13`. This integrated statement is compatible
with opposite coordinate signs at scalar zeros, where the scalar term
vanishes. No lower pointwise relative correction is inferred at those zeros.

## Degree-two algebraic check

This appendix does not add a degree-two numerical Robin row. The pinned
preprojection tangential trace-free bulk source is
`24u² Y HessY + (12u²-2v²) dY tensor dY` modulo the metric trace.
Use `2Y HessY=Hess(Y²)-2dY tensor dY`, the degree-two coefficient
`[Y²]_2=Z2/2`, and the pinned scalar tensor identity
`[(dY tensor dY)^TF]_2=(HessZ2)^TF/2`. The metric terms cancel,
leaving `-v²(HessZ2)^TF`. The linear trace-free row has factor one half,
so its algebraic forcing is `C2=-2v²`. Then the boundary projection is
`SK2=SK/2`, giving independently

\[
E_{N,2}=(C_2)_y-4S_{K,2}
=-4vv_y-2S_K
=24au^2+8puv+\frac53\sigma_2v^2.
\]

This does not supply the degree-zero master reduction or a coupled
second-order solution.

## Reproduction and deliberately incorrect calculations

Run the standalone checker normally and under optimization, then its
separate subprocess control runner:

```sh
python -I -B sol/check_coordinate_shell_sources.py --output sol/SOL_NORMAL.json
python -I -O -B sol/check_coordinate_shell_sources.py --output sol/SOL_OPTIMIZED.json
python -I -B sol/run_sol_controls.py
```

The normal checker independently rebuilds every one of the 113,178 final
scientific CSV endpoints and all exact JSON leaves, display leaves, source
hashes, and scope claims. Seven mathematical controls actually calculate a
wrong force, triangular factor, signed inverse-mu endpoint, resolvent,
radius power, angular moment, or clock measure, and require rejection for
the corresponding mathematical reason in both modes. Two additional controls
privately corrupt the last-cell `Sp_hi` endpoint and the final EN-material
norm endpoint. They must fail at those exact locations in both modes.
The runner never edits the producer or its computed artifacts.
