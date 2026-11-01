# Independent angular law and causal-sign fraction

This calculation uses only the round unit three-sphere and the quadratic
harmonic `Y(n)=n_0^2-n_1^2`. It imports no producer or first-oracle code.
It concerns the causal character of a scalar gradient, not a shell matter
energy or a finite-amplitude statement.

## Uniform measure and the gradient

Outside coordinate sets of measure zero, write

\[
n=(\sqrt z\cos\varphi,\sqrt z\sin\varphi,
     \sqrt{1-z}\cos\chi,\sqrt{1-z}\sin\chi),
\quad 0<z<1,\quad 0\leq\varphi,\chi<2\pi.
\]

The induced metric is diagonal with entries
`1/[4z(1-z)], z, 1-z`. Its determinant is `1/4`; hence the volume element
is `(1/2) dz dφ dχ`. Division by `vol(S^3)=2π²` shows that
`Z=n_0²+n_1²` is uniform on `[0,1]`, while `φ` and `χ` are independent
uniform circle angles. Thus `Y=Z cos(2φ)`.

The ambient gradient is `(2n_0,-2n_1,0,0)`, its squared length is `4Z`,
and its normal component is `n·grad_R4(Y)=2Y`. Orthogonal projection
to the sphere therefore gives the exact pointwise identity

\[
|D_{S^3}Y|^2=4(Z-Y^2).
\]

For reference, the normalized moments are

\[
\langle Y^2\rangle=\frac16,\qquad
\langle |DY|^2\rangle=\frac43.
\]

These follow either from the independent `Z,φ` law or from
`E[n_i⁴]=1/8` and `E[n_i² n_j²]=1/24` for `i≠j`.

## Strictly timelike angular fraction

Let `a=ψ_t` and `b=ψ/cosh(t)`, evaluated at a single event. Up to the
positive factor `2/ρ_b²`, the second derivative of the scalar-gradient
norm has sign equal to

\[
G=-a^2Y^2+b^2|DY|^2
  =4b^2Z-(a^2+4b^2)Z^2\cos^2(2\varphi).
\]

If `a≠0`, define

\[
r=\frac{4b^2}{a^2+4b^2}\in[0,1).
\]

For `Z>0`, strict negativity is equivalent to
`Z cos²(2φ)>r`. Conditional on `φ`, the uniform `Z` law gives probability
`max(0,1-r/cos²(2φ))`. Averaging over `φ` and using the circle symmetries,
with `α=arccos(√r)`, gives

\[
\begin{aligned}
F(r)&=\frac2\pi\int_0^\alpha(1-r\sec^2\theta)\,d\theta\\
&=\frac2\pi\left[\arccos\sqrt r-\sqrt{r(1-r)}\right].
\end{aligned}
\]

The fraction counts strict negativity. All zero sets below have measure
zero except when both coefficients vanish.

* `a≠0, b=0`: `r=0`, `G=-a²Y²`, and `F=1`; the locus `Y=0`
  still has `G=0`, so the sign is not negative at every angle.
* `a=0, b≠0`: `G=b²|DY|²≥0`; the negative fraction is zero.
  The same formula extends continuously to `r=1` and gives zero.
* `a=b=0`: `G=0` at every angle; the negative fraction is zero and
  `r` is undefined.
* If `a≠0` and `b≠0`, then `0<F(r)<1`. Both positive and negative
  values occur: for `n=(1,0,0,0)`, `G=-a²<0`; for
  `n=(1/√2,1/√2,0,0)`, `G=4b²>0`.

Thus no nonzero `b` permits an all-angle negative claim, even when the
negative fraction is arbitrarily close to one. Even at `b=0`, the zero
set prevents strict all-angle negativity. The coefficient does not,
without remainder control, supply a finite-amplitude sign statement.

## Exact rational tail bound

Put `x=√r`. Then

\[
1-F(r)=\frac2\pi[\arcsin x+x\sqrt{1-x^2}].
\]

The bracket vanishes at zero and has derivative `2√(1-x²)≤2` on
`[0,1]`. Therefore

\[
0\leq1-F(r)\leq\frac4\pi\sqrt r <\frac43\sqrt r
\quad(r>0),
\]

using only the classical exact inequality `π>3`. For any nonnegative
rational `δ` with `r≤δ²`, this yields the wholly rational lower bound
`F(r)≥1-4δ/3`. In particular, `r≤10⁻⁶` implies
`F(r)≥749/750`. Strictness can be retained if desired when `r>0`.

## Reproduction

Run `python check_angular_exact.py` from the parent directory.
The checker derives the Hopf metric and sphere projection, evaluates the
moments by an independent exact monomial formula, checks the antiderivative
and tail derivative symbolically, and verifies explicit rational-point
counterexamples to an all-angle sign claim. It uses no numerical integration.
