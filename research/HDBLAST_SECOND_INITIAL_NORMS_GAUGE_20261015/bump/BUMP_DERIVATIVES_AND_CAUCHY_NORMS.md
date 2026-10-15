# Exact bump derivatives and the second Cauchy data

This standalone calculation supplies global derivative bounds for the exact
retained bump. It changes no frozen checkpoint. No sampled maximum is used.

Let

\[
 Q(r)=\exp\!\left(4-\frac1{(r-1)(2-r)}\right),\quad1<r<2,
 \qquad Q=0\text{ elsewhere},\qquad q_\ell(z)=Q(z/\ell),
 \quad\ell=10^{-6}.
\]

The following rational bounds hold on the whole real line, including the
smooth zero extension:

| j | Bound for `sup |Q^(j)|` |
| --- | ---: |
| 0 | 1 |
| 1 | 16 |
| 2 | 256 |
| 3 | 46656/7 |
| 4 | 16777216/49 |
| 5 | 10000000000/343 |
| 6 | 8916100448256/2401 |

In particular `sup |Q''''| < 342393`. These are conservative envelopes,
not claimed optimal constants. The physical radial derivative bound is
`sup |d_z^j q_ell| <= ell^-j B_j`, where `B_j` denotes the table entry.
The physical width factors must be retained before propagating these bounds.

## Recurrence and all-point proof

Put `h=(r-1)(2-r)`, `s=1/h>=4`, and `v=h'=3-2r`. Then

\[
 s'=-vs^2,\qquad v'=-2,\qquad v^2=1-4/s,\qquad |v|\le1.
\]

There is a polynomial `P_n(s)` such that

\[
 Q^{(n)}=e^{4-s}\begin{cases}P_n(s),&n\text{ even},\\
 vP_n(s),&n\text{ odd}.
 \end{cases}
\]

Starting from `P0=1`, differentiation gives the exact recurrence

\[
 P_{n+1}=s^2(P_n-P_n')\quad(n\text{ even}),
\]
\[
 P_{n+1}=(s^2-4s)(P_n-P_n')-2P_n\quad(n\text{ odd}).
\]

The extra `-2P_n` is the derivative of `v`; replacing `v^2` by one would
also give a wrong even derivative. The first higher polynomials are

\[
 P_3=s^6-10s^5+24s^4-12s^3,
\]
\[
 P_4=s^8-20s^7+136s^6-384s^5+420s^4-120s^3,
\]
\[
 P_5=s^{10}-28s^9+276s^8-1200s^7+2340s^6-1800s^5+360s^4,
\]
\[
 P_6=s^{12}-42s^{11}+678s^{10}-5464s^9+23820s^8
            -56400s^7+68040s^6-35280s^5+5040s^4.
\]

For every `n<=6`, both polynomials `s^(2n)+P_n(s)` and
`s^(2n)-P_n(s)` have nonnegative coefficients after substituting `s=u+4`.
The checker performs these exact coefficient calculations. Therefore

\[
 |P_n(s)|\le s^{2n}\quad(s\ge4),\qquad0\le n\le6.
\]

For `d>=0`, the derivative of `s^d e^(4-s)` has the sign of `d-s`, so
its maximum on `s>=4` is at `max(4,d)`. Further,

\[
 e^2>\sum_{j=0}^4\frac{2^j}{j!}=7.
\]

For `n>=2` this gives

\[
 \|Q^{(n)}\|_\infty
 \le (2n)^{2n}e^{-2(n-2)}
 \le \frac{(2n)^{2n}}{7^{n-2}}.
\]

The inequality is strict for `n>2`. Orders zero and one give the displayed
constants directly. This proves the table for all points in the open bump.

The recurrence works at arbitrary finite order and always produces a
polynomial of degree `2n`. Thus every derivative tends to zero at either
endpoint, because `s^m e^-s -> 0` as `s -> infinity` for each fixed m.
For example the positive exponential series gives
`s^m e^-s <= (m+1)!/s`. The zero extension is therefore smooth and flat
to every order, not only through order six.

The verifier also offers a scalable polynomial-envelope calculation. At any
finite order, expand `P_n(u+4)=sum a_j u^j` and choose an integer

\[
 C_n\ge\max_j\frac{|a_j|}{\binom{2n}{j}4^{2n-j}}.
\]

Then `C_n s^(2n) +/- P_n(s)` have nonnegative shifted coefficients, giving
`|P_n|<=C_n s^(2n)`. The algorithm finds `C_n=1` for every order in the
table. Higher orders may use larger multipliers. The same exponential
argument supplies a rational bound without locating numerical extrema.

## What this supplies for the physical Cauchy norm

The second initial constraint construction uses physical variables
`(gamma2,K2,phi2,Pi2)`, with K a covariant tensor and Pi the unit-normal
scalar velocity. It is distinct from an affine spacetime path and from
canonical momentum-density variables. These bump bounds do not change those
conventions or the previously derived lapse conversion.

For the actual initial longitudinal amplitudes,

\[
 u=Cq_\ell,\quad u_t=Ce^{-z}q_\ell,\quad
 v=B(q_{\ell,z}+kq_\ell),\quad
 v_t=Be^{-z}(q_{\ell,z}+(k-1)q_\ell).
\]

Here `C,B,k` are smooth background coefficients and
`d_y=-rho^-1 d_z`. Product bounds can be propagated exactly by Leibniz:

\[
 \|\partial_z^j(Aq_\ell)\|\le
 \sum_{m=0}^j\binom jm A_{j-m}\ell^{-m}B_m,
\]

where `A_j` bounds the corresponding background derivative. Similarly

\[
 \|\partial_z^j(e^{-z}q_\ell)\|\le
 \sum_{m=0}^j\binom jm\ell^{-m}B_m
\]

on `z>=0`. Each radial y derivative must differentiate the `rho^-1`
coefficient when there is more than one derivative. Multiplying the jth
z-derivative bound by `rho_min^-j` alone is not generally valid.

The required bump orders are:

| Second physical datum | Value | One spatial y derivative |
| --- | ---: | ---: |
| gamma2, from Hamiltonian amplitudes f_l | 2 | 3 |
| K2=(F/3)Y^2 gamma0 | 1 | 2 |
| phi2=sum -3f_l,y W_l/(2p) | 3 | 4 |
| Pi2=(G-F_y)Y^2/p | 2 | 3 |

Thus fourth order suffices for the values and one spatial derivative of
all four second physical Cauchy fields. The highest requirement comes from
`phi2_y`, which differentiates the quadratic Hamiltonian source twice.
That source already contains `u_yy` and `v_y^2`, so its second derivative
contains the fourth bump derivative. No fifth or sixth bump derivative is
required for this first-order target; those bounds are reusable additions.

For an adapted-frame covariant C1 norm, background connection terms and the
finite angular harmonic derivative bounds must also be included. They do
not raise the radial bump order. A Sobolev or integrated energy norm additionally
needs the specified norm, volume measure and angular normalization. This
bump calculation alone does not numerically enclose any such full physical
Cauchy norm: the background coefficient jets and all source/product
correlations must be propagated in the separate norm checker.

The corrections remain compact and flat at `z=ell,2ell`, because the
constraint construction uses finite differentiation and multiplication by
smooth coefficients with `p` bounded away from zero. A large but finite
derivative envelope does not change that support statement.

## Reproduction and check scope

```sh
python -I -B check_bump_derivatives.py --output BUMP_DERIVATIVE_BOUNDS.json
python -I -O -B check_bump_derivatives.py --output BUMP_DERIVATIVE_BOUNDS_OPTIMIZED.json
python -I -B check_bump_derivatives.py --run-controls --output BUMP_DERIVATIVE_CONTROLS.json
```

The standalone program imports only the Python standard library and needs no
external input. Its portable report exposes `ell` and
`dimensionless_Q_derivative_sup`, an array of exact rational strings indexed
by derivative order. It also supplies `physical_q_derivative_sup` for the
same indices. Use `--order N` to extend the recurrence beyond the default six.

An independent formal-composition check expands
`h(r+d)=h(r)+(3-2r)d-d^2`, then computes the reciprocal and exponential series
directly. It compares exact derivative-to-Q ratios with the recurrence at
five rational interior points, giving 35 witnesses at the default order.
Those finite witnesses check the identities; the supremum certificate is the
all-point polynomial and exponential argument above.

Six deliberate errors must each fail through an explicit exception under
normal and optimized Python: omit `v'`, replace `v^2` by one, flip the sign
of `s'`, maximize the eighth-degree exponential envelope at four, assert a
false fourth-derivative bound, or omit one physical ell factor. No assertion
statement or floating-point calculation determines acceptance.
