# Independent derivation of the prescribed quadratic material shell terms

15 October 2026 UTC. This note computes specified constitutive terms of the
retained coupled first variation, using the material outward-normal sum laws.
It preserves the earlier conditional construction of a complete finite-angular
second Taylor coefficient through contact. It does not compute that coefficient
quantitatively or give a finite-amplitude Einstein–scalar remainder.

## 1. Fields, convention and source sign

The background is the registered near-`phi=+1` regular branch, at
`delta=1/1000` and `c=96714321425692/161839258930731`. Use derivative
coefficients

\[
 U_\varepsilon=U_0+\varepsilon U_1+\frac{\varepsilon^2}{2}U_2,
 \qquad \psi_\varepsilon=\psi_0+\varepsilon\psi_1+
 \frac{\varepsilon^2}{2}\psi_2,
 \qquad\psi_1=S(t)Y,\quad Y=n_1^2-n_2^2.
\]

Here `psi` is the scalar pulled back to the common material shell. Increasing
the normal coordinate is outward on each copy. The exact scalar law is

\[
 \sum_{i=\pm}p_i+\sigma'(\psi)=0.
\]

Differentiating twice gives the material identity

\[
 \boxed{\sum_i p_{i,2}+\sigma''_b\psi_2
       =Q_\phi:=-\sigma'''_b\psi_1^2=-4S(t)^2Y^2.}       \tag{1}
\]

Thus the known Hessian term on the residual side is `+4 psi1^2`; the
inhomogeneous right side `Qphi` is negative. In the identical-copy sector,
dividing the entire equation by two gives `p2+(sigma''/2)psi2=-2 psi1^2`.
There is no additional `1/2` on either right side when `U2` denotes a
derivative coefficient. The normal, metric and embedding variations are
already included in the material derivatives `p_i,2`. They have not been
computed from the scalar trace alone.

The common tension and stress satisfy

\[
 \sigma_2=\sigma'_b\psi_2+T_{\rm q},\qquad
 T_{\rm q}:=\sigma''_b\psi_1^2,
\]
\[
 S_{2,ab}=-\sigma_bh_{2,ab}-\sigma'_b\psi_2h_{0,ab}
          -2\sigma'_b\psi_1h_{1,ab}-T_{\rm q}h_{0,ab}. \tag{2}
\]

`Tq` is the pure scalar constitutive contribution to the tension derivative.
The corresponding known stress term is `-Tq h0`. In the sum Israel equation
its pure constitutive right side is `+(Tq/3)h0`. These three conventions
must be distinguished before quoting a sign or norm. Equation (2) displays
the remaining terms; `Tq` is not the full second stress observable.

## 2. Unit S3 projections and exact power shares

Work on the unit round sphere in `R4`, with area `2 pi^2`. For ambient
coordinates `x`, set `s=sum x_i^2`, `A=x1^2+x2^2`, `B=x1^2-x2^2`, and

\[
 Z_2=A-s/2,\qquad H_4=B^2-sA/2+s^2/12.
\]

Direct Euclidean differentiation gives

\[
 \Delta B=\Delta Z_2=\Delta H_4=0,\qquad
 B^2=H_4+(s/2)Z_2+s^2/6.
\]

Restricting to `s=1` therefore gives the orthogonal decomposition

\[
 \boxed{Y^2=\frac16+\frac12Z_2+H_4,\qquad
 Z_2=n_1^2+n_2^2-\frac12,\quad
 H_4=Y^2-\frac12(n_1^2+n_2^2)+\frac1{12}.}             \tag{3}
\]

The sphere spectrum is `-l(l+2)`, hence `Delta Y=Delta Z2=-8` times the
same function and `Delta H4=-24 H4`. The ambient dimension is four; using
the `S2` spectrum changes the projections.

For normalized spherical mean, even moments are

\[
 \left\langle\prod_{j=1}^4 n_j^{2a_j}\right\rangle
 =\frac{\prod_j(1/2)_{a_j}}{(2)_{\sum_j a_j}}.
\]

This yields

\[
 \langle Y^2\rangle=\frac16,\quad
 \boxed{\langle Y^4\rangle=\frac3{40},\quad
 \langle Z_2^2\rangle=\frac1{12},\quad
 \langle H_4^2\rangle=\frac{19}{720}.}                 \tag{4}
\]

One independent check uses `J=n1^2+n2^2`, uniform on `[0,1]`, and an angle
uniform in the first two coordinates. Conditional on `J`,
`Y=J cos(2 theta)`; thus `mean Y4=(1/5)(3/8)=3/40`.
The functions `1,Z2,H4` are mutually orthogonal. Their unnormalized norms
are `||Y||^2=pi^2/3`, `||Z2||^2=pi^2/6`, and
`||H4||^2=19 pi^2/360`. In particular, `Y` has not been normalized to one.

At each time, the projections of (1) in these precise, unnormalized bases
are

\[
 Q_{\phi,0}=-\frac23 S^2,\qquad
 Q_{\phi,2}=-2S^2Z_2,\qquad
 Q_{\phi,4}=-4S^2H_4.                               \tag{5}
\]

The pure tension projections are `sigma'' S^2(1/6,Z2/2,H4)`.
The angular squared powers of `Y^2` in degrees zero, two and four are
proportional to `1/36,1/48,19/720`. Dividing by `3/40` gives

\[
 \boxed{(f_0,f_2,f_4)=\left(\frac{10}{27},\frac5{18},
                         \frac{19}{54}\right).}       \tag{6}
\]

These shares sum to one. They are shares of squared L2 power, not of signed
integrals, pointwise magnitude or response energy. Because every projection
has the same time factor, the same shares apply after integration with any
positive time weight. They also apply to `Qphi` and `Tq` on the retained
constant-background shell.

## 3. Proper-time norms

Keep `ell=10^-6`, `0<=t<=3 ell`, and `x=t/ell=tau/(rho_b ell)`. The retained
first-variation reconstruction gives

\[
 \ell S(\ell x)=G\cosh(\ell x)^{-3/2}H(x),\qquad
 G=\frac3{\rho_b^3p_b}<0,\quad
 K=\frac{3\pi^2}{\rho_b^2p_b^2}.
\]

All scalar norms in this note use the positive background material volume

\[
 dV_0=\rho_b^4\cosh^3t\,dt\,d\Omega.
\]

This is proper-time volume: `d tau=rho_b dt`, and the spatial volume factor
is `rho_b^3 cosh^3t`. Set

\[
 C=\int\psi_1^2\,dV_0=\frac K\ell\int_0^3H^2dx,
 \qquad J_4=\int_0^3\cosh(\ell x)^{-3}H(x)^4dx.        \tag{7}
\]

The quadratic cancellation of `cosh^3` is exact. In the quartic integral
one factor `cosh^-3` survives. Substitution and (4) give

\[
 \int\psi_1^4dV_0=\rho_b^4G^4\ell^{-3}\frac{3\pi^2}{20}J_4
               =\frac9{20}KG^2\ell^{-3}J_4.
\]

Consequently,

\[
 \boxed{\|Q_\phi\|_1=4C,\quad \int Q_\phi dV_0=-4C,\quad
 \|Q_\phi\|_2^2=\frac{36}{5}KG^2\ell^{-3}J_4,}        \tag{8}
\]
\[
 \boxed{\|T_{\rm q}\|_1=\sigma''_b C,\quad
 \|T_{\rm q}\|_2^2=\frac9{20}(\sigma''_b)^2KG^2
                              \ell^{-3}J_4.}          \tag{9}
\]

Here `sigma''_b=4+4 eta_b>0` follows from the registered tension. In a
different model, the L1 factor in (9) would be `|sigma''_b|`. Norms of the
scalar Israel amplitude `Tq/3` divide the L1 value by three and the squared
L2 value by nine. Norms of the per-copy scalar right side `Qphi/2` divide
them by two and four. No positive definite tensor norm or perturbed volume
has been selected; (8)–(9) are scalar amplitude norms on the background.

## 4. Continuous post-pulse signs and enclosures

Let `[h_-,h_+]` be the exact actual-H tail bracket in the retained report,
and let `[g_-,g_+]` be its exact negative gain bracket. Define the positive
rational hyperbolic-cosine enclosure `c_-<=cosh(t)^-3<=c_+` on `2 ell<=t<=3 ell`.
Then for every real `x` in `[2,3]`,

\[
 a_-:=\ell^{-2}g_+^2h_-^2c_-\ \le S(\ell x)^2\le\
 a_+:=\ell^{-2}g_-^2h_+^2c_+ .                         \tag{10}
\]

These are interval images of retained continuous enclosures, rather than
sampled signs. They are as sharp as these independent `H,G` brackets and
the specified hyperbolic-cosine enclosure permit; their endpoints are not asserted
to be attained jointly by the unknown exact background or waveform.
Their decimal displays are
`S^2 in [1.0309628629769591e-7,3.3831434489144685e-7]` and
`mean Qphi in [-2.255428965942979e-7,-6.873085753179728e-8]`.
Only the exact rational endpoints in the machine receipt decide inclusion;
these decimal displays are informational.

For a rational enclosure without a transcendental numerical decision, set
`P6(t)=1+t^2/2+t^4/24+t^6/720`. Positivity of the hyperbolic-cosine series and the
successive-tail ratio give

\[
 P_6(t)\le\cosh t\le P_6(t)+
 \frac{t^8}{8!\,[1-t^2/90]},\qquad 0\le t<\sqrt{90}.
\]

Using the lower series at `2 ell` and the upper series at `3 ell`, then
taking reciprocal cubes, supplies `c_-,c_+`. This avoids rounding either
endpoint to one. It improves the much coarser inherited bound
`(1-t^2)^2<=cosh^-3t<=1` by a negligible amount at this width.

Equations (5) and (10) give tail projection coefficient intervals

\[
 Q_{\phi,0}\in[-2a_+/3,-2a_-/3],\quad
 q_{\phi,2}\in[-2a_+,-2a_-],\quad
 q_{\phi,4}\in[-4a_+,-4a_-].                          \tag{11}
\]

The latter two intervals concern coefficients multiplying `Z2,H4`;
the harmonics change sign with angle. The full `Qphi` is nonpositive at
every angle, strictly negative where `n1^2!=n2^2`, and exactly zero on
`n1^2=n2^2`. Its normalized angular mean is strictly negative throughout
the tail. The pure tension mean lies in
`[sigma''_- a_-/6, sigma''_+ a_+/6]` and is strictly positive. The pure
stress amplitude `-Tq` has the opposite sign. The negative physical
quadrupole coefficient `S` does not change these quadratic signs.

## 5. Exact cubic constitutive algebra and its limit

The registered detuning is affine. With `eta=psi-1`,

\[
 \sigma(\eta)=\frac23+2\eta^2+\frac23\eta^3
              +\delta(1+c+c\eta).
\]

For any scalar increment `q`, there is therefore the exact identity

\[
 \boxed{\sigma(\psi_b+q)=\sigma_b+\sigma'_bq+
 \frac12\sigma''_bq^2+\frac23q^3,\qquad
 \sigma'(\psi_b+q)=\sigma'_b+\sigma''_bq+2q^2.}        \tag{12}
\]

On the prescribed scalar slice `q=epsilon psi1`, the tension remainder
after its quadratic Taylor polynomial is exactly
`R_sigma=(2/3)epsilon^3 psi1^3`. The bound

\[
 |R_\sigma|\le\frac23|\varepsilon|^3
             \|\psi_1\|_\infty\,\psi_1^2
\]

also yields the same Lp relative bound against
`(epsilon^2/2)Tq`, for `p=1,2`:

\[
 \|R_\sigma\|_p\le
 \frac{4|\varepsilon|\|\psi_1\|_\infty}{3\sigma''_b}
 \left\|\frac{\varepsilon^2}{2}T_{\rm q}\right\|_p.    \tag{13}
\]

A finite bound on `||psi1||_infinity` follows from all retained continuous
H cells, `|Y|<=1`, and `cosh^-3/2<=1`. Equation (13) is a constitutive
estimate on a prescribed scalar increment. It is not an amplitude interval
for an Einstein solution, a bound on the unknown `psi2`, or a nonlinear
defect-to-error estimate. If `q` itself includes unknown second/higher
response, (12) remains an algebraic identity but those fields remain
unknown. The scalar derivative law being exactly quadratic does not make
the metric, normal, clock or bulk Einstein equations quadratic.

## 6. Independent arithmetic and scope

The accompanying `check_quadratic_material_forcing.py` imports no retained
scientific program. It authenticates portable input bytes against the
recorded source hashes and manifest identities, proves the polynomial,
moment, factorial, copy and norm algebra with exact SymPy expressions,
and independently integrates squared/fourth powers of all 4353 retained
closed-cell H intervals using `Fraction`. It checks gapless `[0,3]`
coverage and all 256 strictly positive tail cells. The fresh exact H2
consequence is checked inside the earlier response certificate. The
weighted H4 bound uses the explicit hyperbolic-cosine enclosure above.
Its decimal display is `J4 in [86.46022616035648,87.42628460767193]`, giving
`||Qphi||_2^2 in [5.778773178542654e21,8.643174843435045e21]`.

The optional `--check-computed REPORT CELLS` mode independently reconstructs
every new projection coefficient, checks the complete new CSV byte sequence,
and compares all twelve numerical report paths for moments, norms,
parameters and tail bounds. It reproduces the producer's coarser inherited
clock bound separately from the sharper clock calculation; the latter's
intervals must lie inside the former's. It imports neither producer.

`SOURCE_PINS.json` identifies every prior scientific source read during
this derivation, including the planning note and earlier complete boundary
system. The packaged runner needs only the copied minimal input files and
their input pins. Its normal and optimized receipts record hashes and exact
rational results. Mathematical corruption modes alter the outward-sum
factor, factorial, sign, sphere dimension, normalization, fourth moment,
harmonic decomposition, proper-time Jacobian, quartic cosine weight and
cubic remainder coefficient; each must fail with an explicit exception.
Four further controls corrupt the expected new cell coefficient, fourth
moment, squared norm and tail sign before independent comparison. A simple
polynomial counterexample verifies that the sign of a forcing alone cannot
fix a response sign: `1-2t^2` and `-1-2t^2` have the same negative second
derivative and opposite signs at zero. It is not proposed as the shell
operator or as a theorem about a particular solved shell response.

The unknown full `psi2`, first and second metric fields, normal/embedding
jets, coupled bulk sources, second Cauchy correction and material proper-time
map still enter the complete second-order solution and stress observable.
Their existence and finite-angular construction retain the earlier
conditional result and hypotheses. A forcing sign does not determine the
sign of that response. These exact local source computations supply no
finite-amplitude waveform, effective Einstein-family epsilon radius,
observed time calibration, gravitational energy, cosmological-origin
inference, external peer review or novelty claim.
