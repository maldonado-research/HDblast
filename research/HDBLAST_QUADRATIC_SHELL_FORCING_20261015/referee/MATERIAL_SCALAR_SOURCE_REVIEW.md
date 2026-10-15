# Independent physical and angular audit of the quadratic scalar-junction source

The proposed calculation is legitimate with one important scope distinction:
it computes the explicit tension nonlinearity in the second differentiated
**material scalar junction**. It does not by itself compute the full quadratic
boundary forcing in arbitrary bulk harmonic coordinates, a second-order
solution, or a second-order shell stress.

All statements below concern coefficients of a conditional perturbative
construction. No finite-amplitude nonlinear solution, quantitative remainder,
novelty claim or physical time calibration follows.

## 1. Exact material source and the factor of two

Use the same material shell, scalar trace `psi`, and outward normal derivatives
`p_+,p_-` as in the retained physical junction construction. Its exact scalar
law is

\[
P= p_++p_-+\sigma'(\psi)=0.
\]

The registered tension, with `eta=phi-1`, is

\[
\sigma=2(1/3+\eta^2+\eta^3/3)
+\delta[1+c(1+\eta)],\qquad
\sigma''=4+4\eta=4\phi,\quad\sigma'''=4.
\]

For the derivative convention
`U_epsilon=U0+epsilon U1+epsilon^2 U2/2`, twice differentiating gives

\[
\boxed{p_{+,2}+p_{-,2}+\sigma''_b\psi_2
=Q_\phi,\qquad Q_\phi=-4\psi_1^2.}
\tag{1}
\]

For identical copies the row divides by two:

\[
p_2+2\phi_b\psi_2=-2\psi_1^2.
\tag{2}
\]

The second derivative source in (1) is `-4 psi1^2`; the corresponding
coefficient of the raw `epsilon^2` term is `-2 psi1^2`. Neither the per-copy
division nor this Taylor-factor distinction can be silently changed.
Also `sigma''_b` is `4 phi_b`, rather than exactly four.

## 2. Why the certified first observable may be used

The preceding physical first variation is

\[
\psi_1(\tau,\Omega)=S(\tau/\rho_b)Y(\Omega),
\quad Y=n_1^2-n_2^2,\quad t=\tau/\rho_b.
\tag{3}
\]

The background material shell scalar is constant. Consequently a first-order
change of intrinsic shell clock or angular labels changes the first scalar
trace by a Lie derivative of a constant, hence by zero. The same first
variation may therefore be used in (1), provided its embedding and harmonic
gauge reconstruction are the retained physical ones. A bare bulk Eulerian
scalar perturbation cannot replace that physical trace without its pullback.

The invariance of the source does not eliminate second-order clock terms.
For a family-dependent intrinsic relabeling with first generator `xi1`,

\[
\widetilde\psi_2=\psi_2+2\mathcal L_{\xi_1}\psi_1,
\quad\widetilde p_{i,2}=p_{i,2}+2\mathcal L_{\xi_1}p_{i,1},
\]

because the background traces are constant. The combination on the left
of (1) changes by `2 L_xi1(P1)=0` on a solution of the first scalar junction.
Thus the complete material row is consistent under the clock/label change;
its individual second jets still depend on that change.

For an ambient collar pullback, the second material derivative of a tensor is

\[
Q^{\mathrm{mat},2}=Q^{(2)}+2\mathcal L_{\xi_1}Q^{(1)}
+(\mathcal L_{\xi_2}+\mathcal L_{\xi_1}^2)Q_0.
\tag{4}
\]

Here `xi2` is the parameter derivative of the generator. Background normal
derivatives, normal normalization, embedding evaluation, tilt and metric
variations are included when (4) constructs the physical material coefficients.
They do not vanish because the intrinsic background trace happens to be constant.

## 3. Coordinate boundary sources that are not determined by (1) alone

The prior boundary-fixed wave chart writes `g_yy=N^2`. With first lapse
`nu1`, first wave-gauge scalar `chi^w`, second scalar `w`, and second lapse
`nu2`, the per-copy scalar row is

\[
w_y+\frac{\sigma''_b}{2}w+\frac{\sigma'_b}{2}\nu_2
=-2(\chi^w)^2-\sigma''_b\chi^w\nu_1.
\tag{5}
\]

The mixed lapse term is indispensable. In the first longitudinal representative,
`g1_yy=-4uY`, `h1_ab=2uY h0_ab`, `chi=vY`. Thus
`nu1=-2uY` and `nu2=k_yy/2-4u^2Y^2`. Substitution into the same physical
row gives a scalar source

\[
\left(-4p_bu^2+2\sigma''_buv-2v^2\right)Y^2,
\qquad p_b=-\sigma'_b/2,
\tag{6}
\]

after the linear `k_yy` term is put on the left. The corresponding extrinsic
curvature row has source
`(sigma''_b v^2/6-12 a_b u^2)Y^2 h0`, in its stated normalization.
Wave-gauge first fields must be obtained from the actual first gauge vector;
the longitudinal and wave-gauge traces cannot be identified by notation.

To build a complete harmonic-coordinate second-order system one must also
retain the metric junction source, the reduced quadratic bulk and gauge
sources, the material/Eulerian conversion, second embedding and normal jets,
proper-clock pullback, induced matching, corner compatibility, constraints
and chosen homogeneous second-order initial data. Formula (1) specifies the
explicit nonlinear tension term in physical scalar-jet variables; it is not
the entire right-hand side after those variables are expanded in bulk fields.

## 4. Exact angular mean, decomposition and full angular norm

Use the unit `S3` of volume `2 pi^2`, normalized angular average `mean`, and
the ambient polynomials

\[
s=\sum_{i=1}^4x_i^2,\quad A=x_1^2+x_2^2,\quad B=x_1^2-x_2^2,
\quad Z_2=A-s/2,
\quad\mathcal H_4=B^2-sA/2+s^2/12.
\]

`B,Z2,H4` are homogeneous harmonic polynomials of degrees two, two and four.
Their restrictions satisfy

\[
Y^2=\mathcal H_4+\frac12Z_2+\frac16,
\quad
Q_\phi=-4S^2\mathcal H_4-2S^2Z_2-\frac23S^2.
\tag{7}
\]

Thus the known scalar source has precisely these degree-zero, degree-two and
degree-four components when `S` is nonzero. This statement supplies no vector
or tensor source projections for the full system.

The normalized moments and orthogonality are

\[
\operatorname{mean}Y^2=\frac16,\quad
\operatorname{mean}Y^4=\frac3{40},\quad
\operatorname{mean}Z_2^2=\frac1{12},\quad
\operatorname{mean}\mathcal H_4^2=\frac{19}{720}.
\]

All distinct displayed degree components are orthogonal. Therefore

\[
\boxed{\operatorname{mean}Q_\phi=-\frac23S^2,\qquad
\int_{S^3}Q_\phi^2\,d\Omega=\frac{12\pi^2}{5}S^4.}
\tag{8}
\]

The degree `0,2,4` contributions to normalized `mean(Q_phi^2)` are respectively
`4S^4/9`, `S^4/3`, `19S^4/45`; they sum to `6S^4/5`.
For a physical spatial section of radius `R= rho_b cosh(t)`, multiply the
squared angular norm in (8) by `R^3`. The per-copy source has half the mean
and half the norm, hence one quarter of the squared norm.

## 5. Proper-shell spacetime identities

Use the same background observation window and positive shell measure as the
previous response coefficient,

\[
d\mu_0=\rho_b^4\cosh^3t\,dt\,d\Omega,\qquad
\mathcal C_{3\ell}=\int\psi_1^2\,d\mu_0.
\]

Since `Q_phi=-4 psi1^2<=0` pointwise,

\[
\boxed{\int Q_\phi\,d\mu_0=-4\mathcal C_{3\ell},
\qquad \|Q_\phi\|_{L^1(d\mu_0)}=4\mathcal C_{3\ell}.}
\tag{9}
\]

The already certified **same wider state** therefore gives
`-2.976e12 < integral Q_phi dmu0 < -2.424e12`. This uses a previously
certified norm by an exact identity, rather than transferring a coefficient
between different widths. It is a signed scalar-source integral, not energy
or an isotropic gravitational force.

With the retained physical gain `G=3/(rho_b^3 p_b)`, normalized trace
`H(x)=ell D(ell x)`, and `K=3pi^2/(rho_b^2p_b^2)`, the full spacetime norm is

\[
\boxed{\|Q_\phi\|_{L^2(d\mu_0)}^2
=\frac{36\mathcal K\mathcal G^2}{5\ell^3}
\int_0^3\cosh^{-3}(\ell x)H(x)^4\,dx.}
\tag{10}
\]

The fourth moment is a new integral and is not fixed by `C`. Unlike the
quadratic response weight, its `cosh^-3` factor remains. Alternatively one
may integrate `(12pi^2/5)rho_b^4 ell^-3 cosh^3(ell x)[ell S(ell x)]^4`.
These two formulas must agree. All measures and clocks here are background
coefficient measures; no exact finite-amplitude family measure is being inferred.

The companion explicit constitutive term is
`T_quad=sigma''_b psi1^2`. Its angular coefficient multipliers are
`sigma''_b S^2 (1/6,1/2,1)`, its signed integral is `sigma''_b C`, and its
squared spacetime norm is `(sigma''_b)^2/16` times (10). It is positive in
the actual retained shell interval. However the full second derivative of
the tension is `sigma'_b psi2+T_quad`, so neither its sign nor the sign of
the full second stress is determined by this positive constitutive part.

## 6. What the sign does not prove

`Q_phi` is nonpositive and is zero on the nodal set of `Y` and whenever `S=0`.
Its mean is strictly negative when `S` is nonzero. This is not strict
negativity at every angular point.

The sign does not determine `psi2`, the complete `U2`, or the second shell
stress. In particular

\[
S_{2,ab}=-\sigma_bh_{2,ab}-\sigma'_b\psi_2h_{0,ab}
-2\sigma'_b\psi_1h_{1,ab}-\sigma''_b\psi_1^2h_{0,ab}
\]

contains second fields and mixed metric terms beyond the scalar source.
The coupled boundary evolution need not have a positivity-preserving inverse,
and free homogeneous second-order data remain to be fixed. The earlier
conditional complete second-coefficient construction is preserved; this
calculation supplies one quantitatively prescribed material source component.
