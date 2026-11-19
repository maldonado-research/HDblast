# An independent test of the frozen parent response

Two measurements of the same oscillator system test the sign of the predicted response. The first solves the parent and the bath together under an applied force. The second eliminates the bath analytically. Matching the static restoring force removes the zero-frequency shift. A small residual change in the slope then measures the additional kinetic coefficient needed to match the frequency dependence. The calculation also gives a way to integrate the remaining frequency curvature accurately when a mode lies close to the spectral threshold. [Numerical results](HAMILTONIAN.json) · [Runnable calculation](hamiltonian.py)

## Proper-time force and normalization

Freeze the scale factor $a$, gate $F$, parent proper momentum $p=q/a$ and mass $m$. Add an external source $-a^3 J\chi$ to the proper-time Hamiltonian, so $J$ is the force per unit parent inertia at zero coupling. Let $b=yF$ and $Z$ be the kinetic coefficient multiplying $a^3\dot\chi^2/2$ in the Lagrangian. Its Hamiltonian term is $\Pi^2/(2a^3 Z)$. Use the bath strengths of the preceding energy study with $C=\sum_jg_j^2/\omega_j^2$ and $D=\sum_jg_j^2/\omega_j^4$. Write $\alpha$ for the multiplier of the static potential term $a b^2 C\chi^2/2$ and $\gamma$ for the multiplier in $Z=1-\gamma b^2D$. A positive $Z$ is required in this realization. [Original finite Hamiltonian](https://github.com/maldonado-research/HDblast/blob/2cbafaefde87652135d8222e50a53da4eb704c41/research/HDBLAST_CAUSAL_YUKAWA_ENERGY_20261119/theory/THEORY.md)

For a force proportional to $\exp(-iws)$, the causal elimination of each oscillator gives $x_j=-ba g_j\chi/(\omega_j^2-a^2w^2)$ in the upper half $w$ plane. The inverse proper-time susceptibility is

$$
\frac{J}{\chi}=p^2+m^2-Zw^2+
\frac{b^2}{a^2}\left[\alpha C-
\sum_j\frac{g_j^2}{\omega_j^2-a^2w^2}\right].
$$

The factors of $a$ follow independently from $\dot x_j=\pi_j/a$ and the source normalization. In terms of $x=(\omega/a)^2$, the continuous proper kernel at Euclidean squared frequency $u=\nu^2$ has measure $(x-p^2)\,dx/(8\pi^2)$. For a fixed proper cutoff $\Omega=\Lambda/a$, changing $a$ while scaling $q,\Lambda$ together leaves the proper response invariant. In an expanding spacetime, holding the comoving cutoff fixed would change its proper value. This check assumes a frozen background.

For the matched case $\alpha=\gamma=1$, the response at $w=i\nu$ satisfies

$$
\frac{J}{\chi}=p^2+m^2+u-b^2R_\Omega(u),\qquad
R_\Omega(u)=\frac{u^2}{8\pi^2}\int_{p^2}^{\Omega^2}\frac{x-p^2}{x^2(x+u)}\,dx.
$$

The positive curvature of the spectral kernel therefore **lowers** the Euclidean inverse response after the two matches. Both the intercept and its slope at $u=0$ equal their free values. These are statements about this quadratic oscillator model for a frozen mode with $p>0$.

## Full canonical evolution

The independent calculation constructs the first-order equations of the parent and 32 bath oscillators from the quadratic Hamiltonian matrix. It solves the driven 66-dimensional system at four imaginary frequencies and a complex probe $w=2.6+0.4i$. It uses $p=1.5$, $m=2.3$, $\Omega=11$, $y=0.4$ and $F=0.8$. The change of $a$ from 1 to 1.7 shifts the coordinate frequencies and canonical momenta but preserves the proper susceptibility to $2.8\times10^{-16}$. The eliminated and full inverse responses agree to $4\times10^{-15}$ over the test.

| Frozen model | Shift of inverse response at zero | Slope with respect to $\nu^2$ at zero |
|---|---:|---:|
| Bare | $-0.14238018$ | 1.00389521 |
| Static compensation | 0 | 1.00389521 |
| Static and kinetic compensation | 0 | 1.00000000 |

As a separate causal check, a 10-oscillator bath starts at rest and is driven in proper time by $J(s)=\exp[\nu(s-s_f)]$ over 32 exponential growth times. With $J(s_f)=1$, the endpoint $\chi(s_f)$ measures the response at $w=i\nu$, up to the transient from the finite start. Direct integration at $\nu=1.15p$ agrees with the three predicted responses to $2.5\times10^{-15}$ in the reported run. All finite-bath results apply to the oscillator Hamiltonian used for the test. Its spectral calibration reproduces the leading fermionic kernel. Extrapolating this finite-coupling control to higher fermionic correlation functions would require additional input.

## The part above the cutoff

Write $X=\Omega^2$. Substituting $t=X/x$ above the cutoff gives an exact expression for the omitted tail,

$$
R_\infty-R_\Omega=\frac{u^2}{8\pi^2 X}
\int_0^1\frac{1-p^2t/X}{1+ut/X}\,dt.
$$

Dividing by the prefactor gives a positive ratio $T$. Expanding the denominator once and integrating the remainder yields

$$
0\le T-\left[1-\frac{p^2+u}{2X}\right]\le\frac{u(p^2+u)}{3X^2}.
$$

This gives an absolute error bound when using the first two tail terms at finite cutoff. Integrating the finite and omitted ranges independently recovers $R_\infty=[(p^2+u)\log(1+u/p^2)-u]/(8\pi^2)$. Across three momenta, four frequency ratios and four cutoffs, the largest missing fraction was 6.70%, at $p=7.98075$, $\nu/p=2$ and $\Omega=60$. That fraction describes a truncation of the continuous spectrum. Errors from resolving the finite integral are evaluated separately below.

## A small quadrature that resolves the threshold

Set $t=p^2/x=(q/\omega)^2$, $k=u/p^2$ and $\beta=p^2/\Omega^2$. The finite remainder becomes

$$
R_\Omega=\frac{u^2}{8\pi^2p^2}\int_\beta^1\frac{1-t}{1+kt}\,dt .
$$

The momentum and cutoff enter through the lower endpoint. The nearest pole stays outside the integration interval for $0<\nu/p\le2$. A small, fixed set of Gauss–Legendre nodes can therefore resolve the continuous curvature across low momenta and large cutoffs. The benchmark compares directly to a 52-digit antiderivative for all 384 combinations of the 24 native momenta, four cutoffs and four probe ratios.

| Gauss nodes | Largest measured relative error | Uniform relative bound in exact arithmetic |
|---:|---:|---:|
| 8 | $5.86\times10^{-7}$ | $2.15\times10^{-6}$ |
| 12 | $2.69\times10^{-10}$ | $9.75\times10^{-10}$ |
| 16 | $1.24\times10^{-13}$ | $4.42\times10^{-13}$ |
| 24 | $1.34\times10^{-15}$ | $9.07\times10^{-20}$ |

The last row reaches floating-point precision. For the bound, map the interval to $v\in[-1,1]$ and write $1+kt=A+Bv$ with $A/B\ge 3/2$. Let $n$ be the quadrature order and $d=2n-1$. The Chebyshev residual $e(v)=T_d(v)/T_d(-A/B)$ equals one at the pole, so $(1-e(v))/(A+Bv)$ is a polynomial of degree $d-1$. After multiplication by the linear numerator, the resulting approximation $P$ has degree at most $d$ and is integrated exactly by the Gauss rule. For the positive integrand $f$ the pointwise bound is $|f-P|\le\delta f$ with $\delta=1/\cosh[d\,\operatorname{arcosh}(3/2)]$. Positivity of the weights gives $|Q-I|\le\delta(Q+I)$ for the numerical and exact integrals $Q,I$. Hence $|Q/I-1|\le2\delta/(1-\delta)$ for any $0<p<\Omega$ in the stated frequency range.

The 24-node result separately checks every entry of the released response tables. Across the 384 probes, the largest relative disagreements with the analytic finite remainder, infinite limit and missing tail are below $6.1\times10^{-15}$. The 2,048-node linear-frequency bath differs by at most $8.51\times10^{-10}$ from the compact integral. [Archived tables](../response/TABLES.npz)

This coordinate substitution targets the continuous, twice-subtracted frozen response. Time evolution of an existing discretized bath retains the frequencies and weights specified in its Hamiltonian. A driving shell that varies in time would also require the stress and scalar source associated with its matched dynamics.

## Replay

From the study root, run

~~~sh
OPENBLAS_NUM_THREADS=1 python -IB audit/hamiltonian.py --native ../HDBLAST_COVARIANCE_FEEDBACK_20261119/native --output /tmp/HAMILTONIAN.json
~~~

The program compares the mapped integral to the sibling `response/TABLES.npz` by default. Pass `--tables` to check a regenerated table. The optional native input supplies the fine-grid parent momenta. With it omitted, the program reconstructs the 24-point Gauss–Legendre parent grid on $(0,8)$.
