# A classical quadrupole generates a spatial mean

A single radial quadrupole of either sign generates a spatially uniform field component in the reconstructed hilltop potential. The sign of this mean reverses between the two archived potentials. At the common early time $x=0.12$, halving the initial perturbation cuts the sign-averaged mean by a factor within 0.00028 of four. This is a deterministic scalar-field result on a prescribed background. [Runs and convergence checks](NONLINEAR.json)

## Initial data and evolution

On a round spatial $S^3$, let $\chi$ be the angle from one pole and use the zonal harmonics $U_n(\cos\chi)=\sin[(n+1)\chi]/\sin\chi$, with their continuous values at the poles. Their eigenvalues under $-\Delta_\Omega$ are $n(n+2)$. They are orthonormal for the **spatial average** $\langle f\rangle=(2/\pi)\int_0^\pi \sin^2\chi\,f(\chi)\,d\chi$. In particular $U_0=1$, $\langle U_2\rangle=0$, $\langle U_2^2\rangle=1$, and $U_2$ ranges from $-1$ to $3$. Initialize $q(0,\chi)=A U_2$, $q_x(0,\chi)=0$ for both signs of $|A|=10^{-3},5\times10^{-4},2.5\times10^{-4},10^{-4},10^{-5}$. Each launch has initial RMS $|A|$, maximum $3|A|$, and exactly zero mean and mean velocity.

Use the prescribed stationary geometry $a=\cosh(hx)/h$, $H=a'/a$, and evolve the full scalar force for $q=\Theta-\Theta_*$. In dimensionless variables the equation is

$$
q''+3Hq'-a^{-2}\Delta_\Omega q=-3u'(\Theta_*+q),
\qquad u=V_E/V_{\rm top}.
$$

We reconstruct the prior spline potential and subtract its residual slope at the stationary root. The removed force has magnitude $3.34\times10^{-16}$ or less. The antiderivative used in the energy check includes the same correction. Projection onto degrees $0,2,\ldots,16$ respects reflection about the equator. A 96-node Gauss–Chebyshev quadrature of the second kind evaluates the full nonlinear force and 257 endpoint-inclusive angular nodes monitor extrema. The force mixes harmonics within this radial symmetry class. The clock ends when $\max_\chi |q|$ reaches 0.04 on the positive-minimum branch or 0.08 on the negative-minimum branch, or when an averaged kinetic, gradient or shifted-potential component reaches $0.1V_*$. Every nonzero launch reaches the field endpoint first. The separate stationary control remains exactly zero through $x=0.8$. A detuning rescales the units of proper time in these equations. The deterministic amplitudes here specify no quantum preparation or detuning dependence.

## Source of the uniform component

Define $\mu^2=-3u''(\Theta_*)>0$. If $q_L=A T_{2,qq}(x)U_2$ is the linear field, the second-order mean $m_2$ has zero initial data and satisfies

$$
m_2''+3H m_2'-\mu^2 m_2
  = -\tfrac32 u'''(\Theta_*)\langle q_L^2\rangle .
$$

The product identity $U_2^2=U_0+U_2+U_4$ explains why the spatial mean and degree four appear together. Evaluating the exact potential on both sides matters numerically because the inferred stationary position is within $2\times10^{-15}$ of a spline knot. A cubic spline is $C^2$ there, and its third derivatives can differ. Expanding the quotient that defines $u$ on each adjacent segment gives the two coefficients $-3u'''/2$, extrapolated from each side to the stationary position. They are (1080.74265, 1080.74426) on the positive branch and (−94.44709, −94.44548) on the negative branch. The comparison averages the appropriate sided coefficient weighted by $U_2^2$ for each sign of $A$. An additional forced solution drives the homogeneous linear response with $\langle-3[u'(\Theta_*+q_L)-u'(\Theta_*)]-\mu^2q_L\rangle$, retaining the full spline force on the linear profile. The nonlinear calculation itself uses the evolved angular field.

We also compute the retarded second-order mean independently at $X=0.12$. For the two analytic monopole transfer functions $R_0=T_{0,qq}$ and $S_0=T_{0,qv}$, its causal kernel is $G_0(X,s)=[S_0(X)R_0(s)-R_0(X)S_0(s)]\cosh^3(hs)$. Integrating $G_0 T_{2,qq}^2$ gives 0.01456939 and 0.01156827. Multiplying these by the appropriate force coefficient and $A^2$ reproduces the integrated second-order equation within $5.3\times10^{-11}$ relative error. [Analytic transfer derivation](../harmonics/HARMONICS.md)

| $\phi_m$ | Amplitude magnitude | Mean at $x=0.12$, $+A$ | Mean at $x=0.12$, $-A$ | Sign-averaged mean divided by $A^2$ |
|---:|---:|---:|---:|---:|
| +0.5 | $10^{-3}$ | $1.58339\times10^{-5}$ | $1.56606\times10^{-5}$ | 15.74724 |
| +0.5 | $10^{-5}$ | $1.57466\times10^{-9}$ | $1.57449\times10^{-9}$ | 15.74578 |
| −0.5 | $10^{-3}$ | $-1.10646\times10^{-6}$ | $-1.07876\times10^{-6}$ | −1.09261 |
| −0.5 | $10^{-5}$ | $-1.09272\times10^{-10}$ | $-1.09244\times10^{-10}$ | −1.09258 |

The sign average isolates the even-in-$A$ part. The full mean has a smaller odd part. On halving $|A|$ from $10^{-3}$ and from $5\times10^{-4}$, the respective even-mean ratios are 4.000279 and 4.000070 on the positive branch, 4.000081 and 4.000020 on the negative branch. For individual signed $10^{-5}$ seeds the second-order prediction reproduces the nonlinear mean to relative differences below 0.006% and 0.013%, respectively, at the probe time.

## Finite field endpoints and energy check

| $\phi_m$ | $A$ | First local 10% force deviation $x$ | Peak-field endpoint $x$ | Mean there | Degree-four coefficient there |
|---:|---:|---:|---:|---:|---:|
| +0.5 | $+10^{-3}$ | 0.17623 | 0.22905 | $+3.113\times10^{-4}$ | $+3.082\times10^{-4}$ |
| +0.5 | $-10^{-3}$ | 0.18648 | 0.24022 | $+3.988\times10^{-4}$ | $+3.730\times10^{-4}$ |
| −0.5 | $+10^{-3}$ | 0.35124 | 0.35278 | $-2.121\times10^{-4}$ | $-2.199\times10^{-4}$ |
| −0.5 | $-10^{-3}$ | Beyond endpoint | 0.34913 | $-1.607\times10^{-4}$ | $-1.346\times10^{-4}$ |

The local 10% event first occurs when any angular sample exits the interval bounded by the two directional force radii. The signed runs reach different endpoints because $U_2$ has unequal positive and negative extrema and the potential has unequal directional force radii. Degree 24 with 144 quadrature nodes and 513 monitoring nodes changes the $|A|=10^{-3},10^{-5}$ endpoints by less than $7\times10^{-13}$ in $x$ and the endpoint means by less than $7\times10^{-11}$ relatively. Degree 8 gives differences below $2.8\times10^{-7}$ in $x$ and $6.3\times10^{-6}$ in relative mean. The data record every retained mode at the bounce, probe and endpoint.

The averaged scalar energy shift divided by $V_*$ is
$$
e=h^{-2}\left[\tfrac16\langle q'^2+a^{-2}|\nabla_\Omega q|^2\rangle
 +\langle u(\Theta_*+q)-u(\Theta_*)-q\,u'(\Theta_*)\rangle\right].
$$
Its continuity equation gives $e'=-Hh^{-2}\langle q'^2+|\nabla_\Omega q|^2/(3a^2)\rangle$. Integrating this expansion loss separately and comparing its accumulated value with the energy difference at 61 times gives a maximum residual below $1.8\times10^{-12}$ of the largest absolute averaged energy component in each nonzero run. Across the launches the largest sampled averaged component stays below $0.0171V_*$. Local stresses can be appreciably larger. On the negative branch the largest endpoint kinetic and absolute potential corrections at individual angles reach $0.134V_*$ and $0.153V_*$ respectively.

This calculation follows a radial classical perturbation only up to the stated field and energy cutoffs. The result measures early generation of a homogeneous field component and higher harmonics. The fixed de Sitter metric cannot respond to their stress. A constraint-satisfying evolution of the geometry and a quantum state would be needed to turn this result into a backreaction calculation or a probability for the cosmological outcome.

## Reproduction

From the study directory, using NumPy and SciPy, run

~~~sh
OPENBLAS_NUM_THREADS=1 python -IB nonlinear/nonlinear.py --prior prior --output /tmp/NONLINEAR.json
~~~
