# A shared source history determines the initial cross-covariance at fixed mode variances

Suppose two homogeneous, isotropic preparations of the closed finite-band scalar share the same initial mode variances $X(q)=\langle\chi_q^2\rangle$ and $P(q)=\langle\Pi_q^2\rangle$ at each momentum. Their symmetrized cross-covariances $C(q)=\langle\{\chi_q,\Pi_q\}/2\rangle$ may differ. On a common smooth shell history their scalar-force histories can agree throughout an interval starting at the preparation only if their cross-covariances also agree almost everywhere in the populated band. On a grid of $N$ distinct positive momentum magnitudes, $N$ odd derivatives of the reduced variance suffice to identify all $N$ covariance differences in exact arithmetic.

Use proper time $s$, conformal time $d\eta=ds/a$, canonical momentum $\Pi=a^3d\chi/ds$ and rescaled field $v=a\chi$. For every real scalar oscillator with squared comoving momentum $\lambda=q^2$, the equation on the shared history is $v''+(\lambda+U)v=0$, where $U=a^2m^2-a''/a$ and primes denote conformal-time derivatives. For the difference $D=\Delta\langle v^2\rangle$ the product equation and the common initial variances give

$$
D'''+4(\lambda+U)D'+2U'D=0,\quad
D_i=0,\quad D'_i=2\Delta C(q),\quad D''_i=4(a'_i/a_i)\Delta C(q).
$$

These identities follow from $v'=(a'/a)v+\Pi/a$ and differentiating the covariance. The scalar force difference is $\Delta j=(\partial_\phi m^2)/(2a^2)\int D\,d\nu$, with $d\nu=q^2dq/(2\pi^2)$. The common instantaneous vacuum subtraction cancels. Assume the factor multiplying the integral is nonzero on an interval immediately after preparation. Source equality there implies that every right derivative of the integrated $D$ at preparation vanishes when the history is smooth on that side.

Differentiating the product equation recursively gives $D^{(2r+1)}_i=P_r(\lambda)\Delta C(q)$ with $P_r$ of degree $r$ and leading coefficient $2(-4)^r$. Indeed, the term $-4\lambda D^{(2r-1)}$ provides the next leading power. The derivatives of $U$ contribute lower powers. The even derivative of order $2r+2$ has degree at most $r$. The resulting triangular system determines $\int \lambda^r\Delta C\,d\nu$ successively from the odd derivatives. A finite momentum band and finite integrable covariance permit differentiation under the integral at each fixed order.

If the integrals of all these monomials vanish, the signed measure $\Delta C\,d\nu$ vanishes. Polynomials in $q^2$ are dense in the continuous functions on $[0,Q]$ since $q\mapsto q^2$ is one-to-one there. For a positive-weight finite quadrature the equivalent statement follows by polynomial interpolation at the distinct nodes. The continuum conclusion follows from smoothness at the fork and compact support alone, by equality of the integrals against every polynomial.

The relation also gives examples that hide a source difference for several derivatives. In units $x=q/Q$, take $\Delta C$ proportional to $S(x)=1-(14/3)x^2+(21/5)x^4$. Its first two weighted moments $\int_0^1x^2S(x)dx$ and $\int_0^1x^4S(x)dx$ vanish, while $\int_0^1x^6S(x)dx$ is positive. The corresponding reduced force first differs at derivative five, with an overall factor set by $Q$ and the covariance scale. These small differences can be realized in an interior covariance disk with sufficient initial uncertainty. At a pure fork with fixed variances admissibility must also be checked against the nodewise uncertainty bound.

This is a uniqueness statement about the closed Gaussian scalar on a hypothetical common geometry. The amplitude of the separation can approach zero, and high-order derivative recovery can be ill conditioned at finite measurement precision. The uniform canonical shear in the native experiment changes $P(q)$ and lies outside this fixed-variance statement. [Exact polynomial and integration check](force_jet.py) · [Output](FORCE_JET.json)
